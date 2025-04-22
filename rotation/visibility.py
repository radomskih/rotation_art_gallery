import math
from shapely.ops import polygonize, unary_union, snap
from shapely.geometry import Polygon, LineString, Point

#Visibility Functions:

def valid_intersections(ray, polygon, origin): #find contributions of ray to origin's visibility polygon
    tolerance = 1e-6
    snapped_ray = snap(ray, polygon.boundary, tolerance) #make sure ray connects with one of the vertices
    intersections = snapped_ray.intersection(polygon.boundary) #collect all intersections with polygon

    if intersections.is_empty:
        return None
    elif isinstance(intersections, Point):   
        intersections = [intersections]
    elif hasattr(intersections, 'geoms'):  
        intersections = [g for g in intersections.geoms if isinstance(g, Point)]
    else:
        intersections = [intersections]

    def far_enough(p):
        return isinstance(p, Point) and origin.distance(p) > 1e-8
    
    filtered = [p for p in intersections if far_enough(p)] #filter out weirdness with edges incident to origin
    if not filtered:
        return []
    
    filtered.sort(key=lambda p: origin.distance(p)) #sort by distance from origin, depths that you can see
    results = []
    
    step_back_distance = 1e-6
    hit_wall = False #track if you have hit an edge, because you will not see past it
         
    for p in filtered:
        if hit_wall:
            break

        direction = (p.x - origin.x, p.y - origin.y)
        length = math.hypot(*direction)
        if length == 0:
            continue  #skip degenerate

        unit = (direction[0] / length, direction[1] / length)
        test_point = Point(p.x - step_back_distance * unit[0], p.y - step_back_distance * unit[1])

        if not polygon.covers(test_point): #extra test, making sure we did not travel outside of polygon to reach this intersection
            hit_wall = True
            break  #skip if not a valid visible point

        
        is_vertex = False
        for index, v in enumerate(polygon.exterior.coords): #check if the point is a vertex of the polygon by comparing with all vertices
            if (round(p.x, 8), round(p.y, 8)) == (round(v[0], 8), round(v[1], 8)):
                is_vertex = True
                if is_reflex_vertex(polygon.exterior.coords, index): #if it is a reflex, there will be another intersection further away
                    results.append(p) 
                else:
                    results.append(p)
                    hit_wall = True  #regular vertex, stop looking for further intersections
                break

        if not is_vertex:
            results.append(p)
            hit_wall = True  #hit an edge, stop looking for further intersections

    if not results:
        return []
    return results



def visible_verts(polygon, i, epsilon=1e-8): #reutrns list of verts defining visibility polygon
    origin = Point(polygon[i])
    n = len(polygon)

    unique_angles = set()
    for j in range(n):
        if j == i:
            continue
        angle = angle_between(polygon[i], polygon[j])  #create list of angles toward each vertex
        unique_angles.add(angle)
        
 
    visible_points = [Point(polygon[(i-1)%n]), Point(polygon[(i+1)%n])] #initialize list with neighbors of origin
    poly = Polygon(polygon)

    for angle in unique_angles: #create ray in each direction, find valid intersections from it
        dx = math.cos(angle)
        dy = math.sin(angle)
        far_point = (origin.x + dx * 1000, origin.y + dy * 1000)
        ray = LineString([origin.coords[0], far_point])
        
        new_points = valid_intersections(ray, poly, origin)

        for pt in new_points:
            visible_points.append(pt) #add to visible verts
    
    sorted_vis_points = sort_ccw(visible_points, origin, poly.exterior.coords) #sort verts counter clockwise
    for index in range(len(sorted_vis_points)):
        pt = sorted_vis_points[index]
        next_pt = sorted_vis_points[(index+1)%len(sorted_vis_points)] #find where to insert origin, between its neighbors in original polygon
        if  pt.equals(Point(polygon[(i-1)%n])) and next_pt.equals(Point(polygon[(i+1)%n])):
            sorted_vis_points.insert(index+1, Point(polygon[i]))
            break
    
    if len(sorted_vis_points) < 4: #shape might be triangle, which would cause issues later
        midpt = midpoint(sorted_vis_points[len(sorted_vis_points)-1], sorted_vis_points[0])
        sorted_vis_points.append(Point(midpt)) #add midpoint of an edge without changing overall shape
        
    return sorted_vis_points


def print_visibility_comparison(vis_polygons): #debugging tool for duplicate polygons
    for i in range(len(vis_polygons)):
        for j in range(i+1, len(vis_polygons)):
            poly1 = vis_polygons[i][1]
            poly2 = vis_polygons[j][1]
            if poly1.equals(poly2):
                print(f"Polygon {vis_polygons[i][0]} and {vis_polygons[j][0]} are exactly equal")
            elif poly1.almost_equals(poly2, decimal=6):
                print(f"Polygon {vis_polygons[i][0]} and {vis_polygons[j][0]} are almost equal")
            else:
                print(f"Polygon {vis_polygons[i][0]} and {vis_polygons[j][0]} differ")


def decompose_visibility(polygon_coords): #take in polygon, return visibility regions
    vis_polygons = []
    
    for i in range(len(polygon_coords)): #get visibility polygon for each vertex
        vis_poly = visible_polygon(polygon_coords, i)
        vis_polygons.append((i, vis_poly))
    
    #normalize and compare polygons
    unique_polygons = {}
    for idx, poly in vis_polygons:
        norm = normalize_polygon(poly)
        
        found_duplicate = False
        for existing_idx, existing_poly in unique_polygons.values():
            if compare_polygons(norm, existing_poly):
                found_duplicate = True
                break
        
        if not found_duplicate:
            unique_polygons[norm] = (idx, poly)
    
    vis_polygons = unique_polygons.values()

    #attempt to build regions
    try:
        polygons = [vis_poly for _, vis_poly in vis_polygons]  #list of polygons from vis_polygons
        edge_set = set()
        for poly in polygons:
            coords = list(poly.exterior.coords)   
            for i in range(len(coords)-1):
                p1 = tuple(round(coord, 8) for coord in coords[i])
                p2 = tuple(round(coord, 8) for coord in coords[i+1])
                edge = tuple(sorted((p1, p2)))  #undirected edge
                edge_set.add(edge)

        edges = [LineString([e[0], e[1]]) for e in edge_set]


        merged_edges = unary_union(edges)
        regions = list(polygonize(merged_edges))  #make sure this is a list


    except Exception as e:
        print("Polygonization failed with error:", e)
        return []

    #if polygonization succeeded but returned nothing, try fallback
    if not regions:
        print("No regions created from polygonize: trying fallback with unary_union...")
        try:
            merged_edges = unary_union(edges)
            regions = list(polygonize(merged_edges))
            print(f"Fallback succeeded. Created {len(regions)} regions.")
        except Exception as e:
            print("Fallback polygonization failed with error:", e)
            return []

    labeled_regions = [] #go through regions, assign vibility sets
    for region in regions:
        visible_by = []
        for idx, vis_poly in vis_polygons:
            centroid = region.centroid
            rounded_centroid = Point(round(centroid.x, 6), round(centroid.y, 6))
            if vis_poly.distance(rounded_centroid) < 1e-8:
                visible_by.append('1')
            else:
                visible_by.append('0') 
        visible_by = ''.join(visible_by)
        visible_by = int(visible_by, 2)
        
        if visible_by ==0:
            continue
        
        labeled_regions.append((region, visible_by))
        
    merged_labeled_regions = merge_regions_if_same_visibility(labeled_regions) #simplify/remove rounding errors

    return merged_labeled_regions

def merge_regions_if_same_visibility(regions, visibility_threshold=0):
    merged_regions = []
    visited = set()  #track which regions have already been merged
    
    for i, (region_i, vis_i) in enumerate(regions):
        if i in visited:
            continue  #skip already merged regions
        
        #initialize list with the current region
        merged_region = region_i
        merged_visibility = vis_i
        visited.add(i)
        
        #check for neighboring regions with the same visibility
        for j, (region_j, vis_j) in enumerate(regions):
            if j != i and j not in visited and vis_i == vis_j:
                #if they have the same visibility, merge the regions
                merged_region = unary_union([merged_region, region_j])
                visited.add(j)
        
        merged_regions.append((merged_region, merged_visibility))
    
    return merged_regions




#Smaller Helper Functions:
def angle_between(p1, p2): #calculate angle between point coords
    return math.atan2(p2[1] - p1[1], p2[0] - p1[0]) % (2*math.pi)

def normalize_polygon(poly, tolerance=1e-6):#normalize/round the coordinates to reduce precision errors
    norm_coords = [tuple(round(coord, 6) for coord in point) for point in poly.exterior.coords]
    normalized_poly = Polygon(norm_coords)
    
    return normalized_poly

def compare_polygons(poly1, poly2, tolerance=1e-6):#compare exterior coordinates of both polygons with a tolerance
    coords1 = list(poly1.exterior.coords)
    coords2 = list(poly2.exterior.coords)
    
    if len(coords1) != len(coords2):
        return False  # different number of coordinates
    
    for c1, c2 in zip(coords1, coords2):
        if not all(abs(a - b) < tolerance for a, b in zip(c1, c2)):
            return False  #coordinates differ by more than tolerance
    
    return True

def is_edge_intersect(point, points): #check if an intersect is on edge, and not vert
    for i in range(len(points) - 1):  #len - 1 because polygon is closed (last == first)
        v1 = points[i]
        v2 = points[(i + 1)%len(points)]
        edge = LineString([v1, v2])

        if edge.distance(point) <= 1e-8: #check if point lies on the edge (with tolerance)
            if not (point.equals_exact(Point(v1), 1e-8) or point.equals_exact(Point(v2), 1e-8)): #make sure not one of the verts of the edge
                return edge #return the edge
        
    return None

def is_on_edge(point, points, edge): #given an edge, check if pt is on it
    if edge.distance(point) <= 1e-8:
        return True
    else:
        return False
 
def angle_and_distance_from(p, center): #get angle and dist between pts
    angle = math.atan2(p.y - center.y, p.x - center.x)
    distance = center.distance(p)
    angle = round(angle,8)
    return (angle, distance)

def sort_ccw(points, center, poly_points): #return verts in order ready to become a polygon
    preliminary_sorted = sorted(points, key=lambda p: angle_and_distance_from(p, center)) #sort by angle, use distance as secondary

    final_sorted = []
    for i, p in enumerate(preliminary_sorted): #now go through to make sure two intersects on same ray are handled well
        if p in final_sorted: #if it was added as next_p before, skip
            continue
        angle_i = math.atan2(p.y - center.y, p.x - center.x) 
        next_p = preliminary_sorted[(i+1)%len(preliminary_sorted)]
        angle_next = math.atan2(next_p.y - center.y, next_p.x - center.x)
        prev_p = preliminary_sorted[(i-1)%len(preliminary_sorted)]
        if abs(angle_i-angle_next) < 1e-8: #checking for cases with two intersections on same ray
            edge = is_edge_intersect(p, poly_points) #if the first intersect is on an edge (not a vertex), check its neighbor on the other side
            if edge:
                if not is_on_edge(prev_p, poly_points, edge): #if the neighbor on its other side is not on the same edge, they do not belong next to each other; swap intersect order
                    final_sorted.append(next_p)
                    final_sorted.append(p)
                else:
                    final_sorted.append(p) #if the other neighbor is on the same edge, they belong together
                    final_sorted.append(next_p)
            else:
                edge = is_edge_intersect(next_p, poly_points) #do the same if the next point is the edge point
                if edge:
                    if is_on_edge(prev_p, poly_points, edge):
                        final_sorted.append(next_p)
                        final_sorted.append(p)
                    else:
                        final_sorted.append(p)
                        final_sorted.append(next_p)            
        else:
            final_sorted.append(p) #if no special relationship, add and move on
    return final_sorted

def midpoint(p1, p2):#calculate the midpoint of the line segment joining p1 and p2
    mx = (p1.x + p2.x) / 2
    my = (p1.y + p2.y) / 2
    return (mx, my)

def visible_polygon(polygon_coords, i): #return visibility polygon of vert i
    return Polygon(visible_verts(polygon_coords, i))

def is_reflex_vertex(vertices, i):
    #remove duplicate last point
    if vertices[0] == vertices[-1]:
        vertices = vertices[:-1]

    a = Point(vertices[i - 1]) #check witth prev and next pts
    b = Point(vertices[i])
    c = Point(vertices[(i + 1) % (len(vertices))])

    ab = (a.x - b.x, a.y - b.y)
    cb = (c.x - b.x, c.y - b.y)

    cross = ab[0] * cb[1] - ab[1] * cb[0] #get angle between three pts
    return cross > 0  #reflex if angle turns right (assuming CCW polygon)
