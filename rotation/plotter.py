import matplotlib.pyplot as plt

def plot_polygon(ax, polygon, facecolor='lightgrey', edgecolor='black', alpha=0.4, label=None):
    if polygon.is_empty:
        return
    if polygon.geom_type == 'Polygon':
        x, y = polygon.exterior.xy
        ax.fill(x, y, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha, label=label)
        for interior in polygon.interiors:
            x, y = interior.xy
            ax.fill(x, y, facecolor='white', edgecolor=edgecolor, alpha=1.0)
    elif polygon.geom_type == 'MultiPolygon':
        for p in polygon.geoms:
            plot_polygon(ax, p, facecolor, edgecolor, alpha)

def plot_visibility_debug(outer_polygon, guards, visibility_polygons, regions, region_visibility=None):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')

    # Plot base polygon
    plot_polygon(ax, outer_polygon, facecolor='lightgrey', edgecolor='black', alpha=0.4, label='Main Polygon')

    # Plot guard positions
    for i, guard in enumerate(guards):
        ax.plot(guard.x, guard.y, 'ro')
        ax.text(guard.x, guard.y, f"G{i}", fontsize=9, ha='center', va='bottom', color='red')

    # Plot visibility polygons
    for i, vis_poly in enumerate(visibility_polygons):
        plot_polygon(ax, vis_poly, facecolor='lightblue', edgecolor='blue', alpha=0.3)
        centroid = vis_poly.centroid
        ax.text(centroid.x, centroid.y, f"V{i}", fontsize=9, color='blue')
    
    # Plot decomposed visibility regions
    for i, region in enumerate(regions):
        plot_polygon(ax, region, facecolor='green', edgecolor='darkgreen', alpha=0.3)
        centroid = region.centroid
        ax.plot(centroid.x, centroid.y, 'kx')
        label = f"R{i}"
        if region_visibility is not None:
            label += f"\n{region_visibility[i]}"
        ax.text(centroid.x, centroid.y, label, fontsize=8, ha='left', va='center', color='black')
    
    plt.title("Visibility Regions and Guards")
    plt.legend()
    plt.grid(True)
    plt.show()


def extract(method, results, sizes):
    guards = []
    times = []
    for n in sizes:
        entries = results[n][method]
        # Take average over all trials
        if entries and all(e[0] is not None for e in entries):
            g = sum(e[0] for e in entries) / len(entries)
            t = sum(e[1] for e in entries) / len(entries)
        else:
            g, t = None, None
        guards.append(g)
        times.append(t)
    return guards, times


def plot_results(results, sizes):
    
    guards_region, times_region = extract("greedy_region", results, sizes)
    guards_guard, times_guard = extract("greedy_guard", results, sizes)
    guards_brute, times_brute = extract("brute_force", results, sizes)
    guards_hybrid, times_hybrid = extract("hybrid", results, sizes)
    
    plt.figure(figsize=(10, 5)) #plotting solution size
    plt.plot(sizes, guards_region, 'o-', label='Greedy Region Covering')
    plt.plot(sizes, guards_guard, 's-', label='Greedy Guard Removal')
    plt.plot(sizes, guards_brute, 'x-', label='Brute Force (Optimal)')
    plt.plot(sizes, guards_hybrid, '*-', label='Hybrid Greedy')
    plt.xlabel("Polygon Size (n)")
    plt.ylabel("Number of Guards")
    plt.title("Number of Guards vs Polygon Size")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("solution_size")
    plt.show()

    plt.figure(figsize=(10, 5)) #plotting execution time
    plt.plot(sizes, times_region, 'o-', label='Greedy Region Addition')
    plt.plot(sizes, times_guard, 's-', label='Greedy Guard Removal')
    plt.plot(sizes, times_brute, 'x-', label='Brute Force (Optimal)')
    plt.plot(sizes, times_hybrid, '*-', label='Greedy Hybrid')
    plt.xlabel("Polygon Size (n)")
    plt.ylabel("Execution Time (s)")
    plt.title("Execution Time vs Polygon Size")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("computation_times")
    plt.show()

def plot_regions(region_info, sizes):
    ave_num_regions = [] #plot individual values as pts, and averages on line
    for s in sizes:
        total = 0
        num = 0
        for i, r in enumerate(region_info[0]):
            if region_info[1][i] == s:
                total += r
                num += 1
        ave_num_regions.append(total/num)

    plt.plot(sizes, ave_num_regions, '-', label='Average Number of Regions')           
    plt.scatter(region_info[1], region_info[0], marker='x')
    plt.xlabel("Polygon Size (n)")
    plt.ylabel("Number of Regions")
    plt.title("Number of Visibility Regions vs Polygon Size")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("num_regions")
    plt.show()