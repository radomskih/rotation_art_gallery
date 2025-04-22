
from polygenerator import random_polygon
from itertools import combinations
import time
from plotter import *
from visibility import *


def is_complete(regions, configuration): #checks if visibility over time is complete for all regions
    config_mask = 0
    for i, b in enumerate(configuration):
        if b == 1:
            config_mask |= (1 << i)

    for (_, VT) in regions:
        if (VT & config_mask) == 0:
            return False
    return True


def greedy_region(regions, n):   
    configuration = [1]+[0]*(n-1) #start with one guard
    VTs = []
    for (region, vis_set) in regions: #get visibility over timewith one guard
        VTs.append(vis_set)
    return place_guards_greedy(regions, VTs, configuration) #get the final guard configuration
        

def place_guards_greedy(regions, VTs, configuration):  
    complete = is_complete(regions, configuration)
    while not complete:
        n = len(regions)   
        incomplete_regions = []
        incomplete_VTs = []
        VT_index = []
        for i, VTr in enumerate(VTs): #first check if a new guard is needed
            if VTr != 2**n - 1: 
                incomplete_regions.append(regions[i]) #make a list of regions that need more guards
                incomplete_VTs.append(VTr)
                VT_index.append(i)
   
        max_diff = 0
        best_guard = 0
        for s, b in enumerate(configuration): #then try the guard positions that are still empty
            if b == 0:
                total_diff = 0    
                for i, (region, vis_set) in enumerate(incomplete_regions):
                    VT_extended = (incomplete_VTs[i] << n) |  incomplete_VTs[i] #extend to 2 cycles for overlap
                    VT_proposed = (vis_set << s) | VT_extended 
                    VT_diff = VT_extended ^ VT_proposed
                    total_diff += bin(VT_diff).count("1") #sum up the changed bits
                    if total_diff > max_diff: #find the guard placement with the most improvement
                        max_diff = total_diff
                        best_guard = s
    
        configuration[best_guard] = 1 #add the best guard to the configuration
        for i,(region, vis_set) in enumerate(incomplete_regions): #update all the VTs that werent already complete
            VT_extended = (incomplete_VTs[i] << n) | incomplete_VTs[i]
            new_VT = (vis_set << best_guard | VT_extended)
            index = VT_index[i]
            VTs[index] = new_VT
        complete = is_complete(regions, configuration)
    return configuration
  

def build_visibility_matrix(regions, n): #set up matrix based on region visibility
    matrix = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(len(regions))]

    for r, (_, VT) in enumerate(regions):
        for g in range(n):
            rotated = ((VT << g) | (VT >> (n - g))) & ((1 << n) - 1)
            for t in range(n):
                if (rotated >> t) & 1:
                    matrix[r][g][t] = 1  # Mark visibility for guard g at time t in region r
    return matrix


def guard_redundancy_score(matrix, configuration, guard_index):
    n = len(configuration)
    score = 0
    is_covered = True

    for r in range(len(matrix)):
        for t in range(n):
            if matrix[r][guard_index][t] == 1:
                temp = score #note the redundancy score before checking other guards
                # Check if other guards are covering this (r, t)
                for g in range(n):
                    if g != guard_index and configuration[g] == 1 and matrix[r][g][t] == 1:
                        score+=1
                if temp == score: #if the score is unchanged, no other guards see this region at this time, so no redundancy
                    is_covered = False
                    score = -1
                    return score
    return score


def most_redundant_guard(matrix, configuration): #find the guard with the most overlapping visibility
    max_score = -1
    most_redundant = -1
    for g in range(len(configuration)):
        if configuration[g] == 1:
            score = guard_redundancy_score(matrix, configuration, g)
            if score > max_score:
                max_score = score
                most_redundant = g
    return most_redundant, max_score


def remove_guards_greedy(matrix, regions, configuration):
    while True:
        g, score = most_redundant_guard(matrix, configuration) #find most redundant guard
        if g == -1:
            break
        #temporarily remove it
        configuration[g] = 0
        if not is_complete(regions, configuration):
            configuration[g] = 1  #restore if coverage is broken
            break
    return configuration


def greedy_guard(regions, n):    
    configuration = [1]*(n) #start with guards at every vertex
    VTs = []
    matrix = build_visibility_matrix(regions, n)
    return remove_guards_greedy(matrix, regions, configuration) #get the final guard configuration


def hybrid_greedy(regions, n): #use both approaches
    region_results = greedy_region(regions, n)    
    matrix = build_visibility_matrix(regions, n)
    return remove_guards_greedy(matrix, regions, region_results)

def brute_force(regions, n): #check every possible configuration until one works
    for r in range(1, n+1):
        for indices in combinations(range(n), r):
            configuration = [0]*n
            for i in indices:
                configuration[i] = 1   
            if is_complete(regions, configuration):
                return configuration
    return None


def run_tests(sizes, trials):
    results = {n: {"greedy_region": [], "greedy_guard": [], "brute_force": [], "hybrid": []} for n in sizes}
    num_regions = []
    indiv_sizes = []
    for n in sizes: #run trials for different sized polygons
        for trial in range(trials):    
            polygon = random_polygon(num_points=n)
            #debugging:
            #polygon = [(0,0), (4,0), (4,4), (0,4)]
            #polygon = [(0, 0.8615804574574811), (0.674019340409448, 0.7927078730944127), (0.3881692020560385, 0.3775882112632486), (0.8600753885663877, 0), (1, 0.8950553778651332), (0.3370343302471497, 1)]

            regions = decompose_visibility(polygon) #pre-process visibility regions
            num_regions.append(len(regions))
            indiv_sizes.append(n)
            
            '''
            debugging purposes:
            outer_polygon = Polygon(polygon)
            guards = [Point(p) for p in polygon]
            visibility_polygons = [visible_polygon(polygon, i) for i in range(len(polygon))]
            region_geometries = [region for region, _ in regions]
            region_visibility = [bin(vis)[2:].zfill(len(polygon)) for _, vis in regions]
            plot_visibility_debug(outer_polygon, guards, visibility_polygons, region_geometries, region_visibility)
            '''

            start = time.perf_counter() #run algorithms and store results and time spent
            greedy_region_soln = greedy_region(regions, n)
            results[n]["greedy_region"].append((sum(greedy_region_soln), time.perf_counter()-start))
            #print("found greedy region soln, size " , n, " trial ", trial)

            start = time.perf_counter()
            greedy_guard_soln = greedy_guard(regions, n)
            results[n]["greedy_guard"].append((sum(greedy_guard_soln), time.perf_counter()-start))
            #print("found greedy guard soln, size " , n, " trial ", trial)
            

            start = time.perf_counter()
            brute_force_soln = brute_force(regions, n)
            if brute_force_soln:
                results[n]["brute_force"].append((sum(brute_force_soln), time.perf_counter()-start))
            else:
                results[n]["brute_force"].append((None, time.perf_counter()-start))   
            #print("found brute force soln, size " , n, " trial ", trial)
            #print("number guards: ", sum(brute_force_soln))    
            

            start = time.perf_counter()
            hybrid_soln = hybrid_greedy(regions, n)
            if hybrid_soln:
                results[n]["hybrid"].append((sum(hybrid_soln), time.perf_counter()-start))
            else:
                results[n]["hybrid"].append((None, time.perf_counter()-start))
            #print("found hybrid soln, size ", n, " trial ", trial)
    region_info = [num_regions, indiv_sizes]       
    return results, region_info



if __name__ == "__main__":
    sizes = [4, 6, 8, 10, 12] #set test sizes
    results, region_info = run_tests(sizes, 5) #run tests and collect results
    plot_results(results, sizes) #plot
    plot_regions(region_info, sizes)
    
        
            