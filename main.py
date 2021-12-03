#!/usr/bin/env python
# coding: utf-8

import sys
import math
import random
import copy
import itertools
import time

import Levenshtein # pip install python-Levenshtein # Levenshtein.distance(s1,s2) - fastest
#import editdistance # pip install editdistance # editdistance.eval(s1,s2) - pretty fast

# levenshtein edit distance implementation
# Heikki H., "Explaining and extending the bit-parallel approximate string matching algorithm of Myers", (2001).
def edit_distance(s1, s2):
    m=len(s1)+1
    n=len(s2)+1

    tbl = {}
    for i in range(m): tbl[i,0]=i
    for j in range(n): tbl[0,j]=j
    for i in range(1, m):
        for j in range(1, n):
            cost = 0 if s1[i-1] == s2[j-1] else 1
            tbl[i,j] = min(tbl[i, j-1]+1, tbl[i-1, j]+1, tbl[i-1, j-1]+cost)

    return tbl[i,j]

# Modified version of itertools.product (with C G consideration)
def xproduct(*args, repeat=1, CG_min=0.45, CG_max=0.55):
    at_max = (1-CG_min) * repeat
    cg_max = CG_max * repeat
    pools = [tuple(pool) for pool in args] * repeat
    result = [[]]
    for i,pool in enumerate(pools):
        #result = [x+[y] for x in result for y in pool]
        tmp_result = []
        for x in result:
            for y in pool:
                check_cg = True
                # CG_min check
                if i >= at_max:
                    a_count = ''.join(x+[y]).count('A')
                    t_count = ''.join(x+[y]).count('T')
                    if a_count+t_count > math.floor(at_max):
                        check_cg = False
                # CG_max check
                if i >= cg_max:
                    c_count = ''.join(x+[y]).count('C')
                    g_count = ''.join(x+[y]).count('G')
                    if c_count+g_count > math.floor(cg_max):
                        check_cg = False
                if check_cg:
                    tmp_result += [x+[y]]
        result = tmp_result
    for prod in result:
        yield tuple(prod)


# Sorting strands based on sum of edit distance of the rest
def strands_sorting_ed_sum(strands_arr):
    ed_arr = []
    for i in range(len(strands_arr)):
        total = 0
        for j in range(len(strands_arr)):
            if i != j:
                total += Levenshtein.distance(strands_arr[i], strands_arr[j])
        ed_arr.append(total)

    # Sort the array based on edit distance - order by least similar (to the rest) first
    ed_arr, strands_arr = zip(*sorted(zip(ed_arr, strands_arr),reverse=True))
    strands_arr = list(strands_arr)
    # Return sorted strands
    return strands_arr

# Filter out those strands that is within ED threshold limit - evaluating from the first sorted item
def final_strands_ed_filter(strands_arr, D, L):
    final_strands_arr = []
    while (len(strands_arr) > 0):
        #n = random.randint(0, len(strands_arr)-1)
        #n=math.floor(len(strands_arr)/2.0)
        n=0
        strand = strands_arr[n]
        strands_arr.pop(n)
        #print(strand)
        # 6. Check if this strand has lower Edit Distance (similar) than the Threshold - otherwise remove
        ed_limit = True
        for s in final_strands_arr:
            ed_score = Levenshtein.distance(strand, s)
            if ed_score < D*L:
                ed_limit = False
                break
        if ed_limit:
            final_strands_arr.append(strand) 
    return final_strands_arr

# check if strands is within CG limit threshold
def check_strand_CG(s, CG_min, CG_max, L):
    if float(s.count('C')+s.count('G'))/L >= CG_min and float(s.count('C')+s.count('G'))/L <= CG_max:
        return True
    return False
        
    
# Check if strand is within threshold limit with the stored items
def check_strand_ed(strands_arr, s, D, L):
    for strand in strands_arr:
        ed_score = Levenshtein.distance(strand, s)
        if ed_score < D*L:
            return False
    return True


# # Method 1: Brute Force
# Method 1: Brute Force
# Method 1: Brute Force
def method_1(char_arr, L, CG_min, CG_max, D, sort=False):
    # 1. Permutate all possible number of samples of length 
    # O(4^L)
    strands_arr = [''.join(x) for x in list(itertools.product(char_arr, repeat=L))]
    print('Number of possible strands:', len(strands_arr))
    # 2. Remove those that does not fall under the CG limits
    # O(4^L)
    strands_arr = [x for x in strands_arr if float(x.count('C')+x.count('G'))/L >= CG_min and float(x.count('C')+x.count('G'))/L <= CG_max]
    print('Number of strands after (removing CG threshold):', len(strands_arr))

    # 3. Need to Pick intelligently the sequence of strands to put in the safe_lists (most disimilar)
    # Measure the sum of edit distance for 1 strand to the rest of the sample in population
    # Reason being: I want to evaluate those strands from the most similar
    # O(S^2)
    # 4. Sort the array based on edit distance - order by least similar (to the rest) first
    if sort:
        strands_arr = strands_sorting_ed_sum(strands_arr)
    #print(strands_arr)
    
    # 5. Iterate whole item in the list and put in safe_lists 
    # O(S * Log S)
    final_strands_arr = final_strands_ed_filter(strands_arr, D,L)

    # 6. Return the final list
    print('Number of strands after (removing ed limit):', len(final_strands_arr))
    return final_strands_arr


# # Method 2: Brute Force*
# Method 2: Brute Force*
def method_2(char_arr, L, CG_min, CG_max, D, sort=False):
    # 1. Permutate selectively based on CG count
    # O(L log(L))
    strands_arr = [''.join(x) for x in list(xproduct(char_arr, repeat=L, CG_min=CG_min, CG_max=CG_max))] 
    print('Number of possible strands (with CG threshold):', len(strands_arr))
    
    # 2. Need to Pick intelligently the sequence of strands to put in the safe_lists (most disimilar)
    # Measure the sum of edit distance for 1 strand to the rest of the sample in population
    # Reason being: I want to evaluate those strands from the most similar
    # O(S^2)
    # 3. Sort the array based on edit distance - order by least similar (to the rest) first
    if sort:
        strands_arr = strands_sorting_ed_sum(strands_arr)
        
    # 4. Iterate whole item in the list and put in safe_lists 
    # O(S * Log S)
    final_strands_arr = final_strands_ed_filter(strands_arr, D,L)
    

    # 3. Return the final list
    print('Number of strands after (removing ed limit):', len(final_strands_arr))
    return final_strands_arr


# # Method 3: Population Generation
# Method 3: population Generation
def method_3(char_arr, L, CG_min, CG_max, D):
    final_strands_arr = []
    # 1. Pick first set of strands: that preserve the CG_min, CG_max
    subL = int(((CG_max + CG_min) / 2.0) * L)
    # 2. First part is only for C and G, Second part is only for A and T
    char_div2 = math.floor(len(char_arr)*1.0/2)
    sub_char1_arr = char_arr[0:char_div2]
    sub_char2_arr = char_arr[char_div2:]
    rchar = {sub_char1_arr[0]:sub_char1_arr[1], sub_char1_arr[1]:sub_char1_arr[0], sub_char2_arr[0]:sub_char2_arr[1], sub_char2_arr[1]:sub_char2_arr[0] }
    # 3. Generate those strands that has lowest sum of ED score
    # Initial Population
    for s1 in sub_char1_arr:
        for s2 in sub_char2_arr:
            strand1 = ''.join([s1 for x in range(subL)]+[s2 for x in range(subL)])
            strand2 = ''.join([s2 for x in range(subL)]+[s1 for x in range(subL)])
            if (check_strand_ed(final_strands_arr, strand1, D, L)) and check_strand_CG(strand1, CG_min, CG_max, L):
                final_strands_arr.append(strand1)
            if (check_strand_ed(final_strands_arr, strand2, D, L)) and check_strand_CG(strand2, CG_min, CG_max, L):
                final_strands_arr.append(strand2)
            
    if L == 8:
        sub1_strands_arr = method_1(char_arr, 4, CG_min, CG_max, D)
        for s1 in sub1_strands_arr:
            for s2 in sub1_strands_arr:
                tmp = s1+s2
                if (check_strand_ed(final_strands_arr, tmp, D, L)) and check_strand_CG(tmp, CG_min, CG_max, L):
                    final_strands_arr.append(tmp)
    
    if L == 10:
        sub1_strands_arr = method_1(char_arr, 6, CG_min, CG_max, D)
        sub2_strands_arr = method_1(char_arr, 4, CG_min, CG_max, D)
        for s1 in sub1_strands_arr:
            for s2 in sub2_strands_arr:
                tmp = s1+s2
                if (check_strand_ed(final_strands_arr, tmp, D, L)) and check_strand_CG(tmp, CG_min, CG_max, L):
                    final_strands_arr.append(tmp)
                tmp = s2+s1
                if (check_strand_ed(final_strands_arr, tmp, D, L)) and check_strand_CG(tmp, CG_min, CG_max, L):
                    final_strands_arr.append(tmp)
                    
    if L == 20:
        sub1_strands_arr = method_1(char_arr, 10, CG_min, CG_max, D)
        for s1 in sub1_strands_arr:
            for s2 in sub1_strands_arr:
                    tmp = s1+s2
                    if (check_strand_ed(final_strands_arr, tmp, D, L)) and check_strand_CG(tmp, CG_min, CG_max, L):
                        final_strands_arr.append(tmp)
    
    # Generation
    n_gen = 10
    n_shuffle = L*10 #1000
    n_mutate = L*10 #1000
    n_remove = 0.20 #0.25
    n_crossover = L*10 #1000
    
    n_pop = []
    best_generation = []
    print('Start Generation')
    for i in range(n_gen):
        strands = final_strands_arr
        print('-----------------')
        print('Generation ',i+1, ' Population: ', len(final_strands_arr))
        # I. Cross-Over
        for j in range(n_crossover):
            for s1 in strands:
                n = len(strands)
                s2 =  strands[random.randint(0, n-1)]
                # Randomly pick the cross-over point and cross-over length
                co_start = random.randint(0, L-1)
                co_end = random.randint(co_start, L-1)
                s1_char = s1[co_start:co_end]
                s2_char = s2[co_start:co_end]
                s1 = s1[:co_start] + s2_char + s1[co_end:]
                s2 = s2[:co_start] + s1_char + s2[co_end:]
                if (check_strand_ed(final_strands_arr, s1, D, L)) and check_strand_CG(s1, CG_min, CG_max, L):
                    final_strands_arr.append(s1)
                if (check_strand_ed(final_strands_arr, s2, D, L)) and check_strand_CG(s2, CG_min, CG_max, L):
                    final_strands_arr.append(s2)
        print('Population After Cross-over: ', len(final_strands_arr))
        
        # II. Shifting
        for j in range(1, L):
            for strand in strands:
                tmp = strand[j:L]+strand[0:j]
                if (check_strand_ed(final_strands_arr, tmp, D, L)) and check_strand_CG(tmp, CG_min, CG_max, L):
                    final_strands_arr.append(tmp)
        print('Population After Shift: ', len(final_strands_arr))
        
        
        # III. Random Mutation
        for j in range(n_mutate):
            for strand in strands:
                # Randomly pick the mutation severity
                tmp = list(strand)
                rand = random.randint(int(L*D), L-1)
                for j in range(rand):
                    # Randomly pick location of mutation
                    idx = random.randint(0, L-1)
                    tmp[idx] = random.choice(char_arr)
                tmp = ''.join(tmp)
                if (check_strand_ed(final_strands_arr, tmp, D, L)) and check_strand_CG(tmp, CG_min, CG_max, L):
                    final_strands_arr.append(tmp)
        print('Population After Mutation: ', len(final_strands_arr))
        
        # IV. Shuffling
        for j in range(n_shuffle):
            for strand in strands:
                tmp = list(strand)
                random.shuffle(tmp)
                tmp = ''.join(tmp)
                if (check_strand_ed(final_strands_arr, tmp, D, L)) and check_strand_CG(tmp, CG_min, CG_max, L):
                    final_strands_arr.append(tmp)
        print('Population After Shuffling: ', len(final_strands_arr))


        # X. Save best number of generation
        n_pop.append(len(final_strands_arr))
        if len(best_generation) < len(final_strands_arr):
            best_generation = copy.deepcopy(final_strands_arr)
            
        # L. Remove strand that has high ED
        final_strands_arr = strands_sorting_ed_sum(final_strands_arr)
        #final_strands_arr = final_strands_arr[math.floor(n_remove*len(final_strands_arr)):] # - WORST
        final_strands_arr = final_strands_arr[0:math.floor((1-n_remove)*len(final_strands_arr))] # - BETTER
        # Random removal of items in the list
        #for j in range(math.floor((n_remove)*len(final_strands_arr))):
        #    final_strands_arr.pop(random.randint(0, len(final_strands_arr)-1))
        print('Population After Removal: ', len(final_strands_arr))
        

    print('-------------')
    print(n_pop)
    print('-------------')
    if len(best_generation) > len(final_strands_arr):
        return best_generation
    return final_strands_arr



char_arr = ['A','T','C','G']
L = 8
CG_min = 0.45
CG_max = 0.55
D = 0.4
M = 1


if __name__ == "__main__":
    for i in range(1,len(sys.argv),2):
        if sys.argv[i] == '-l':
            L = int(sys.argv[i+1])
        elif sys.argv[i] == '-d':
            D = float(sys.argv[i+1])
        elif sys.argv[i] == '-min':
            CG_min = float(sys.argv[i+1])
        elif sys.argv[i] == '-max':
            CG_max = float(sys.argv[i+1])
        elif sys.argv[i] == '-m':
            M = int(sys.argv[i+1])
            
print('Length: ', L)
print('CG_min: ', CG_min)
print('CG_max: ', CG_max)
print('D: ', D)

start = time.time()
if M == 1:
    print('Run Method 1: Brute Force')
    final_strands_arr = method_1(char_arr, L, CG_min, CG_max, D)
if M == 2:
    print('Run Method 2: Brute Force*')
    final_strands_arr = method_2(char_arr, L, CG_min, CG_max, D, sort=True)
if M == 3:
    print('Run Method 3: Population Generation')
    final_strands_arr = method_3(char_arr, L, CG_min, CG_max, D)
    
end = time.time()
print('time: ', end - start)
print('Final # Strands:', len(final_strands_arr)) 

with open('dna_M'+str(M)+'_L'+str(L), 'w') as f:
    for s in final_strands_arr:
        f.write("%s\n" % s)
