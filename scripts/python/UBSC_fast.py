import math
from itertools import combinations
from itertools import product
from scipy.special import comb
import numpy
from collections import defaultdict

def binom(n, k, p):
    return comb(n, k, exact = True)  * (p ** k) * (1 - p) ** (n - k)  

def generate_freqs_with_prob(counts, p):
    
    def generate(counts, p):
        if len(counts) == 0:
            return [[[], 1]]
        else:
            results = []

            for first_freq in range(counts[0] + 1):
                freqs_with_prob =  generate(counts[1:], p)

                for freqs, prob in freqs_with_prob:
                    results.append( ([first_freq] + freqs, prob *  binom(counts[0], first_freq ,p)))

            return results

    return [ (tuple(freqs), prob) for freqs, prob  in generate(counts, p) ] 


def binom_from(size, freq, p):
    return sum(binom(size, freq, p) for freq in range(freq, size + 1 ))


def UBCS_fast(shifts, p):

    n = len(shifts) // 2 + 1
    
    if shifts[n - 1] == 0:
        return 0
    
    prob = 0
    
    first_cluster_size = sum(shifts[:n])
    
    if first_cluster_size > 4:
        prob = binom_from(first_cluster_size, math.ceil(0.8 * first_cluster_size), p)
        
    mem = defaultdict(lambda: 1)
        
    for k in range(1, n):   
    
        conditional_counts = shifts[k : k + n - 1]
        cluster_size = sum(shifts[k : k + n])
        prev_cluster_size = sum(shifts[k - 1 : k + n - 1])
        
        prev_mem = mem
        mem = defaultdict(int)
        
        for conditional_freqs, conditional_freqs_prob in generate_freqs_with_prob(conditional_counts, p):
            conditional_freqs_size = sum(conditional_freqs)
            
            n_a = 0 
            a = 0
            
            a = binom_from(shifts[k + n - 1], math.ceil(0.8 * cluster_size) -  conditional_freqs_size, p)
            
            upper_bound =  min( math.ceil(0.8 * prev_cluster_size) - conditional_freqs_size , shifts[k - 1] + 1)
            
            if prev_cluster_size < 5:
                upper_bound = shifts[k - 1] + 1
            
            for freq in range(0 , upper_bound):
                n_a += binom(shifts[k - 1], freq, p)  *  prev_mem[ (freq, ) + conditional_freqs[:-1] ]

            if cluster_size > 4:
                prob += a * n_a * conditional_freqs_prob
                
            mem[conditional_freqs] = n_a
            
    return prob
