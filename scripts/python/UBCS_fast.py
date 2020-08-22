import math
from itertools import combinations
from itertools import product
from scipy.special import comb
import numpy
from collections import defaultdict
from consts import CLUSTERED_SUBSTITUTIONS_MIN_SIZE 
from consts import WTS_MIN_PERCENT

def compress_snd_window(window):
    n = len(window) // 2 + 1
    
    window = [ int(snd is not None) for snd in window]

    starts = []
    indexes = set()
    
    for i in range(n):
        
        if sum(window[i:i + n]) >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
            
            ind = tuple( i + j for j, v in enumerate(window[i:i + n]) if v == 1)
            
            if ind not in indexes:
                
                indexes.add(ind)
                starts.append(i)
                
    starts = starts + [ start + n for start in starts ]

    return  [ sum(window[start: end]) for start, end in zip(starts, starts[1:])]

def get_snds_in_window(index, snds, WINDOW_SIZE):
    start =  end = index
    middle_coord = snds[index].coord
    
    while start >= 0 and middle_coord - snds[start].coord < WINDOW_SIZE:
        start -= 1
   
    while end < len(snds) and snds[end].coord - middle_coord < WINDOW_SIZE:
        end += 1
    
    return snds[start + 1: end]
    
    
def get_snd_window(index, snds, WINDOW_SIZE):
 
    window = [None] * (2 * WINDOW_SIZE - 1)
    middle_coord = snds[index].coord
    
    for snd in  get_snds_in_window(index, snds, WINDOW_SIZE):
        window[snd.coord - middle_coord + WINDOW_SIZE - 1] = snd
    
    return window

def contains_cluster(window):

    WINDOW_SIZE = len(window) // 2 + 1

    for i in range(WINDOW_SIZE):
        if len([snd for  snd in window[i:i + WINDOW_SIZE] if snd]) >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
            return True
        
    return False

def is_biased_clustered(window):
    WINDOW_SIZE = len(window) // 2 + 1

    for i in range(WINDOW_SIZE):
        number_of_snds = len([snd for  snd in window[i:i + WINDOW_SIZE] if snd])
        number_biased = len([snd for  snd in window[i:i + WINDOW_SIZE] if snd and snd.biased])
        if  number_of_snds >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE and number_biased / number_of_snds >= WTS_MIN_PERCENT:
            return True
        
    return False


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
    
    if first_cluster_size >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
        prob = binom_from(first_cluster_size, math.ceil(WTS_MIN_PERCENT * first_cluster_size), p)
        
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
            
            a = binom_from(shifts[k + n - 1], math.ceil(WTS_MIN_PERCENT * cluster_size) -  conditional_freqs_size, p)
            
            upper_bound =  min( math.ceil(WTS_MIN_PERCENT * prev_cluster_size) - conditional_freqs_size , shifts[k - 1] + 1)
            
            if prev_cluster_size < CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
                upper_bound = shifts[k - 1] + 1
            
            for freq in range(0 , upper_bound):
                n_a += binom(shifts[k - 1], freq, p)  *  prev_mem[ (freq, ) + conditional_freqs[:-1] ]

            if cluster_size >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
                prob += a * n_a * conditional_freqs_prob
                
            mem[conditional_freqs] = n_a
            
    return prob
