from collections import defaultdict
from consts import CLUSTERED_SUBSTITUTIONS_MIN_SIZE
from consts import WTS_MIN_PERCENT
import math
from itertools import product
from scipy.special import comb
import numpy as np


def binom(n, k, p):
    return comb(n, k, exact = True)  * (p ** k) * (1 - p) ** (n - k)

def generate_freqs(counts):

    ranges = []

    for i, count in enumerate(counts):
        ranges.append(range(count + 1))

    return list(product(*ranges))

def binom_from(size, freq, p):
    return sum(binom(size, freq, p) for freq in range(freq, size + 1 ))

class SNDWindow:

    def __init__(self, index, snds, window_size, number_of_bins):
        self.window_size = window_size
        self.number_of_bins = number_of_bins
        self.bin_size = window_size // self.number_of_bins
        self.index = index
        self.window = defaultdict(list)
        self.window_counts = [0] * (2 * number_of_bins - 1)
        self.window_biased_counts = [0] * (2 * number_of_bins - 1)
        self.compressed_window_counts = None
    
        self.middle_coord = snds[index].target_coord() 
        
        self.middle_coord_start =  self.middle_coord // self.bin_size * self.bin_size
        self.middle_coord_bin = self.middle_coord // self.bin_size 
          
        start_window_coord =  max(0, self.middle_coord_start - self.window_size + self.bin_size)
        end_window_coord =  self.middle_coord_start + self.window_size - 1

        start = self.index
        end = self.index
        
        #print(start_window_coord, end_window_coord,  self.middle_coord )
        
        while start >=0 and snds[start].target_coord()  >= start_window_coord:
            start -= 1

        while end < len(snds) and snds[end].target_coord()  <= end_window_coord:
            end += 1
            
        self.number_of_snds = end - start - 1
        
        if self.number_of_snds < CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
            return

        for snd in snds[start + 1: end]:
            bin_coord = snd.target_coord()  // self.bin_size - self.middle_coord_bin + self.number_of_bins - 1
            self.window[bin_coord].append(snd)
            self.window_counts[bin_coord] += 1
            
            if snd.biased():
                self.window_biased_counts[bin_coord] += 1
                
    def get_number_of_snds_in_cluster(self, i):
        return sum([count for count in self.window_counts[i:i + self.number_of_bins]])
    
    def get_number_of_biased_snds_in_cluster(self, i):
        return sum([count for count in self.window_biased_counts[i:i + self.number_of_bins]])
    
    def is_clustered(self):
        
        if self.number_of_snds < CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
            return False

        for i in range(self.number_of_bins):
            if self.get_number_of_snds_in_cluster(i) >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
                return True

        return False

    def is_biased_clustered(self):
        
        if self.number_of_snds < CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
            return False

        for i in range(self.number_of_bins):
            number_of_snds = self.get_number_of_snds_in_cluster(i)
            number_of_biased_snds = self.get_number_of_biased_snds_in_cluster(i)
            
            if number_of_snds >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE and number_of_biased_snds / number_of_snds >= WTS_MIN_PERCENT:
                return True
        return False
    
    def get_middle_coord(self):
        return self.middle_coord

    def get_number_of_bins(self):
        return self.number_of_bins

    def compress_window_counts(self):

        if self.number_of_snds < CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
            self.compressed_window_counts = []
            self.compressed_number_of_bins = 0
            return 

        starts = []
        indexes = set()

        for i in range(self.number_of_bins):

            if sum(self.window_counts[i:i + self.number_of_bins]) >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE:

                ind = tuple( i + j for j, v in enumerate(self.window_counts[i:i + self.number_of_bins]) if v > 0)

                if ind not in indexes:
                    indexes.add(ind)
                    starts.append(i)

        starts = starts + [ start + self.number_of_bins for start in starts ]

        self.compressed_window_counts =  [ sum(self.window_counts[start: end]) for start, end in zip(starts, starts[1:])]
        self.compressed_number_of_bins = len(self.compressed_window_counts) // 2 + 1

    def get_compressed_window_counts(self):

        if not self.compressed_window_counts:
            self.compress_window_counts()
        
        return self.compressed_window_counts   
  
    def get_max_frequency_of_cluster(self):
        self.get_compressed_window_counts()
       
        if self.compressed_number_of_bins:
            return max([sum(self.compressed_window_counts[i:i + self.compressed_number_of_bins]) for i in range(self.compressed_number_of_bins)])
        else:
            return 0

    def get_prob_of_bcs(self, p):

        if self.bin_size == 1:
            bins = self.get_compressed_window_counts()
        else:
            bins = self.window_counts

        if not bins:
            return 0

        n = len(bins) // 2 + 1

        first_cluster_size = sum(bins[:n])
	
        prob = 0

        if p == 0:
            print(bins)

        if first_cluster_size >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
            prob = binom_from(first_cluster_size, math.ceil(WTS_MIN_PERCENT * first_cluster_size), p)

        mem = defaultdict(lambda: 1)

        for k in range(1, n):

            conditional_counts = bins[k : k + n - 1]
            cluster_size = sum(bins[k : k + n])
            prev_cluster_size = sum(bins[k - 1 : k + n - 1])

            prev_mem = mem
            mem = defaultdict(int)

            for conditional_freqs in generate_freqs(conditional_counts):

                conditional_freqs_size = sum(conditional_freqs)
                conditional_freqs_prob = np.prod([binom(count, freq, p) for count, freq in zip(conditional_counts, conditional_freqs)])

                n_a = 0
           
                a = binom_from(bins[k + n - 1], max(math.ceil(WTS_MIN_PERCENT * cluster_size) -  conditional_freqs_size, 0), p)

                upper_bound =  min( math.ceil(WTS_MIN_PERCENT * prev_cluster_size) - conditional_freqs_size, bins[k - 1] + 1)

                if prev_cluster_size < CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
                    upper_bound = bins[k - 1] + 1

                for freq in range(0 , upper_bound):
                    n_a += binom(bins[k - 1], freq, p)  *  prev_mem[ (freq, ) + conditional_freqs[:-1] ]

                if cluster_size >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
                    prob += a * n_a * conditional_freqs_prob

                mem[conditional_freqs] = n_a

        return prob
                                                                                                                  
