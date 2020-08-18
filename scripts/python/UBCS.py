import math
from itertools import combinations
from itertools import product
from scipy.special import comb
import numpy

from consts import STW_MAX_PERCENT

def generate_unbiased_freqs_list(bins_lengths , threshold):
    if len(bins_lengths) == 0:
        return [[]]
    else:
        results = []

        for first_bin_content in range(bins_lengths[0] + 1):

            if first_bin_content <= threshold:
                remining_intervals = generate_unbiased_freqs_list(bins_lengths[1:], threshold - first_bin_content)

                for interval in remining_intervals:
                    results.append([first_bin_content] + interval)
        
        return results


def get_expected_BCS(intervals, p):
    intervals.sort()
    n = len(intervals)
    exp_value = 0

    for i in range(len(intervals)):
        j = i

        while j < n and intervals[j][0] <= intervals[i][1]:
            j += 1

        for k in range(j - i + 1):
            int_combinations = combinations(intervals[i + 1:j], k)

            a = 0

            for int_combination in int_combinations:
                int_exp_value = get_interception_expected_BCS([intervals[i]] + list(int_combination), p)
                a += 1

                if k % 2 == 0:
                    exp_value += int_exp_value
                else:
                    exp_value -= int_exp_value

    return exp_value


def get_interception_expected_BCS(intervals, p):
    intervals.sort()

    int_thresholds = [ math.floor(STW_MAX_PERCENT * (interval[1] - interval[0] + 1) ) for interval in intervals ]
    int_num = len(intervals)
 
    max_left = intervals[-1][0]
    min_right = intervals[0][1]

    if min_right < max_left: return 0
 
    bins_starts = [ interval[0] for interval in intervals] + [ interval[1] + 1 for interval in intervals] 
    bins_starts = list(set(bins_starts))
    bins_starts.sort()
    
    bins_lengths = [ end - start for start, end in zip( bins_starts, bins_starts[1:])]
    
    first_interval_unbiased_freqs_list = generate_unbiased_freqs_list(bins_lengths[:int_num], int_thresholds[0])
    last_interval_unbiased_freqs_list = generate_unbiased_freqs_list(bins_lengths[int_num - 1:], int_thresholds[-1])

    prob_intervals = []

    exp_value = 0

    for first_interval_unbiased_freqs in first_interval_unbiased_freqs_list:
        for last_interval_unbiased_freqs in last_interval_unbiased_freqs_list:
            broken_flag = False
            
            if first_interval_unbiased_freqs[-1] == last_interval_unbiased_freqs[0]:
                unbiased_freqs = first_interval_unbiased_freqs + last_interval_unbiased_freqs[1:]
 
                for i, int_to_check in enumerate(intervals):
                    if sum(unbiased_freqs[i:i + int_num]) > int_thresholds[i]:
                        broken_flag = True
                        break
                
                if broken_flag: break

                prob = 1
                
                for bin_length, unbiased_freq in zip(bins_lengths, unbiased_freqs): 
                    prob *= comb(bin_length, unbiased_freq, exact = True) * ((1 - p) ** unbiased_freq) * p ** (bin_length - unbiased_freq) 

                exp_value += bins_lengths[len(bins_lengths) // 2] * prob

    return exp_value

def get_expected_value_BCS_naively(intervals, p):
    intervals.sort()

    start = intervals[0][0]
    end = intervals[-1][1]
    
    intervals = [ (i1 - start, i2 - start) for i1, i2 in intervals] 
    int_length = end - start + 1

    outcomes_list = list(product([0,1], repeat = int_length))
    
    exp_value = 0

    for outcome in outcomes_list:
        prob = (1 - p) ** sum(outcome) * p ** (int_length - sum(outcome))
        BCSes = set()
        for interval in intervals:
            if sum(outcome[interval[0]:interval[1] + 1]) <= math.floor(STW_MAX_PERCENT * (interval[1] - interval[0] + 1)):
                BCSes = BCSes.union(list(range(interval[0], interval[1] + 1)))
            exp_value += prob * len(BCSes)
    
    return exp_value
