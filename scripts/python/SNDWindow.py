from collections import defaultdict
from consts import CLUSTERED_SUBSTITUTIONS_MIN_SIZE
from consts import WTS_MIN_PERCENT


class SNDWindow:

    def __init__(self, index, snds, window_size, number_of_bins):
        self.window_size = window_size
        self.number_of_bins = number_of_bins
        self.bin_size = window_size // self.number_of_bins
        self.index = index
        self.window = defaultdict(list)
        self.window_counts = [0] * (2 * number_of_bins - 1)
        self.window_biased_counts = [0] * (2 * number_of_bins - 1)
    
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
            
            if snd.biased:
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

