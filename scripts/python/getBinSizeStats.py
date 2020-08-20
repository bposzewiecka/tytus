BINS = 1, 2, 3, 4, 5, 10, 300

derived_in_type = snakemake.wildcards['derived']
chrom = snakemake.wildcards['chromosome']

liftover_fasta_file = snakemake.input['liftover_fasta_file']
stats_tsv_file = snakemake.output['stats_tsv_file']

from SND import get_SNDs_derived_in_target
from collections import defaultdict
from collections import Counter
from consts import CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH 

def get_clustered_substitutions_coords(snds, number_of_bins):
    coords = set()

    bin_size = CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH  / number_of_bins
    previous = -CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH 

    n = len(snds)

    for shift_no in range(number_of_bins):
        
        i = 0
        window_no = -1
        substs_in_window = []
        previous_size =  0
     
        while i < n:
            target_coord = snds[i].target_coord()
            curr_window = (target_coord + shift_no * bin_size) // CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH 
            
            if curr_window  == window_no:
                substs_in_window.append(target_coord)
                
            else:
                if len(substs_in_window) > 4:
                    coords.update(substs_in_window)

                substs_in_window = [ target_coord ]
                window_no = curr_window 
            i += 1  
            
    return coords

def get_bin_size_statistics(number_of_bins_list):
    
    with open(stats_tsv_file, 'w') as f_out:
    
        if derived_in_type == 'target':
            snds = get_SNDs_derived_in_target(liftover_fasta_file)
        else:
            snds = get_SNDs_derived_in_query(liftover_fasta_file)

        for number_of_bins in number_of_bins_list:
            coords = get_clustered_substitutions_coords(snds, number_of_bins)

            counter = Counter([ coord // 1000 // 1000 for coord in sorted(coords)])

            for megabase, value in sorted(counter.items()):
                o = chrom, megabase * 1000 * 1000, (megabase + 1) * 1000 * 1000, number_of_bins, value
                f_out.write('\t'.join(map(str, o)))
                f_out.write('\n')

get_bin_size_statistics(BINS)
