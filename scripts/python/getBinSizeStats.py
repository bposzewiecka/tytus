derived_in_type = snakemake.wildcards['derived']
chrom = snakemake.wildcards['chromosome']

liftover_fasta_file = snakemake.input['liftover_fasta_file']
stats_tsv_file = snakemake.output['stats_tsv_file']

from SND import get_SNDs_derived_in_target
from SND import get_SNDs_derived_in_query

from collections import defaultdict
from consts import CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH 

from SNDWindow import SNDWindow

NUMBER_OF_BINS_LIST = 1, 2, 3, 4, 5, 10, 300

def get_bin_size_statistics():
    
    with open(stats_tsv_file, 'w') as f_out:
    
        if derived_in_type == 'target':
            snds = get_SNDs_derived_in_target(liftover_fasta_file)
            for snd in snds:
                snd.biased = snd.biased_in_target()
        else:
            snds = get_SNDs_derived_in_query(liftover_fasta_file)
            for snd in snds:
                snd.biased = snd.biased_in_query()

        for number_of_bins in NUMBER_OF_BINS_LIST:

            number_of_clustered = defaultdict(int)
            number_of_biased_clustered = defaultdict(int)
                
            for index in range(len(snds)):
                snd_window = SNDWindow(index, snds, CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH, number_of_bins)
                chrom_bin = snd_window.get_middle_coord() // 1000 // 1000
                is_clustered = snd_window.is_clustered() 
                is_biased_clustered = snd_window.is_biased_clustered() 

                if is_clustered:
                    number_of_clustered[chrom_bin] += 1

                if is_biased_clustered:
                    number_of_biased_clustered[chrom_bin] += 1  


            for region in sorted(number_of_clustered):
                o = chrom, region * 1000 * 1000, (region + 1) * 1000 * 1000 - 1, number_of_bins, number_of_clustered[region],  number_of_biased_clustered[region], number_of_biased_clustered[region] / number_of_clustered[region] 
                f_out.write('\t'.join(map(str, o)))
                f_out.write('\n')


get_bin_size_statistics()
