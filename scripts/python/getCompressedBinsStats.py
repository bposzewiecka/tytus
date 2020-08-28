WINDOW_SIZE =  int(snakemake.wildcards['window_size'])

derived_in_type = snakemake.wildcards['derived']
chrom = snakemake.wildcards['chromosome']

liftover_fasta_file = snakemake.input['liftover_fasta_file']
stats_tsv_file = snakemake.output['stats_tsv_file']

from SND import get_SNDs_derived_in_target
from SND import get_SNDs_derived_in_query
from collections import Counter
from SNDWindow import SNDWindow

def get_compressed_bins_stats():

    with open(stats_tsv_file, 'w') as f_out:

        if derived_in_type == 'target':
            snds = get_SNDs_derived_in_target(liftover_fasta_file)
        else:
            snds = get_SNDs_derived_in_query(liftover_fasta_file)

        compressed_window_stats = []

        for index in range(len(snds)):
            snd_window = SNDWindow(index, snds,  WINDOW_SIZE, WINDOW_SIZE)
            compresed_window = snd_window.get_compressed_window_counts()
            
            compressed_window_stats.append((len(compresed_window), snd_window.get_max_frequency_of_cluster()))
        
        counter = Counter(compressed_window_stats)
        
        for (window_len, max_window_freq), value in sorted(counter.items()):
            o = chrom, window_len, max_window_freq, value
            f_out.write('\t'.join(map(str, o)))
            f_out.write('\n')
            
get_compressed_bins_stats()     

