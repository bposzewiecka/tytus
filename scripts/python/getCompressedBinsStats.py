WINDOW_SIZE =  snakemake.wildcards['window_size']

derived_in_type = snakemake.wildcards['derived']
chrom = snakemake.wildcards['chromosome']

liftover_fasta_file = snakemake.input['liftover_fasta_file']
stats_tsv_file = snakemake.output['stats_tsv_file']

from SND import get_SNDs_derived_in_target
from SND import get_SNDs_derived_in_query
from collections import Counter
from consts import CLUSTERED_SUBSTITUTIONS_MIN_SIZE

def compress_substitutions_occurences(substs):
    
    n = len(substs) // 2 + 1

    starts = []
    indexes = set()
    
    for i in range(n):
        
        if sum(substs[i:i + n]) >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
            
            ind = tuple( i + j for j, v in enumerate(substs[i:i + n]) if v == 1)
            
            if ind not in indexes:
                
                indexes.add(ind)
                starts.append(i)
                
    starts = starts + [ start + n for start in starts ]

    return  [ sum(substs[start: end]) for start, end in zip(starts, starts[1:])]

def get_substitutions_occurrences(coord_index, coords):
    start =  end = coord_index
    middle_coord = coords[coord_index]
    
    while start >= 0 and middle_coord - coords[start] < WINDOW_SIZE:
        start -= 1
   
    while end < len(coords) and coords[end] - middle_coord < WINDOW_SIZE:
        end += 1
    
    if end - start - 1 < CLUSTERED_SUBSTITUTIONS_MIN_SIZE:
        return []
    
    occurrences = [0] * (2 * WINDOW_SIZE - 1)
    
    for coord in coords[start + 1: end]:
        occurrences[coord - middle_coord + WINDOW_SIZE - 1] = 1
    
    return occurrences
    

def get_compressed_bins_stats():
    
    
    with open(stats_tsv_file, 'w') as f_out:

        if derived_in_type == 'target':
            snds = get_SNDs_derived_in_target(liftover_fasta_file)
        else:
            snds = get_SNDs_derived_in_query(liftover_fasta_file)
    
        coords = [ snd.target_coord() for snd in snds]
        
        compressed_bins_sizes = []

        for coord_index in range(len(coords)):
            occurrences = get_substitutions_occurrences(coord_index, coords)
            
            compressed_bins = compress_substitutions_occurences(occurrences)
            compressed_bins_sizes.append(len(compressed_bins))
        
        counter = Counter(compressed_bins_sizes)

        for compressed_bins_size, value in sorted(counter.items()):
            o = chrom, compressed_bins_size,  value
            f_out.write('\t'.join(map(str, o)))
            f_out.write('\n')

get_compressed_bins_stats()    
