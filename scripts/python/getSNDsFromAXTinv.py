from utils import get_chrom_sizes
from consts import QUALITY_WINDOW_LENGTH 

from Bio.Seq import Seq

alignment_file = snakemake.input['alignment_file']
snd_file = snakemake.output['snd_file']
chrom_sizes_file = snakemake.input['query_chrom_sizes_file']
chrom_sizes = get_chrom_sizes(chrom_sizes_file)

query_name = snakemake.wildcards['query']
target_name = snakemake.wildcards['target']
output_window_len = int(snakemake.wildcards['window'])
ttype = snakemake.wildcards['type']

snp_file = snakemake.input['snp_file']

snps = {}

with open(snp_file) as f:
    
    for line in f:
        
        line = line.split()

        start_coord = int(line[1]) + 1
        end_coord = int(line[2]) + 1

        for i in range(start_coord, end_coord):
            base = line[4]
            snps[i] = (base, line)


with open(alignment_file) as f_in:
    
    snd_id, snp_id = 1, 1
    lines = {}
 
    while True:
        header = f_in.readline()
        if not header: break

        header = header.split()
        t_chrom = header[4]
        t_start = int(header[5])

        q_strand = header[7]
        q_chrom = header[1] 
        q_start = int(header[2])

        query_al_seq= f_in.readline().rstrip().upper()
        target_al_seq =  f_in.readline().rstrip().upper()

        t_count = 0
        q_count = 0

        for i, (t_chr, q_chr) in enumerate(zip(target_al_seq, query_al_seq)):

            if t_chr != '-':
                t_count += 1

            if q_chr != '-':
                q_count += 1

            if  t_chr == q_chr:
                continue

            t_window = target_al_seq[i - QUALITY_WINDOW_LENGTH: i + QUALITY_WINDOW_LENGTH + 1]
            q_window = query_al_seq[i - QUALITY_WINDOW_LENGTH: i + QUALITY_WINDOW_LENGTH + 1]

            if len(t_window) != 2 * QUALITY_WINDOW_LENGTH + 1:  continue

            #print('wystarczajaco miejsca')

            if '-' in t_window + q_window: continue

            #print('bez myslnika')

            if len([ 1  for t_win_chr, q_win_chr  in zip(t_window, q_window) if t_win_chr != q_win_chr]) > 2: continue

            #print('dozwolona ilosc bledow')

            t_output_window = target_al_seq[i - output_window_len: i + output_window_len + 1]
            q_output_window = query_al_seq[i - output_window_len: i + output_window_len + 1]

            if q_strand == '-':
                t_output_window =  str(Seq(t_output_window).reverse_complement())
                q_output_window =  str(Seq(q_output_window).reverse_complement())

            if q_strand == '+':
                t_subst_start = t_start + t_count - 2
            else:
                t_subst_start = chrom_sizes[t_chrom] - (t_start + t_count - 1)

            t_coord_start = str(t_subst_start - output_window_len)
            t_coord_end =  str(t_subst_start + output_window_len + 1)

            t_coords = t_chrom + ':' + t_coord_start + '-' + t_coord_end

            if int(t_coord_start) +  QUALITY_WINDOW_LENGTH + 1 in snps:
                snp_id += 1
                #if t_output_window[5] != snps[int(t_coord_start)+ 6][0]:
                #    print(t_output_window, t_output_window[5], snps[int(t_coord_start) + 6])
                continue

            q_subst_start = q_start + q_count - 2
            q_coord_start = str(q_subst_start - output_window_len)
            q_coord_end = str(q_subst_start + output_window_len + 1)

            q_coords = q_chrom + ':' + q_coord_start + '-' + q_coord_end
            
            name = ('SND_ID:'  + str(snd_id) , 'QUERY:' + query_name, 'TARGET:'+ target_name, 
                    'TARGET_COORDS:' + t_coords, 'QUERY_COORDS:' + q_coords, 
                    'CHANGE:' +  t_output_window + '>' + q_output_window, 'TYPE:' + ttype)

            name = ';'.join(name)

            line = '\t'.join([t_chrom, t_coord_start, t_coord_end, name, '0', q_strand ])

            lines[int(t_coord_start)] = line

            snd_id += 1

        f_in.readline()

with open(snd_file, 'w') as f_out:

    for t_coord_start in sorted(lines):
        line = lines[t_coord_start]
        f_out.write(line)
        f_out.write('\n')

print(snd_id, snp_id)