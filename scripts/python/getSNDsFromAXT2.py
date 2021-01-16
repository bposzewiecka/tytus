from utils import get_chrom_sizes
from consts import QUALITY_WINDOW_LENGTH

from Bio.Seq import Seq

alignment_file = snakemake.input['alignment_file']
snd_file = snakemake.output['snd_file']
chrom_sizes_file = snakemake.input['query_chrom_sizes_file']

chrom_sizes = get_chrom_sizes(chrom_sizes_file)

query_name = snakemake.wildcards['query']
target_name = snakemake.wildcards['target']
ffrom = snakemake.params['ffrom']

snp_file = snakemake.input['snp_file']

def get_snps(snp_file):
   
    snps = {}

    with open(snp_file) as f:

        for line in f:

            line = line.split()

            start_coord = int(line[1]) + 1
            end_coord = int(line[2]) + 1

            for i in range(start_coord, end_coord):
                base = line[4]
                snps[i] = (base, line)

    return snps


def get_axt_entry(f, ffrom):

    header = f.readline()

    if not header: return None

    header = header.split()

    strand = header[7]

    if ffrom == 'target':

        t_chrom = header[1]
        t_start = int(header[2])

        q_chrom = header[4]
        q_start = int(header[5])
        q_end = int(header[6])

        target_al_seq = f.readline().rstrip().upper()
        query_al_seq = f.readline().rstrip().upper()

    else:

        t_chrom = header[4]
        t_start = int(header[5])

        q_chrom = header[1]
        q_start = int(header[2])

        query_al_seq = f.readline().rstrip().upper()
        target_al_seq =  f.readline().rstrip().upper()           

    f.readline()

    return t_chrom, t_start, q_chrom, q_start, strand, target_al_seq, query_al_seq 

def add_snds_from_axt_entry(t_chrom, t_start, q_chrom, q_start, strand, target_al_seq, query_al_seq, snd_bed_entries):

    t_count = 0
    q_count = 0

    for i, (t_chr, q_chr) in enumerate(zip(target_al_seq, query_al_seq)):

        if t_chr != '-':
            t_count += 1

        if q_chr != '-':
            q_count += 1

        if t_chr == q_chr:
            continue

        t_window = target_al_seq[i - QUALITY_WINDOW_LENGTH: i + QUALITY_WINDOW_LENGTH + 1]
        q_window = query_al_seq[i - QUALITY_WINDOW_LENGTH: i + QUALITY_WINDOW_LENGTH + 1]

        if len(t_window) != 2 * QUALITY_WINDOW_LENGTH + 1:  continue

        #print('enough place')

        if '-' in t_window + q_window: continue

        #print('without deletions and insertions')

        if len([ 1  for t_win_chr, q_win_chr  in zip(t_window, q_window) if t_win_chr != q_win_chr]) > 2: continue

        #print('number of  mismatches below

        t_subst_start = t_start + t_count - 2

        if strand == '-' and ffrom == 'query':
            t_subst_start = chrom_sizes[t_chrom] - (t_start + t_count - 1)
            t_window = str(Seq(t_window).reverse_complement())
            q_window = str(Seq(q_window).reverse_complement())

        t_coord_start = str(t_subst_start - QUALITY_WINDOW_LENGTH)
        t_coord_end =  str(t_subst_start + QUALITY_WINDOW_LENGTH + 1)

        t_coords = t_chrom + ':' + t_coord_start + '-' + t_coord_end

        q_subst_start = q_start + q_count - 2

        if strand == '-' and ffrom == 'target':
            q_subst_start = chrom_sizes[q_chrom] - (q_start + q_count - 1)

        q_coord_start = str(q_subst_start - QUALITY_WINDOW_LENGTH)
        q_coord_end = str(q_subst_start + QUALITY_WINDOW_LENGTH + 1)

        q_coords = q_chrom + ':' + q_coord_start + '-' + q_coord_end

        t_snp =  int(t_coord_start) +  QUALITY_WINDOW_LENGTH + 1 in snps

        #if t_output_window[5] != snps[int(t_coord_start)+ 6][0]:
        #    print(t_output_window, t_output_window[5], snps[int(t_coord_start) + 6]

        name = ('SND_ID:'  + str(len(snd_bed_entries)) , 'QUERY:' + query_name, 'TARGET:'+ target_name,
                'TARGET_COORDS:' + t_coords, 'QUERY_COORDS:' + q_coords,
                'CHANGE:' +  t_window + '>' + q_window, 'FROM:' + ffrom, 'SNP:' + str(t_snp))

        name = ';'.join(name)

        snd_bed_entry = '\t'.join([t_chrom, t_coord_start, t_coord_end, name, '0', strand ])

        snd_bed_entries[int(t_coord_start)] = snd_bed_entry

def get_snd_bed_entries(alignment_file, snps):

    with open(alignment_file) as f:

        snd_bed_entries = {}

        while True:

            axt_entry = get_axt_entry(f, ffrom)
         
            if axt_entry is None:
                break

            t_chrom, t_start, q_chrom, q_start, strand, target_al_seq, query_al_seq = axt_entry
            
            add_snds_from_axt_entry(t_chrom, t_start, q_chrom, q_start, strand, target_al_seq, query_al_seq, snd_bed_entries)
	
        return snd_bed_entries

def save_snd_bed_entries(snd_bed_entries_file, snd_bed_entries):

    with open(snd_file, 'w') as f_out:
 
        for t_coord_start in sorted(snd_bed_entries):
            snd_bed_entry = snd_bed_entries[t_coord_start]
            f_out.write(snd_bed_entry)
            f_out.write('\n')

snps = get_snps(snp_file)
snd_bed_entries = get_snd_bed_entries(alignment_file, snps)
save_snd_bed_entries(snd_file, snd_bed_entries)
