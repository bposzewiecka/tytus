import gzip

## CONSTANTS TO BE SET
chr_number =  snakemake.wildcards['chromosome']
query_name = snakemake.wildcards['query']
target_name = snakemake.wildcards['target']
window = int(snakemake.wildcards['window'])
ttype = snakemake.wildcards['type']


## compressed VCF filename to be analyzed
vcf_file = snakemake.input['vcf_file']
## BED output filename 
snd_file = snakemake.output['snd_file']

def get_ref_alt(pos_info, pos):

    ref = ''
    alt = ''

    if pos_info[pos][1] == '.': return None, None

    diffs = 0
    
    for coord in range(pos - window, pos + window + 1):

        if coord in pos_info:

           ref_base, alt_base = pos_info[coord][0], pos_info[coord][1]

           if len(ref_base) != 1 or len(alt_base) != 1: return None, None

           ref = ref + ref_base
           alt = alt + (alt_base if alt_base != '.' else ref_base)

           if  alt_base != '.': diffs += 1

        else:
           return None, None

        if diffs - 1 > 2: return None, None

    return ref, alt

def print_snd_bed_entry(f, pos_info, pos, snp_count):
    
    ref, alt = get_ref_alt(pos_info, pos)

    if ref and alt:
        snp_count   = snp_count + 1

        start_pos   = str(pos - window - 1)
        end_pos     = str(pos + window)

        chr_name    = 'chr' + chr_number
        snd_id      = 'SND_ID:' + str(snp_count) + ';'
        query       = 'QUERY:' + query_name + ';'
        target      = 'TARGET:' + target_name  + ';'
        coords      = 'chr' + chr_number + ':' + start_pos + '-'+ end_pos + ';'
        target_crd  = 'TARGET_COORDS:' + coords
        query_crd   = 'QUERY_COORDS:' + coords
        change      = 'CHANGE:' + ref + '>' + alt + ';'
        info        = 'TYPE:'+ ttype  + ';'

        line_out = chr_name + '\t' + start_pos + '\t' + end_pos + '\t' + snd_id + query + target + target_crd + query_crd + change + info + '\t0\t+\n'

        f.write(line_out)

        return 1
    else:
        return 0


def read_pos_info():

    pos_info = {}

    i = 0

    with gzip.open(vcf_file, 'rt') as f:

        for line in f:

            # skipping header
            if line.startswith('#'): continue

            line_tab =  line.strip().split('\t')

            pos_info[int(line_tab[1])] = line_tab[3], line_tab[4] 

            i += 1

    return pos_info


def save_snd_info(pos_info):

    with open(snd_file, 'w') as f:

        snp_count = 1

        for pos in sorted(pos_info.keys()):
            snp_count += print_snd_bed_entry(f, pos_info, pos, snp_count)


pos_info = read_pos_info()
save_snd_info(pos_info)
