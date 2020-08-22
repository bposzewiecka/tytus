
from collections import Counter
from collections import defaultdict
from SND import get_SNDs_derived_in_target
from SND import get_SNDs_derived_in_query
from consts import CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH
from consts import CLUSTERED_SUBSTITUTIONS_MIN_SIZE
from consts import WTS_MIN_PERCENT
from UBCS_fast import UBCS_fast
from UBCS_fast import compress_snd_window
from UBCS_fast import get_snds_in_window
from UBCS_fast import get_snd_window
from UBCS_fast import contains_cluster
from UBCS_fast import is_biased_clustered

liftover_fasta_file = snakemake.input['liftover_fasta_file']
ubcs_stats_file = snakemake.output['ubcs_stats_file']
derived_in_type = snakemake.wildcards['derived']
ttype = snakemake.wildcards['type']

REGION_LEN = int(snakemake.wildcards['region'])
chrom = snakemake.wildcards['chromosome']
WINDOW_SIZE =  int(snakemake.wildcards['window_size'])
target = snakemake.wildcards['target']
query = snakemake.wildcards['query']
outgroup = snakemake.wildcards['outgroup']


def get_chrom_SNDs(liftover_fasta_file):

    if derived_in_type == 'target':
        SNDs = get_SNDs_derived_in_target(liftover_fasta_file)
    else:
        SNDs = get_SNDs_derived_in_query(liftover_fasta_file)

    for SND in SNDs:
        if ttype == 'target':
            SND.chrom = SND.target_chrom
            SND.start = SND.target_start
            SND.end = SND.target_end
            SND.coord = SND.target_coord()
        else:
            SND.chrom = SND.query_chrom
            SND.start = SND.query_start
            SND.end = SND.query_end
            SND.coord = SND.query_coord()

        if derived_in_type == 'target':
            SND.biased = SND.biased_in_target()
        else:
            SND.biased = SND.biased_in_query()


    return SNDs


def get_UBCS_stats():

    with  open(ubcs_stats_file, 'w') as f_out:

        SNDs = get_chrom_SNDs(liftover_fasta_file)
        region_snds_number = defaultdict(int)
        region_biased_snds_number = defaultdict(int)
        region_p = defaultdict(int)
        region_expected_BCS = defaultdict(int)
        region_actual_BCS = defaultdict(int)

        for SND in SNDs:

            region = SND.coord // REGION_LEN
            region_snds_number[region] += 1

            if SND.biased:
                region_biased_snds_number[region] += 1

        for region in  region_snds_number:
            region_p[region] = region_biased_snds_number[region] / region_snds_number[region]

        for index, SND in enumerate(SNDs):    
           
            region = SND.coord // REGION_LEN
            window = get_snd_window(index, SNDs, WINDOW_SIZE)

            if contains_cluster(window):

                compressed_snd_window = compress_snd_window(window)
                if len(compressed_snd_window) > 65:
                    continue
                region_actual_BCS[region] += int(is_biased_clustered(window))
                region_expected_BCS[region] += UBCS_fast(compressed_snd_window, region_p[region])


        for region in region_p:
            p = region_p[region]
            expected_BCS = region_expected_BCS[region] 
            actual_BCS =  region_actual_BCS[region]
            o = chrom, region * REGION_LEN, (region + 1) * REGION_LEN -1, ttype, derived_in_type, target, query, outgroup, p, expected_BCS, actual_BCS, actual_BCS - expected_BCS
            f_out.write('\t'.join(map(str, o)))
            f_out.write('\n')

get_UBCS_stats()
