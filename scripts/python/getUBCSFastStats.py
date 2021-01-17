from collections import defaultdict
from SND import get_SNDs_derived_in_target
from SND import get_SNDs_derived_in_query
from SNDWindow import SNDWindow

liftover_fasta_file = snakemake.input['liftover_fasta_file']
ubcs_stats_file = snakemake.output['ubcs_stats_file']
ubcs_stats_details_file = snakemake.output['ubcs_stats_details_file']
derived_in_type = snakemake.wildcards['derived']
query = snakemake.wildcards['query']
target = snakemake.wildcards['target']
outgroup = snakemake.wildcards['outgroup']


REGION_LEN = int(snakemake.wildcards['region'])
chrom = snakemake.wildcards['chromosome']
WINDOW_SIZE =  int(snakemake.wildcards['window_size'])
NUMBER_OF_BINS = int(snakemake.wildcards['number_of_bins'])

def get_UBCS_stats():

    with  open(ubcs_stats_file, 'w') as f_out, open(ubcs_stats_details_file, 'w') as f_details_out:        

        region_snds_number = defaultdict(int)
        region_biased_snds_number = defaultdict(int)
        region_p = defaultdict(int)
        region_expected_BCS = defaultdict(int)
        region_actual_BCS = defaultdict(int)

        if derived_in_type == 'target':
            SNDs = get_SNDs_derived_in_target(liftover_fasta_file)
        else:
            SNDs = get_SNDs_derived_in_query(liftover_fasta_file)

        for SND in SNDs:

            region = SND.target_coord() // REGION_LEN
            region_snds_number[region] += 1

            if SND.biased():
                region_biased_snds_number[region] += 1

        for region in  region_snds_number:
            region_p[region] = region_biased_snds_number[region] / region_snds_number[region]

        for index, SND in enumerate(SNDs):
            coord = SND.target_coord()
         
            region = coord // REGION_LEN
            window = SNDWindow(index, snds = SNDs , window_size = WINDOW_SIZE, number_of_bins = NUMBER_OF_BINS)

            if WINDOW_SIZE == NUMBER_OF_BINS and window.get_max_frequency_of_cluster() >= 23:
                window = SNDWindow(index, snds = SNDs , window_size = WINDOW_SIZE, number_of_bins = 20)

            is_biased = SND.biased()
            is_clustered = window.is_clustered()
            is_biased_clustered = window.is_biased_clustered()
            
            p =  region_p[region]
            prob_of_bcs = window.get_prob_of_bcs(p)
            region_actual_BCS[region] += int(is_biased_clustered)
            region_expected_BCS[region] += prob_of_bcs
            o = chrom,  REGION_LEN, NUMBER_OF_BINS, ttype, derived_in_type, target, query, outgroup, coord, p, is_biased, is_clustered, is_biased_clustered,  prob_of_bcs
            f_details_out.write('\t'.join(map(str, o)))
            f_details_out.write('\n')


        for region in region_p:
            p = region_p[region]
            expected_BCS = region_expected_BCS[region]
            actual_BCS =  region_actual_BCS[region]
            o = chrom, region * REGION_LEN, (region + 1) * REGION_LEN -1, WINDOW_SIZE, NUMBER_OF_BINS, ttype, derived_in_type, target, query, outgroup, p, expected_BCS, actual_BCS, actual_BCS - expected_BCS
            f_out.write('\t'.join(map(str, o)))
            f_out.write('\n')

get_UBCS_stats()               
