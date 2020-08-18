from collections import Counter
from SND import get_SNDs_derived_in_target
from SND import get_SNDs_derived_in_query
from SND import get_SNDs
from UBCS import get_expected_BCS
from consts import CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH 
from consts import CLUSTERED_SUBSTITUTIONS_MIN_SIZE
from consts import WTS_MIN_PERCENT 

SHIFT_LENGTH = CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH // 2

liftover_fasta_file = snakemake.input['liftover_fasta_file']
ubcs_stats_file = snakemake.output['ubcs_stats_file']
derived_in_type = snakemake.wildcards['derived']
ttype = snakemake.wildcards['type']

REGION_LEN = int(snakemake.wildcards['region'])

def get_chrom_SNDs(liftover_fasta_file):
    
    chrom_SNDs = {}

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

        if SND.chrom in chrom_SNDs:
            chrom_SNDs[SND.chrom].append(SND)
        else:
            chrom_SNDs[SND.chrom] = [SND]

    return chrom_SNDs


def get_SNDs_clusters(liftover_fasta_file):
    chrom_SNDs = get_chrom_SNDs(liftover_fasta_file)

    SNDs_clusters = {}

    for chrom in chrom_SNDs:
        SNDs = chrom_SNDs[chrom]
        SNDs.sort(key = lambda x: x.coord)

        n = len(SNDs)

        clusters = {}

        for i, snd in enumerate(SNDs):
            for shift in [0, SHIFT_LENGTH]:

                cluster_no = (snd.coord + shift) // SHIFT_LENGTH

                if cluster_no in clusters:
                     clusters[cluster_no][1] = i
                else:
                     clusters[cluster_no] = [i, i]

        clusters = set(tuple(cluster) for cluster in clusters.values())
        clusters = [ (start, end)  for start, end in clusters if end - start + 1 >= CLUSTERED_SUBSTITUTIONS_MIN_SIZE]

        SNDs_clusters[chrom] = clusters, SNDs
        
    return SNDs_clusters


#def get_SNDs_clusters(liftover_fasta_file):
#    chrom_SNDs = get_chrom_SNDs(liftover_fasta_file)
#
#    print('SND' + str(region_no))
#
#    SNDs_clusters = {}
#
#    for chrom in chrom_SNDs:
#        SNDs = chrom_SNDs[chrom]
#        SNDs.sort(key = lambda x: (x.start, x.end))
#
#        n = len(SNDs)
#
#        clusters = []
#
#        for i, snd in enumerate(SNDs):
#            
#            j = i
#
#            while j < n and snd.start + CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH  > SNDs[j].start:
#                j += 1
#
#            cluster_len = j - i
#            if cluster_len < CLUSTERED_SUBSTITUTIONS_MIN_SIZE: continue
#
#            if not clusters or clusters[-1][1] != j - 1:
#                clusters.append((i, j - 1))
#
#        SNDs_clusters[chrom] = (clusters, SNDs)
#
#    return SNDs_clusters

def get_actual_BCS(clusters, SNDs):
    BCS = set()

    for cluster in clusters:
        
        BCS_in_cluster =  sum([ 1 for i in range(cluster[0], cluster[1] + 1) if SNDs[i].biased ])
        
        if BCS_in_cluster >= WTS_MIN_PERCENT * (cluster[1] - cluster[0] + 1):
            BCS.update(list(range(cluster[0], cluster[1] + 1)))

    return len(BCS)


def get_UBCS_stats(liftover_fasta_file, ubcs_stats_file):

    with  open(ubcs_stats_file, 'w') as f_out:
    
        for chrom, (clusters, SNDs) in get_SNDs_clusters(liftover_fasta_file).items():
            last_coord = SNDs[-1].end 
        
            regions = { i : { 'clusters':  [] , 'SNDs': [] }  for i in range(last_coord // REGION_LEN + 1) }

            for cluster in clusters:
                end_coord = SNDs[cluster[1]].end
                regions[end_coord // REGION_LEN]['clusters'].append(cluster)

            for SND in  SNDs:
                end_coord = SND.end
                regions[end_coord // REGION_LEN]['SNDs'].append(SND)

            for region_no in regions:

                region = regions[region_no]

                clusters = region['clusters']
                clusterSNDs = region['SNDs']

                print(region_no, clusters)

                if not clusters: continue
            
                start_SND = clusters[0][0]
                end_SND = clusters[-1][1]

                p = sum([ 1  for snd in clusterSNDs if snd.biased]) / len(clusterSNDs)

                expected_BCS = int(get_expected_BCS(clusters, p))

                actual_BCS = get_actual_BCS(clusters, SNDs)

                o = chrom, region_no * REGION_LEN, (region_no + 1) * REGION_LEN -1, ttype, derived_in_type, snakemake.wildcards['target'], snakemake.wildcards['query'], snakemake.wildcards['outgroup'], p, expected_BCS, actual_BCS, actual_BCS - expected_BCS
                f_out.write('\t'.join(map(str, o)))
                f_out.write('\n')
        
    
get_UBCS_stats(liftover_fasta_file, ubcs_stats_file)
