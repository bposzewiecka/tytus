apes = [ 'panTro6',  'papAnu4', 'ponAbe3', 'panPan2', 'panTro5', 'nomLeu3', 'gorGor5', 'gorGor4'] 
chroms =  [ 'chr' + str(name) for name in list(range(1,23)) + ['X', 'Y']]


rule main:
    input:	
        expand('data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.w{window}.derived_in_{derived}.outgroup.{outgroup}.bin_size_stats.tsv'
, target = 'hg38', query = 'panTro5', window = 5, outgroup = 'rheMac8',  derived = 'target', type= 'target', chromosome = chroms),


rule generate_bin_statistics:
    input:
        liftover_fasta_file = 'data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.w{window}.outgroup.{outgroup}.rbest.snd.liftover.fasta'
    output:
        stats_tsv_file = 'data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.{chromosome}.w{window}.derived_in_{derived}.outgroup.{outgroup}.bin_size_stats.tsv'
    script:
        'scripts/python/getBinSizeStats.py'


rule join_bin_statistics:
    input:
        stats_tsv_file = lambda wildcards: expand('data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.{chromosome}.w{window}.derived_in_{derived}.outgroup.{outgroup}.bin_size_stats.tsv', chromosome = chroms,  **wildcards)
    output:
        stats_tsv_file = 'data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.w{window}.derived_in_{derived}.outgroup.{outgroup}.bin_size_stats.tsv'
    shell:
        'cat {input} > {output}'	 
	
