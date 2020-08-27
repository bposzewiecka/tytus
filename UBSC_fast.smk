apes = [ 'ponAbe3', 'panPan2', 'panTro5', 'nomLeu3', 'gorGor5'] # papAnu4
chroms =  [ 'chr' + str(name) for name in list(range(1,23)) + ['X', 'Y']]


rule main:
    input:	
        expand('data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.w{window}.derived_in_{derived}.outgroup.{outgroup}.bin_size_stats.tsv'
, target = 'hg38', query = 'panTro5', window = 5, outgroup = 'rheMac8',  derived = 'target', type= 'target', chromosome = chroms),
	#expand('data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.w{window}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.compressed_bins_stats.tsv', target = 'hg38', query = apes, window = 5, outgroup = 'rheMac8',  derived = 'target', type= 'target', chromosome = chroms, window_size = 300),
	#expand('data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}/{type}.{target}.{query}.{chromosome}.w{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.rbest.snd.liftover.ubcs.tsv', target = 'hg38', query = ['panTro5', 'panPan2', 'gorGor5'], window = 5, outgroup = 'rheMac8',  derived = 'target', type= 'target', chromosome = chroms, window_size = 300, region = 1000 * 1000)


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


rule generate_compressed_bins_statistics:
    input:
        liftover_fasta_file = 'data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.w{window}.outgroup.{outgroup}.rbest.snd.liftover.fasta'
    output:
        stats_tsv_file = 'data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.{chromosome}.w{window}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.compressed_bins_stats.tsv'
    script:
        'scripts/python/getCompressedBinsStats.py'


rule join_compressed_bins_statistics:
    input:
        stats_tsv_file = lambda wildcards: expand('data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.{chromosome}.w{window}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.compressed_bins_stats.tsv', chromosome = chroms,  **wildcards)
    output:
        stats_tsv_file = 'data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.w{window}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.compressed_bins_stats.tsv'
    shell:
        'cat {input} > {output}'


rule generate_ubcs_statistics:
    input:
        liftover_fasta_file = 'data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.w{window}.outgroup.{outgroup}.rbest.snd.liftover.fasta'
    output:
        ubcs_stats_file = 'data/{query}/{target}.{query}/{type}/w{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}/{type}.{target}.{query}.{chromosome}.w{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.rbest.snd.liftover.ubcs.tsv'
    script:
        'scripts/python/getUBCSFastStats.py'

	
