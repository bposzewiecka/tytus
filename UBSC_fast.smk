apes = [ 'ponAbe3', 'panPan2', 'panTro5', 'nomLeu3', 'gorGor5'] # papAnu4
chroms =  [ 'chr' + str(name) for name in list(range(1,23)) + ['X', 'Y']]

rule main:
    input:
        expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}.window_size.{window_size}.number_of_bins.{number_of_bins}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.number_of_bins.{number_of_bins}.rbest.from.{axt_dir}.snd.liftover.ubcs.summary.tsv', target = 'hg38', query = ['panTro6', 'panPan3'], window = 5, outgroup = 'rheMac10',  derived = ['target', 'query'], type= 'target', chromosome = chroms, window_size = 300, number_of_bins = [300, 2], region = 1000000, axt_dir = 'query'),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.v{window}.derived_in_{derived}.outgroup.{outgroup}.bin_size_stats.tsv', target = 'hg38', query = 'panTro5', window = 5, outgroup = 'rheMac8',  derived = 'target', type= 'target', chromosome = chroms),
	#expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.v{window}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.compressed_bins_stats.tsv', target = 'hg38', query = 'panTro5', window = 5, outgroup = 'rheMac8',  derived = ['target', 'query'], type= 'target', chromosome = chroms, window_size = [300, 240, 180]),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}.window_size.{window_size}.number_of_bins.{number_of_bins}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.number_of_bins.{number_of_bins}.rbest.snd.liftover.ubcs.{file_type}.tsv', target = 'hg38', query = ['panTro5'], window = 5, outgroup = 'rheMac8',  derived = ['target', 'query'], type= 'target', chromosome = chroms, window_size = [180, 240, 300], region = [ 250 * 1000, 500 * 1000, 1000 * 1000], number_of_bins = [1, 2, 4,5,10], file_type = ['summary', 'details'] ),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}.window_size.{window_size}.number_of_bins.{number_of_bins}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.number_of_bins.{number_of_bins}.rbest.snd.liftover.ubcs.{file_type}.tsv', target = 'hg38', query = ['panTro5'], window = 5, outgroup = 'rheMac8',  derived = ['target', 'query'], type= 'target', chromosome = chroms, window_size = [180], region = [ 250 * 1000, 500 * 1000, 1000 * 1000], number_of_bins = [180], file_type = ['summary', 'details'] ),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}.window_size.{window_size}.number_of_bins.{number_of_bins}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.number_of_bins.{number_of_bins}.rbest.snd.liftover.ubcs.{file_type}.tsv', target = 'hg38', query = ['panTro5'], window = 5, outgroup = 'rheMac8',  derived = ['target', 'query'], type= 'target', chromosome = chroms, window_size = [240], region = [ 250 * 1000, 500 * 1000, 1000 * 1000], number_of_bins = [240], file_type = ['summary', 'details'] ),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}.window_size.{window_size}.number_of_bins.{number_of_bins}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.number_of_bins.{number_of_bins}.rbest.snd.liftover.ubcs.{file_type}.tsv', target = 'hg38', query = ['panTro5'], window = 5, outgroup = 'rheMac8',  derived = ['target', 'query'], type= 'target', chromosome = chroms, window_size = [300], region = [ 250 * 1000, 500 * 1000, 1000 * 1000], number_of_bins = [300], file_type = ['summary', 'details'] ),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}.window_size.{window_size}.number_of_bins.{number_of_bins}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.number_of_bins.{number_of_bins}.rbest.snd.liftover.ubcs.{file_type}.tsv', target = 'hg38', query = apes, window = 5, outgroup = 'rheMac8',  derived = ['target', 'query'], type= 'target', chromosome = chroms, window_size = [300], region = [ 1000 * 1000], number_of_bins = [300], file_type = ['summary', 'details'] ),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}.window_size.{window_size}.number_of_bins.{number_of_bins}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.number_of_bins.{number_of_bins}.rbest.snd.liftover.ubcs.{file_type}.tsv', target = 'hg38', query = 'gorGor5', window = 5, outgroup = 'rheMac8',  derived = ['target', 'query'], type= 'target', chromosome = chroms, window_size = [300], region = [ 1000 * 1000, 100 * 1000, 50 * 1000], number_of_bins = [2], file_type = ['summary', 'details'] )


rule generate_bin_statistics:
    input:
        liftover_fasta_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.rbest.snd.liftover.fasta'
    output:
        stats_tsv_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.{chromosome}.v{window}.derived_in_{derived}.outgroup.{outgroup}.bin_size_stats.tsv'
    script:
        'scripts/python/getBinSizeStats.py'

rule join_bin_statistics:
    input:
        stats_tsv_file = lambda wildcards: expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.{chromosome}.v{window}.derived_in_{derived}.outgroup.{outgroup}.bin_size_stats.tsv', chromosome = chroms,  **wildcards)
    output:
        stats_tsv_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.v{window}.derived_in_{derived}.outgroup.{outgroup}.bin_size_stats.tsv'
    shell:
        'cat {input} > {output}'	 


rule generate_compressed_bins_statistics:
    input:
        liftover_fasta_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.rbest.snd.liftover.fasta'
    output:
        stats_tsv_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.{chromosome}.v{window}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.compressed_bins_stats.tsv'
    script:
        'scripts/python/getCompressedBinsStats.py'


rule join_compressed_bins_statistics:
    input:
        stats_tsv_file = lambda wildcards: expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.{chromosome}.v{window}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.compressed_bins_stats.tsv', chromosome = chroms,  **wildcards)
    output:
        stats_tsv_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/{type}.{target}.{query}.v{window}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.compressed_bins_stats.tsv'
    shell:
        'cat {input} > {output}'

rule generate_ubcs_statistics:
    input:
        liftover_fasta_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.rbest.from.{axt_dir}.snd.liftover.fasta'
    output:
        ubcs_stats_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}.window_size.{window_size}.number_of_bins.{number_of_bins}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.number_of_bins.{number_of_bins}.rbest.from.{axt_dir}.snd.liftover.ubcs.summary.tsv',
        ubcs_stats_details_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/ubcs_fast/derived_in_{derived}/r{region}.window_size.{window_size}.number_of_bins.{number_of_bins}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.window_size.{window_size}.number_of_bins.{number_of_bins}.rbest.from.{axt_dir}.snd.liftover.ubcs.details.tsv'
    script:
        'scripts/python/getUBCSFastStats.py'

	
