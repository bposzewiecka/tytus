configfile: "config.yaml"

DATA_FOLDER = config['data_folder']

chromosomes = {}
chromosomes['hg38'] = [ 'chr' + str(chrom) for chrom in list(range(1,23)) + ['X', 'Y']]

rule main:
    input:
        [ expand( DATA_FOLDER +  '/{query}/outgroup.{outgroup}/{chromosome}.{target}.{query}.outgroup.{outgroup}.from.{ffrom}.rbest.fasta', chromosome = chromosomes[sample['target']], **sample) for sample in config['samples']],
        [ expand(DATA_FOLDER + '/{query}/outgroup.{outgroup}/derived_in_{derived}/region.{region}.window.{window_size}.bins.{number_of_bins}.with_snps.{with_snps}/{chromosome}.{target}.{query}.outgroup.{outgroup}.derived_in_{derived}.region.{region}.window.{window_size}.bins.{number_of_bins}.with_snps.{with_snps}.from.{ffrom}.ubcs.summary.tsv', chromosome = chromosomes[sample['target']], **sample, **ubcs_params) for sample in config['samples'] for ubcs_params in config['ubcs_params'] if sample['query'] == 'panTro6'],


rule generate_ubcs_statistics:
    input:
        liftover_fasta_file = DATA_FOLDER + '/{query}/outgroup.{outgroup}/{chromosome}.{target}.{query}.outgroup.{outgroup}.from.{ffrom}.rbest.fasta'
    output:
        ubcs_stats_file = DATA_FOLDER + '/{query}/outgroup.{outgroup}/derived_in_{derived}/region.{region}.window.{window_size}.bins.{number_of_bins}.with_snps.{with_snps}/{chromosome}.{target}.{query}.outgroup.{outgroup}.derived_in_{derived}.region.{region}.window.{window_size}.bins.{number_of_bins}.with_snps.{with_snps}.from.{ffrom}.ubcs.summary.tsv',
        ubcs_stats_details_file = DATA_FOLDER + '/{query}/outgroup.{outgroup}/derived_in_{derived}/region.{region}.window.{window_size}.bins.{number_of_bins}.with_snps.{with_snps}/{chromosome}.{target}.{query}.outgroup.{outgroup}.derived_in_{derived}.region.{region}.window.{window_size}.bins.{number_of_bins}.with_snps.{with_snps}.from.{ffrom}.ubcs.details.tsv',
	ubcs_stats_log_file = DATA_FOLDER + '/{query}/outgroup.{outgroup}/derived_in_{derived}/region.{region}.window.{window_size}.bins.{number_of_bins}.with_snps.{with_snps}/{chromosome}.{target}.{query}.outgroup.{outgroup}.derived_in_{derived}.region.{region}.window.{window_size}.bins.{number_of_bins}.with_snps.{with_snps}.from.{ffrom}.ubcs.log.tsv'
    script:
        'scripts/python/getUBCSFastStats.py'

rule liftover:
    input:
        SNDs_file = DATA_FOLDER + '/{query}/{chromosome}.{target}.{query}.from.{ffrom}.rbest.bed',
        chain_file = DATA_FOLDER + '/{outgroup}/{target}.{outgroup}.rbest.chain'
    output:
        liftover_file = DATA_FOLDER + '/{query}/outgroup.{outgroup}/{chromosome}.{target}.{query}.outgroup.{outgroup}.from.{ffrom}.rbest.bed'
    script:
        'scripts/R/liftover.R'

rule get_fasta_from_liftover:
    input:
        liftover_file =  DATA_FOLDER + '/{query}/outgroup.{outgroup}/{chromosome}.{target}.{query}.outgroup.{outgroup}.from.{ffrom}.rbest.bed',
        ref_file = DATA_FOLDER + '/{outgroup}/{outgroup}.fa',
	ref_file_fai =  DATA_FOLDER + '/{outgroup}/{outgroup}.fa.fai'
    output:
        DATA_FOLDER +  '/{query}/outgroup.{outgroup}/{chromosome}.{target}.{query}.outgroup.{outgroup}.from.{ffrom}.rbest.fasta'
    shell:
        'bedtools getfasta -s -name -fi {input.ref_file} -bed {input.liftover_file} > {output}'

rule bed_snd_from_axt_target:
    input:
        alignment_file = DATA_FOLDER + '/{query}/{target}.{query}/{chromosome}.{target}.{query}.from.target.rbest.axt',
        query_chrom_sizes_file = DATA_FOLDER + '/{query}/{query}.chrom.sizes',
	snp_file = DATA_FOLDER + '/{target}/snps/{chromosome}.{target}.snps.bed'
    output:
        snd_file = DATA_FOLDER + '/{query}/{chromosome}.{target}.{query}.from.target.rbest.bed'
    params:
        ffrom = 'target'	    
    script:
        'scripts/python/getSNDsFromAXT.py'

rule bed_snd_from_axt_query:
    input:
        alignment_file = DATA_FOLDER + '/{query}/{query}.{target}/{chromosome}.{query}.{target}.from.query.rbest.axt',
        query_chrom_sizes_file = DATA_FOLDER + '/{target}/{target}.chrom.sizes',
        snp_file = DATA_FOLDER + '/{target}/snps/{chromosome}.{target}.snps.bed'
    output:
        snd_file = DATA_FOLDER + '/{query}/{chromosome}.{target}.{query}.from.query.rbest.bed'
    params:
        ffrom = 'query'
    script:
        'scripts/python/getSNDsFromAXT.py'

rule create_fai_index:
    input:
        DATA_FOLDER + '/{outgroup}/{outgroup}.fa'
    output:
        DATA_FOLDER + '/{outgroup}/{outgroup}.fa.fai' 	    
    shell:
        'samtools faidx {input}'	    
