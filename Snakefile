#configfile: 'config.yaml'

DATA_FOLDER = 'data1'

hg_chroms = [ 'chr' + str(chrom) for chrom in list(range(1,23)) + ['X', 'Y']]

rule main:
    input:
        expand(DATA_FOLDER + '/{query}/{chromosome}.{target}.{query}.from.target.rbest.bed', chromosome = 'chr1', query = 'nomLeu3', target = 'hg38'),
        expand(DATA_FOLDER + '/{query}/{chromosome}.{target}.{query}.from.query.rbest.bed', chromosome = 'chr1', query = 'panTro6', target = 'hg38'),

rule get_fasta_from_liftover:
    input:
        bed_file = DATA_FOLDER + '/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.{file_type}.from.{axt_dir}.snd.liftover.bed',
        ref_file = DATA_FOLDER + '/{outgroup}/{outgroup}.fa'
    output:
        'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.{file_type}.from.{axt_dir}.snd.liftover.fasta'
    shell:
        'bedtools getfasta -s -name -fi {input.ref_file} -bed {input.bed_file} > {output}'

rule liftover:
    input:
        SNDs_file = 'data/{query}/{target}.{query}/{type}/v{window}/{type}.{target}.{query}.{chromosome}.v{window}.from.{axt_dir}.rbest.snd.bed',
        chain_file = lambda wildcards : 'data/{outgroup}/{target}.{outgroup}.rbest.chain'
    output:
        liftover_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.{file_type}.from.{axt_dir}.snd.liftover.bed'
    script:
        'scripts/R/liftover.R'

rule bed_snd_from_anx_target:
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

rule bed_snd_from_anx_query:
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

