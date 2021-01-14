from scripts.python.splitAXT import split_AXT

DATA_FOLDER = 'data1'

rule all:
    input:
         expand(DATA_FOLDER + '/{target}/snps/{target}.{chromosome}.snps.bed', target = 'hg38', chromosome = 'chr1'),
         expand(DATA_FOLDER + '/{query}/{query}.{target}.rbest.axt.gz', query = [ 'panTro6', 'panPan3', 'gorGor6', 'ponAbe3' ], target = 'hg38'),
         expand(DATA_FOLDER + '/{query}/{target}.{query}.rbest.axt.gz', query = [ 'nomLeu3' ], target = 'hg38'),
         expand(DATA_FOLDER + '/{query}/{chromosome}.{query}.{target}.from.query.rbest.axt', query = [ 'panTro6', 'panPan3', 'gorGor6', 'ponAbe3' ], target = 'hg38', chromosome = 'chr1'),
         expand(DATA_FOLDER + '/{query}/{chromosome}.{target}.{query}.from.target.rbest.axt', target = 'hg38', query = 'nomLeu3', chromosome = 'chr1')
	 

rule split_snps:
    input:
        snps_file = DATA_FOLDER + '/{target}/{target}.snps.bed'
    output:
        DATA_FOLDER + '/{target}/snps/{target}.{chromosomes}.snps.bed'
    params:
        data_folder = DATA_FOLDER
    shell:
        """   awk '{{ file = sprintf("{params.data_folder}/{wildcards.target}/snps/{wildcards.target}.%s.snps.bed", $1, ".bed"); print >> file }}' {input}  """


rule download_snps:
    output:
        DATA_FOLDER + '/{target}/{target}.snps.bb' 
    params:
        data_folder = DATA_FOLDER
    shell:
        ' wget -q -P  {params.data_folder}/{wildcards.target} http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153Common.bb; '
        ' mv {params.data_folder}/{wildcards.target}/dbSnp153Common.bb {params.data_folder}/{wildcards.target}/{wildcards.target}.snps.bb'

rule download_bigBedToBed:
    output:
        DATA_FOLDER + '/bigBedToBed'
    params:
        data_folder = DATA_FOLDER
    shell:
        ' wget -q -P {params.data_folder} http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed ; '
        ' chmod +x {params.data_folder}/bigBedToBed'


rule convert_BigBedtoBed:
    input:
        bb = DATA_FOLDER + '/{target}/{target}.snps.bb',
        bigBedToBed = DATA_FOLDER + '/bigBedToBed'
    output:
        DATA_FOLDER + '/{target}/{target}.snps.bed'
    shell:
        './{input.bigBedToBed} {input.bb} {output}'

rule download_chrom_sizes:
    output:
        DATA_FOLDER + '/{query}/{query}.chrom.sizes'
    params:
        data_folder = DATA_FOLDER
    shell:
        'wget -q -P {params.data_folder}/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.query}/bigZips/{wildcards.query}.chrom.sizes '

rule download_ref:
    output:
        DATA_FOLDER + '/{query}/{query}.fa.gz'
    params:
        data_folder = DATA_FOLDER
    shell:
        'wget -q -P {params.data_folder}/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.query}/bigZips/{wildcards.query}.fa.gz'

uppercase_first_letter = lambda name:  lambda wildcards:  wildcards[name][0].upper() + wildcards[name][1:]

rule download_query_target_axt:
    output:
        DATA_FOLDER + '/{query}/{query}.{target}.rbest.axt.gz'
    params:
        query_upper_first_letter = uppercase_first_letter('target'),
        data_folder = DATA_FOLDER 
    shell:
        'wget -q -P {params.data_folder}/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.query}/vs{params.query_upper_first_letter}/reciprocalBest/axtRBestNet/{wildcards.query}.{wildcards.target}.rbest.axt.gz '

rule download_target_query_axt:
    output:
        DATA_FOLDER + '/{query}/{target}.{query}.rbest.axt.gz'
    params:
        query_upper_first_letter = uppercase_first_letter('query'),
        data_folder = DATA_FOLDER
    shell:
        'wget -q -P {params.data_folder}/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.target}/vs{params.query_upper_first_letter}/reciprocalBest/axtRBestNet/{wildcards.target}.{wildcards.query}.rbest.axt.gz '


rule split_alignments_from_target:
    input:
        alignment_file = DATA_FOLDER + '/{query}/{target}.{query}.rbest.axt.gz'
    output:
        alignment_file = DATA_FOLDER + '/{query}/{chromosome}.{target}.{query}.from.target.rbest.axt'
    run:
        split_AXT(input.alignment_file, wildcards.target, wildcards.query, 'target', DATA_FOLDER + '/{query}/{chromosome}.{target}.{query}.from.target.rbest.axt' )


rule split_alignments_from_query:
    input:
        alignment_file = DATA_FOLDER + '/{query}/{query}.{target}.rbest.axt.gz'
    output:
        alignment_file = DATA_FOLDER + '/{query}/{chromosome}.{query}.{target}.from.query.rbest.axt'
    run:
        split_AXT(input.alignment_file, wildcards.target, wildcards.query, 'query', DATA_FOLDER + '/{query}/{chromosome}.{query}.{target}.from.query.rbest.axt')
