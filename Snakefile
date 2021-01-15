#configfile: 'config.yaml'

apes = ['panTro6',  'papAnu4', 'ponAbe3', 'panPan2', 'panTro5', 'nomLeu3', 'gorGor5', 'gorGor4']

hg_chroms = [ 'chr' + str(chrom) for chrom in list(range(1,23)) + ['X', 'Y']]

rule get_fasta_from_liftover:
    input:
        bed_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.{file_type}.from.{axt_dir}.snd.liftover.bed',
        ref_file = 'data/{outgroup}/{outgroup}.fa'
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

rule create_snd_chrom_file_rbest:
    input:
        alignment_file = 'data/{query}/{target}.{query}/{type}/{type}.{target}.{query}.{chromosome}.rbest.axt',
        query_chrom_sizes_file = 'data/{query}/{query}.chrom.sizes',
	snp_file = 'data/{target}/snps/{chromosome}.{target}.snps.bed'
    output:
        snd_file = 'data/{query}/{target}.{query}/{type}/v{window}/{type}.{target}.{query}.{chromosome}.v{window}.rbest.snd.bed'
    script:
        'scripts/python/getSNDsFromAXT.py'

