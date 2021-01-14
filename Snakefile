# wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153Common.bb
# wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
# ./bigBedToBed dbSnp153Common.bb dbSnp153Common.bed

#configfile: 'config.yaml'

#from scripts.python.consts import CHROMOSOMES
from scripts.python.utils import get_chrom_sizes

MIN_CHROM_SIZE = 1000 * 1000

apes = ['panTro6',  'papAnu4', 'ponAbe3', 'panPan2', 'panTro5', 'nomLeu3', 'gorGor5', 'gorGor4'] # 'chlSab2', 'macFas5', 'rheMac8', 'nasLar1', 'rhiRox1', 'calJac3', 'saiBol1', 'tarSyr2', 'micMur2', 'otoGar3'] 

hg_chroms = [ 'chr' + str(chrom) for chrom in list(range(1,23)) + ['X', 'Y']]

arch_homini_chroms = list(range(1,23)) + ['X', 'Y']
arch_homini = ['AltaiNea', 'Denisova']

def get_list(o):
     if isinstance(o, list): 
         return o
     else:
         return [o]
    
def expand_chroms(pattern, type, query, target, **kwargs):
    queries = get_list(query)
    targets = get_list(target)
    ttypes = get_list(type)

    tmp = []
   
    for ttype in ttypes:
        for target in targets:
            for query in queries:

                ref = target if ttype == 'target' else query
                chrom_sizes_file = 'data/{ref}/{ref}.chrom.sizes'.format(ref = ref)
       	        chromosomes = [ chrom for chrom, chrom_size in get_chrom_sizes(chrom_sizes_file).items() if chrom_size > MIN_CHROM_SIZE ]
           
                tmp += expand(pattern, chromosome = chromosomes, type = ttype, query = query , target = target , **kwargs)
    return tmp

rule all:
    input:
        expand('data/{target}/snps/split.log', target = 'hg38'),	    
        expand('data/{query}/{target}.{query}/{type}/query.alignment.split.log', target = 'hg38', query = ['panPan3', 'panTro6', 'gorGor6', 'ponAbe3'], type ='target'),
        expand('data/{query}/{target}.{query}/{type}/v{window}/{type}.{target}.{query}.{chromosome}.v{window}.from.{axt_dir}.rbest.snd.bed',  target = 'hg38', query = ['panPan3', 'panTro6', 'gorGor6', 'ponAbe3'], type ='target', chromosome = hg_chroms, window = 5, axt_dir = 'query' ),
	expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.rbest.from.{axt_dir}.snd.liftover.bed', target = 'hg38', query = ['panPan3', 'panTro6', 'gorGor6', 'ponAbe3'], type ='target', chromosome = hg_chroms, window = 5, axt_dir = 'query', outgroup ='rheMac10' ),
	expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.rbest.from.{axt_dir}.snd.liftover.fasta', target = 'hg38', query = ['panPan3', 'panTro6', 'gorGor6', 'ponAbe3'], type ='target', chromosome = hg_chroms, window = 5, axt_dir = 'query', outgroup ='rheMac10' ),

        #expand('data/{query}/{target}.{query}/{type}/split.log', query = ['panTro6'], target = 'hg38', type = 'target'),
	#expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.rbest.snd.liftover.fasta', target = 'hg38', query = 'panTro6', outgroup = 'rheMac10',type ='target', window=5, chromosome = hg_chroms) 
        #expand('data/{query}/{target}.{query}/{type}/split.log', query = ['panTro4', 'hg38'], target = 'panPan2', type = ['target', 'query']),
	#expand('data/{query}/{target}.{query}/all/v{window}/outgroup.{outgroup}/all.{target}.{query}.v{window}.r{region}.outgroup.{outgroup}.rbest.snd.liftover.ubcs.tsv', target = 'hg38', query = ['panTro6'], window = 5, outgroup= 'rheMac10', region = 1000 * 1000 ),
	#expand('data/{target}/snps/split.log', target = 'hg38')
        #expand('data/{query}/{target}.{query}/{type}/split.log', query = 'panTro4', target = 'panPan2', type = ['target', 'query']),
        #expand('data/{query}/{target}.{query}/all/v{window}/outgroup.{outgroup}/all.{target}.{query}.v{window}.r{region}.outgroup.{outgroup}.rbest.snd.liftover.ubcs.tsv', target = 'panPan2', query = 'panTro4',window = 5, outgroup= 'gorGor5', region = 1000 * 1000 ),
        #expand('data/{query}/{target}.{query}/{type}/split.log', query = apes, target = 'panPan2', type = ['target', 'query']),
        #expand('data/{target}/all/v{window}/outgroup.{outgroup}/all.{target}.all.v{window}.r{region}.outgroup.{outgroup}.rbest.snd.liftover.ubcs.tsv', target = 'hg38', window = 5, outgroup= 'rheMac8', region = 1000 * 1000  ),
        #expand('data/{species}/{species}.hg19_1000g.{chromosome}.snd.ac2.vcf', chromosome =  arch_homini_chroms, species = arch_homini),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.chr{chromosome}.v{window}.outgroup.{outgroup}.vcf.ac2.snd.liftover.fasta', query = arch_homini, target = 'hg19', type = 'target', window = 5, chromosome = arch_homini_chroms, outgroup = 'panTro6'),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/derived_in_{derived}/r{region}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.vcf.ac2.snd.liftover.ubcs.tsv', target = 'hg19', query = arch_homini, window = 5, outgroup= 'panTro6', region = 8000 * 1000, type = ['target'], derived = ['target', 'query'], chromosome = hg_chroms),
        #expand('data/{query}/{target}.{query}/{type}/split.log', query = apes, target = 'hg38', type = ['target', 'query']),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/derived_in_{derived}/r{region}/{type}.{target}.{query}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.rbest.snd.liftover.ubcs.tsv', target = 'hg38', query = apes, window = 5, outgroup= 'rheMac8', region = 1000 * 1000, type = ['target', 'query'], derived = ['target', 'query'])
        #expand('data/{species}/{species}.hg19_1000g.{chromosome}.snd.ac2.vcf', chromosome =  arch_homini_chroms, species = arch_homini),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.vcf.snd.liftover.fasta', query = 'Denisova', target = 'hg19', type = 'target', window = 5, chromosome = CHROMOSOMES['hg38'], outgroup='panTro6'),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/derrived_in_{derived}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.vcf.snd.liftover.ubcs.tsv', query = arch_homini, target = 'hg19', type = 'target', window = 5, chromosome = 'chr2', outgroup='panTro6', derived = 'query', region = 1000 * 1000),
        #expand('data/{query}/{query}.hg38.rbest.{ftype}.gz', query = apes, ftype = ['chain', 'net']),
        #expand('data/{species}/{species}.hg19_1000g.{chromosome}.snd.vcf', chromosome =  arch_homini_chroms, species = arch_homini),
        #expand('data/{query}/{target}.{query}/{type}/{type}.{target}.{query}.chr1.rbest.axt', query = apes, target = 'hg38', type = ['target']),
        #expand_chroms('data/{query}/{target}.{query}/{type}/v{window}/{type}.{target}.{query}.{chromosome}.v{window}.rbest.snd.bed', query = apes, target = 'hg38', window = 5, ttype = 'target' ),
        #expand_chroms('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.rbest.snd.liftover.bed', query = apes, target = 'hg38', window = 5, ttype = 'target', outgroup = 'rheMac8'),
        #expand_chroms('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.rbest.snd.liftover.fasta', query = apes, target = 'hg38', window = 5, ttype = 'target', outgroup = 'rheMac8'),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/{type}.{target}.{query}.chr{chromosome}.v{window}.vcf.snd.bed', query = 'Denisova', target = 'hg19', type = 'target', window = 5, chromosome = arch_homini_chroms),
        #expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.vcf.snd.liftover.bed', query = 'Denisova', target = 'hg19', type = 'target', window = 5, chromosome = CHROMOSOMES['hg38'], outgroup='panTro6')


rule get_UBCS_stats_all_species:
    input:
        lambda wildcards: expand('data/{query}/{target}.{query}/all/v{window}/outgroup.{outgroup}/all.{target}.{query}.v{window}.r{region}.outgroup.{outgroup}.{file_type}.snd.liftover.ubcs.tsv',  **wildcards, query=apes)
    output:
        'data/{target}/all/v{window}/outgroup.{outgroup}/all.{target}.all.v{window}.r{region}.outgroup.{outgroup}.{file_type}.snd.liftover.ubcs.tsv'
    shell:
        'cat {input} > {output}'


rule get_UBCS_stats_all_types:
    input:
        lambda wildcards: expand('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/derived_in_{derived}/r{region}/{type}.{target}.{query}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.{file_type}.snd.liftover.ubcs.tsv', **wildcards, type = ['target', 'query'], derived = ['target', 'query'] )
    output:
        'data/{query}/{target}.{query}/all/v{window}/outgroup.{outgroup}/all.{target}.{query}.v{window}.r{region}.outgroup.{outgroup}.{file_type}.snd.liftover.ubcs.tsv'
    shell:
        'cat {input} > {output}'

rule get_UBCS_stats_all_chroms:
    input:
        lambda wildcards: expand_chroms('data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/derived_in_{derived}/r{region}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.{file_type}.snd.liftover.ubcs.tsv', **wildcards)
    output:
        'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/derived_in_{derived}/r{region}/{type}.{target}.{query}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.{file_type}.snd.liftover.ubcs.tsv'
    shell:
        'cat {input} > {output}'

rule get_UBCS_stats:
    input:
        liftover_fasta_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.{chromosome}.v{window}.outgroup.{outgroup}.{file_type}.snd.liftover.fasta',
        chrom_sizes_file = 'data/{target}/{target}.chrom.sizes'
    output:
        ubcs_stats_file =   'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/derived_in_{derived}/r{region}/{type}.{target}.{query}.{chromosome}.v{window}.r{region}.derived_in_{derived}.outgroup.{outgroup}.{file_type}.snd.liftover.ubcs.tsv'
    script:
        'scripts/python/getUBCSStats.py'

rule create_snd_chrom_file_vcf:
    input:
        vcf_file = 'data/{query}/{query}.{target}_1000g.{chromosome}.mod.vcf.gz'
    output:
        snd_file = 'data/{query}/{target}.{query}/{type}/v{window}/{type}.{target}.{query}.chr{chromosome}.v{window}.vcf.snd.bed'
    script:
        'scripts/python/getSNDsFromVCF.py'

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

rule split_alignments1:
    input:
        alignment_file = 'data/{query}/{target}.{query}.rbest.axt',
        chromosome_file = lambda wildcards : 'data/{ref}/{ref}.chrom.sizes'.format(ref = wildcards.target if wildcards.type == 'target' else wildcards.query)
    output:
        'data/{query}/{target}.{query}/{type}/target.alignment.split.log'
    params:
        axt_dir = 'target'	    
    script: 
        'scripts/python/splitAXT.py'

rule split_alignments2:
    input:
        alignment_file = 'data/{query}/{query}.{target}.rbest.axt',
        chromosome_file = lambda wildcards : 'data/{ref}/{ref}.chrom.sizes'.format(ref = wildcards.target if wildcards.type == 'target' else wildcards.query)
    output:
        'data/{query}/{target}.{query}/{type}/query.alignment.split.log'
    params:
        axt_dir = 'query'
    script:
        'scripts/python/splitAXT.py'

rule create_snd_chrom_file_rbest_inv:
    input:
        alignment_file = 'data/{query}/{target}.{query}/{type}/{type}.{target}.{query}.{chromosome}.from.{axt_dir}.rbest.axt',
	query_chrom_sizes_file = 'data/{target}/{target}.chrom.sizes',
        snp_file = 'data/{target}/snps/{chromosome}.{target}.snps.bed'
    output:
        snd_file = 'data/{query}/{target}.{query}/{type}/v{window}/{type}.{target}.{query}.{chromosome}.v{window}.from.{axt_dir}.rbest.snd.bed'
    script:
        'scripts/python/getSNDsFromAXTinv.py'


rule split_snps:
    input:
        snps_file = 'data/{target}/{target}.snps.bed'
    output:
        'data/{target}/snps/split.log'
    shell:
        """
	   awk '{{ file = sprintf("data/{wildcards.target}/snps/%s.%s.snps.bed", $1, "{wildcards.target}", ".bed"); print >> file }} END {{ echo "END" }} ' {input} > {output} 
	"""

rule filter_allele_count2:
    input:
        vcf_ac2_file = 'data/{query}/{query}.hg19_1000g.{chromosome}.snd.ac2.vcf',
        liftover_fasta_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.chr{chromosome}.v{window}.outgroup.{outgroup}.vcf.snd.liftover.fasta'
    output:
        liftover_fasta_file = 'data/{query}/{target}.{query}/{type}/v{window}/outgroup.{outgroup}/{type}.{target}.{query}.chr{chromosome}.v{window}.outgroup.{outgroup}.vcf.ac2.snd.liftover.fasta'
    shell:
        """ awk '{{ if (substr($0,1 ,1) != "#") print "chr" $1  ":" $2 - 6  "-" $2 + 5  }}' {input.vcf_ac2_file}  | """ 
        """ grep -A 1  -f - {input.liftover_fasta_file}  |  """
        """ awk '{{ if ($0 != "--") print $0 }}' > {output} """

rule create_snd_vcf:
    input:
        'data/{species}/{species}.hg19_1000g.{chromosome}.mod.vcf.gz'
    output:
        'data/{species}/{species}.hg19_1000g.{chromosome}.snd.vcf'
    shell:
        """ gunzip -cd {input}  | """
        """ awk '{{ if (substr($0,1,1) == "#"  || ($5 != "." && length($4) == 1 && length($5) == 1)) print $0 }}' > {output} """

rule create_snd_ac2_vcf:
    input:
        'data/{species}/{species}.hg19_1000g.{chromosome}.snd.vcf'
    output:
        'data/{species}/{species}.hg19_1000g.{chromosome}.snd.ac2.vcf'
    shell:
        """ awk '{{ if (substr($0,1,1) == "#"  || ($5 != "." && length($4) == 1 && length($5) == 1 && substr($8 , 1, 4) == "AC=2")) print $0 }}' > {output} """

rule download_ref:
    output:
        'data/{query}/{query}.fa.gz'
    log:
        out = 'data/{query}/out/{query}.fa.out',
        err = 'data/{query}/err/{query}.fa.err'
    shell:
        'wget -P data/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.query}/bigZips/{wildcards.query}.fa.gz '
        ' > {log.out} 2> {log.err}'

uppercase_first_letter = lambda wildcards: wildcards['query'][0].upper() + wildcards['query'][1:]

rule download_chain_file:
    output:
        'data/{query}/{target}.{query}.rbest.net.gz'
    params:
        query_upper_first_letter = uppercase_first_letter
    log:
        out = 'data/{query}/out/{target}.{query}.rbest.chain.out',
        err = 'data/{query}/err/{target}.{query}.rbest.chain.err'
    shell:
        'wget -P data/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.target}/vs{params.query_upper_first_letter}/reciprocalBest/{wildcards.target}.{wildcards.query}.rbest.net.gz'
        ' > {log.out} 2> {log.err}'

rule download_chrom_sizes:
    output:
        'data/{query}/{query}.chrom.sizes'
    log:
        out = 'data/{query}/out/{query}.chrom.sizes.out',
        err = 'data/{query}/err/{query}.chrom.sizes.err'
    shell:
        'wget -P data/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.query}/bigZips/{wildcards.query}.chrom.sizes '
        ' > {log.out} 2> {log.err}'

ruleorder: download_axt > download_hg_ref
ruleorder: download_axt > download_ref_hg

rule download_hg_ref:
    output:
        'data/{query}/hg38.{query}.rbest.{ftype}.gz'
    params:
        query_upper_first_letter = uppercase_first_letter
    log:
        out = 'data/{query}/out/hg38.{query}.rbest.{ftype}.out',
        err = 'data/{query}/err/hg38.{query}.rbest.{ftype}.err'
    shell:
        'wget -P data/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vs{params.query_upper_first_letter}/reciprocalBest/hg38.{wildcards.query}.rbest.{wildcards.ftype}.gz '
        ' > {log.out} 2> {log.err}'


rule download_hg19_ref:
    output:
        'data/{query}/hg19.{query}.rbest.{ftype}.gz'
    params:
        query_upper_first_letter = uppercase_first_letter
    log:
        out = 'data/{query}/out/hg19.{query}.rbest.{ftype}.out',
        err = 'data/{query}/err/hg19.{query}.rbest.{ftype}.err'
    shell:
        'wget -P data/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/hg19/vs{params.query_upper_first_letter}/reciprocalBest/hg19.{wildcards.query}.rbest.{wildcards.ftype}.gz '
        ' > {log.out} 2> {log.err}'
     
rule download_ref_hg:
    output:
        'data/{query}/{query}.hg38.rbest.{ftype}.gz'
    params:
        query_upper_first_letter = uppercase_first_letter
    log:
        out = 'data/{query}/out/{query}.hg38.rbest.{ftype}.out',
        err = 'data/{query}/err/{query}.hg38.rbest.{ftype}.err'
    shell:
        'wget -P data/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vs{params.query_upper_first_letter}/reciprocalBest/{wildcards.query}.hg38.rbest.{wildcards.ftype}.gz '
        ' > {log.out} 2> {log.err}'

rule download_axt:
    output:
        'data/{query}/{target}.{query}.rbest.axt.gz'
    params:
        query_upper_first_letter = uppercase_first_letter
    log:
        out = 'data/{query}/out/{target}.{query}.rbest.axt.out',
        err = 'data/{query}/err/{target}.{query}.rbest.axt.err'
    shell:
        'wget -P data/{wildcards.query} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.target}/vs{params.query_upper_first_letter}/reciprocalBest/axtRBestNet/{wildcards.target}.{wildcards.query}.rbest.axt.gz '
        ' > {log.out} 2> {log.err}'

rule download_altai_neandertal_vcf:
    output:
        'data/AltaiNea/AltaiNea.hg19_1000g.{chromosome}.mod.vcf.gz'
    log:
        out = 'data/AltaiNea/out/AltaiNea.hg19_1000g.{chromosome}.mod.vcf.out',
        err = 'data/AltaiNea/err/AltaiNea.hg19_1000g.{chromosome}.mod.vcf.err'
    shell:
        'wget -P data/AltaiNea http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.{wildcards.chromosome}.mod.vcf.gz'
        ' > {log.out} 2> {log.err}'

rule download_denisova_vcf:
    output:
        'data/Denisova/Denisova.hg19_1000g.{chromosome}.mod.vcf.gz'
    log:
        out = 'data/Denisova/out/Denisova.hg19_1000g.{chromosome}.mod.vcf.out',
        err = 'data/Denisova/err/Denisova.hg19_1000g.{chromosome}.mod.vcf.err'
    shell:
        'wget -P data/Denisova http://cdna.eva.mpg.de/denisova/VCF/hg19_1000g/T_hg19_1000g.{wildcards.chromosome}.mod.vcf.gz'
        ' > {log.out} 2> {log.err}; '
        'mv data/Denisova/T_hg19_1000g.{wildcards.chromosome}.mod.vcf.gz data/Denisova/Denisova.hg19_1000g.{wildcards.chromosome}.mod.vcf.gz '

ruleorder: download_hg_ref > gunzip
ruleorder: download_hg19_ref > gunzip
ruleorder: download_ref_hg > gunzip
ruleorder: liftover > gunzip
ruleorder: get_UBCS_stats > gunzip
#ruleorder: split_alignments > gunzip

rule gunzip:
    input:
        '{name}.gz'
    output:
        '{name}1'
    shell:
        'gunzip -k {input}'
