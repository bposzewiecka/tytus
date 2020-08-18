import os
import consts

from utils import get_chroms

alignment_file = snakemake.input['alignment_file']
chromosome_file = snakemake.input['chromosome_file']

log_file = snakemake.output[0]

target = snakemake.wildcards['target']
query = snakemake.wildcards['query']
ttype = snakemake.wildcards['type']

directory = os.path.join('data', query, target + '.' + query, ttype)

if not os.path.exists(directory):
    os.makedirs(directory)

chromosomes =  { chromosome: []  for chromosome in get_chroms(chromosome_file)}

with open(alignment_file) as f:
    
    while True:
        header = f.readline()

        if not header:
            break

        if ttype == 'target':
            p_chrom = header.split()[1]
        else:
            p_chrom = header.split()[4]

        chrom_list = chromosomes[p_chrom]

        chrom_list.append((header, f.readline(), f.readline(), f.readline()))
        

for chromosome in chromosomes:
    desc_name = 'data/{query}/{target}.{query}/{ttype}/{ttype}.{target}.{query}.{chromosome}.rbest.axt'.format(directory = directory, chromosome = chromosome, target = target, query = query, ttype = ttype)
    with open(desc_name, 'w') as f:

        for alignment in chromosomes[chromosome]:
            for line in alignment:
                f.write(line)

with open(log_file, 'w') as f:
    for chromosome in chromosomes:
        f.write(chromosome + ' ' + str(len(chromosomes[chromosome])) + '\n')


