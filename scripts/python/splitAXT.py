from scripts.python.utils import get_chroms
from collections import defaultdict
import gzip

def split_AXT(alignment_file, target, query, ffrom, fn_pattern):

    chromosomes =  defaultdict(list)

    with gzip.open(alignment_file, 'rt') as f:
    
        while True:
            header = f.readline()

            if not header:
                break

            if ffrom == 'target':
                p_chrom = header.split()[1]
            else:
                p_chrom = header.split()[4]

            chrom_list = chromosomes[p_chrom]

            chrom_list.append((header, f.readline(), f.readline(), f.readline()))
        

    for chromosome in chromosomes:

        file_name = fn_pattern.format(target = target, query = query, ffrom = ffrom, chromosome = chromosome)

        with open(file_name, 'w') as f:

            for alignment in chromosomes[chromosome]:
                for line in alignment:
                    f.write(line)
