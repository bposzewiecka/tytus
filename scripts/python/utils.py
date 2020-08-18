def get_chrom_sizes(chrom_sizes_file):
    
    with open(chrom_sizes_file) as f_in:

        chrom_sizes = {}

        for line in f_in:
            line = line.rstrip().split()
            chrom_name = line[0]
            chrom_size = int(line[1])

            chrom_sizes[chrom_name] = chrom_size

        return chrom_sizes

def get_chroms(chrom_sizes_file):

    with open(chrom_sizes_file) as f_in:

        chroms = []

        for line in f_in:
            line = line.rstrip().split()
            chroms.append(line[0])

        return chroms

