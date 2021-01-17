from consts import WEAK_NUCLEOTIDES 
from consts import STRONG_NUCLEOTIDES

window_length = 5

def get_coords(coords_str):
    
    coords_str = coords_str.split(':')

    chrom = coords_str[0]

    coords_str = coords_str[1].split('-')

    start = int(coords_str[0])
    end = int(coords_str[1])

    return chrom, start, end

class SND:

    def target_base(self):
        return self.target_seq[window_length]

    def query_base(self):
        return self.query_seq[window_length]

    def outgroup_base(self):
        return self.outgroup_seq[window_length]

    def __init__(self, name, outgroup_seq):

        attrs = name.split('::')[0].split(';')
        attrs = {  attr[:attr.find(':')]:  attr[attr.find(':') + 1:]  for attr in  attrs }

        self.target = attrs['TARGET']
        self.query = attrs['QUERY']
        self.target_seq = attrs['CHANGE'][:2 * window_length + 1] 
        self.query_seq =  attrs['CHANGE'][2 * window_length + 2:]

        self.target_chrom, self.target_start, self.target_end = get_coords(attrs['TARGET_COORDS'])

        self.query_chrom, self.query_start, self.query_end = get_coords(attrs['QUERY_COORDS'])
        
        self.snp = attrs['SNP'] == 'True'

        self.outgroup_seq = outgroup_seq.upper()

    def target_coord(self):
        return self.target_start + window_length

    def query_coord(self):
        return self.query_start + window_length
    
    def derrived_in_query(self):
        return self.target_base() == self.outgroup_base()

    def derrived_in_target(self):
        return self.query_base() == self.outgroup_base()

    def biased_in_query(self):
        return self.derrived_in_query() and self.outgroup_base() in WEAK_NUCLEOTIDES and self.query_base() in STRONG_NUCLEOTIDES

    def biased_in_target(self):
        return self.derrived_in_target() and self.outgroup_base() in WEAK_NUCLEOTIDES and self.target_base() in STRONG_NUCLEOTIDES

    def biased(self):
        if self.derrived_in_query():
            return  self.outgroup_base() in WEAK_NUCLEOTIDES and self.query_base() in STRONG_NUCLEOTIDES
        elif self.derrived_in_target():
            return self.outgroup_base() in WEAK_NUCLEOTIDES and self.target_base() in STRONG_NUCLEOTIDES 

        return False 

    def is_snp(self):
        return self.snp


def get_SNDs(file_name):

    SNDs = []

    with open(file_name) as f_in:

        for i, line in enumerate(f_in):
            line = line.rstrip()

            if line.startswith('>'):
                name = line[1:-1]
            else:
                outgroup_seq = line

                if len(outgroup_seq) != 2 * window_length + 1:  continue

                snd = SND(name, outgroup_seq)
                SNDs.append(snd)

        return SNDs

def get_SNDs_derived_in_target(file_name):

    SNDs = get_SNDs(file_name)


    return [ SND  for SND in SNDs if SND.derrived_in_target() ]

def get_SNDs_derived_in_query(file_name):

    SNDs = get_SNDs(file_name)

    return [ SND  for SND in SNDs if SND.derrived_in_query() ]

