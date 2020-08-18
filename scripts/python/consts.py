QUALITY_WINDOW_LENGTH = 5
CLUSTERED_SUBSTITUTIONS_WINDOW_LENGTH = 300
CLUSTERED_SUBSTITUTIONS_MIN_SIZE = 5
WTS_MIN_PERCENT = 0.8
STW_MAX_PERCENT = 1 - WTS_MIN_PERCENT

WEAK_NUCLEOTIDES = ['A', 'T']
STRONG_NUCLEOTIDES = ['G', 'C']

#from gorGor5 import gorGor5_chroms

CHROMOSOMES = { 'hg38': [ 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'],
                'panTro5': ['chr1', 'chr2A', 'chr2B', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'],
                'panTro6': ['chr1', 'chr2A', 'chr2B', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'],
                'panPan2': ['chr1', 'chr2A', 'chr2B', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrUn_NW_014024332v1', 'chrUn_NW_014024333v1'],
                'papAnu4': ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chrX'],
                'ponAbe3': ['chr1', 'chr2A', 'chr2B', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chr16_NW_019937288v1_random','chr7_NW_019937272v1_random', 'chr7_NW_019937273v1_random', 'chrX_NW_019937303v1_random', 'chr17_NW_019937293v1_random', 'chr3_NW_019937269v1_random', 'chr16_NW_019937290v1_random', 'chrUn_NW_019937321v1', 'chr16_NW_019937289v1_random'],
                'nomLeu3': [ 'chr1a',  'chr1b', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7b', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22a', 'chr23', 'chr24', 'chr25', 'chrX'],
 #               'gorGor5': gorGor5_chroms
                
                } 



