suppressMessages(library('liftOver'))
suppressMessages(library('rtracklayer'))

chain_file <- snakemake@input[['chain_file']]
SNDs_file <- snakemake@input[['SNDs_file']]
liftover_file <- snakemake@output[['liftover_file']]

chain <- import.chain(chain_file)
SNDs <- import(SNDs_file, format = 'BED')
strand(SNDs) <- '+'

after_liftover <- liftOver(SNDs , chain)
after_liftover <- unlist(after_liftover)

export(after_liftover, liftover_file, format = 'BED')
