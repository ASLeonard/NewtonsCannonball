

include: 'snakepit/coverage.smk'
include: 'snakepit/SRF.smk'


rule all:
    input:
        'coverage/regions.csv',
        'coverage/genome.csv',
        expand('satellites/{sample}.abundance.txt',sample=config['samples']),
        expand('satellites/{sample}.srf.fa.full_tbl',sample=config['samples'])
