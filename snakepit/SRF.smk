from pathlib import PurePath

rule all:
    input:
        expand('satellites/{sample}.abundance.txt',sample=config['samples'])

rule samtools_fastq:
    input:
        '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long_reads/{sample}/alignment/{sample}.mm2.cram'
    output:
        'reads/{sample}.fa.gz'
    threads: 4
    resources:
        mem_mb = 2500
    shell:
        '''
        samtools cat -@ {threads} {input} |\
        samtools fasta -@ {threads} --reference {config[reference]} - | pigz -p {threads} -c > {output}
        '''

rule estimate_coverage:
    input:
        rules.samtools_fastq.output
    output:
        'kmers/{sample}.coverage'
    #localrules: True
    shell:
        '''
        zgrep -v ">" {input} | awk '{{ l+=length($1) }} END {{ print l/2759153975 }}' > {output}
        '''

rule KMC_count:
    input:
        reads = rules.samtools_fastq.output,
        coverage = rules.estimate_coverage.output
    output:
        multiext('kmers/{sample}.kmc','.kmc_pre','.kmc_suf')
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix(''),
        threshold = lambda wildcards, input: int(float(open(input.coverage[0]).read())*10)
    threads: 8
    resources:
        mem_mb = 8000
    shell:
        '''
        kmc -k151 -t{threads} -ci{params.threshold} -cs100000 -m20 -fa {input.reads} {params.prefix} $TMPDIR
        '''

rule KMC_dump:
    input:
        rules.KMC_count.output
    output:
        'kmers/{sample}.kmers'
    params:
        prefix = lambda wildcards, input: PurePath(input[0]).with_suffix('')
    localrule: True
    shell:
        '''
        kmc_dump {params.prefix} {output}
        '''

rule SRF:
    input:
        rules.KMC_dump.output
    output:
        'satellites/{sample}.srf.fa'
    resources:
        mem_mb = 2500
    shell:
        '''
        srf -p {wildcards.sample} {input} > {output}
        '''

rule minimap2_align:
    input:
        satellites = rules.SRF.output,
        reads = rules.samtools_fastq.output
    output:
        'satellites/{sample}.realign.paf'
    threads: 16
    resources:
        mem_mb = 3000,
        walltime = '24h'
    shell:
        '''
        minimap2 -c -t {threads} -N1000000 -f1000 -r100,100 <(srfutils.js enlong {input.satellites}) {input.reads} > {output}
        '''
 
rule srfutils:
    input:
        rules.minimap2_align.output
    output:
        bed = 'satellites/{sample}.abundance.bed',
        abun = 'satellites/{sample}.abundance.txt'
    localrule: True
    shell:
        '''
        srfutils.js paf2bed -l 1000 {input} > {output.bed}   # filter and extract non-overlapping regions
        srfutils.js bed2abun -g 2700000000 {output.bed} > {output.abun}
        '''

