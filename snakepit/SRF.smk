

rule KMC_count:
    input:
        ''
    output:
        ''
    shell:
        '''
        kmc -fq -k 151 -t {threads} -ci 100 -cs 100000 @fofn.txt {params.prefix} $TMPDIR
        '''

rule KMC_dump:
    input:
        rules.KMC_count.output
    output:
        ''
    params:
        prefix = ''
    localrule: True
    shell:
        '''
        kmc_dump {params.prefix} {output}
        '''

rule SRF:
    input:
        rules.KMC_dump.output
    output:
        ''
    shell:
        '''
        srf -p {params.prefix} {input} > {output}
        '''

rule minimap2_align:
    input:
        satellites = '',
        reads = ''
    output:
        ''
    threads: 4
    shell:
        '''
        minimap2 -c -t {threads} -N1000000 -f1000 -r100,100 <(./srfutils.js enlong {input.satellites}) {input.reads} > {output}
        '''
 
rule srfutils:
    input:
        rules.minimap2_align.output
    output:
        ''
    shell:
        '''
        srfutils.js paf2bed {input} > {output.bed}   # filter and extract non-overlapping regions
        srfutils.js bed2abun {output.bed} > {output.len}
        '''

