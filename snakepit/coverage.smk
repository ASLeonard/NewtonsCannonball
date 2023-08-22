
ruleorder: genome_coverage > region_coverage

rule genome_coverage:
    input:
        bam = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long_reads/{sample}/alignment/{sample}.mm2.cram',
        reference = config['reference']
    output:
        'coverage/{sample}.genome.csv'
    threads: 1
    resources:
        mem_mb = 25000
    shell:
        '''
        samtools coverage --reference {input.reference} {input.bam} |\
        grep -P "^(\d|X|Y|MT)" |\
        awk '{{print "{wildcards.sample}",$1,$6,$7,$9}}' > {output}
        '''

rule gather_genome:
    input:
        expand(rules.genome_coverage.output,sample=config['samples'])
    output:
        'coverage/genome.csv'
    localrule: True
    shell:
        '''
        echo "sample chromosome coverage meandepth meanmapq" > {output}
        cat {input} >> {output}
        '''

rule region_coverage:
    input:
        bam = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long_reads/{sample}/alignment/{sample}.mm2.cram',
        reference = config['reference']
    output:
        'coverage/{sample}.{region}.{filtering}.csv'
    params:
        region = lambda wildcards: config['regions'][wildcards.region],
        flags = lambda wildcards: {'default':'0','secondary':'256'}[wildcards.filtering]
    threads: 1
    resources:
        mem_mb = 15000
    shell:
        '''
        samtools bedcov -c -g {params.flags} --reference {input.reference} <(echo {params.region} |\
        sed s'/[^0-9]/\\t/g') {input.bam} |\
        awk '{{print "{wildcards.sample}","{wildcards.region}","{wildcards.filtering}",$4,$5}}' > {output}
        '''

rule gather_region:
    input:
        expand(rules.region_coverage.output,sample=config['samples'],region=config['regions'],filtering=('default','secondary'))
    output:
        'coverage/regions.csv'
    localrule: True
    shell:
        '''
        echo "sample region filtering bases reads" > {output}
        cat {input} >> {output}
        '''
