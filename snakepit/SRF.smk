from pathlib import PurePath

rule all:
    input:
        expand('satellites/{sample}.realign.paf',sample=config['samples'])

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
        mem_mb = 9000
    shell:
        '''
        kmc -k151 -t{threads} -ci{params.threshold} -cs1000000000 -m20 -fa {input.reads} {params.prefix} $TMPDIR
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

rule TRF:
    input:
        rules.SRF.output
    output:
        TRF = 'satellites/{sample}.TRF',
        monomers = 'satellites/{sample}.monomers'
    params:
        min_span = 0.75
    shell:
        '''
        TRF {input} 2 5 7 80 10 50 2000 -h -ngs | tee {output.TRF} |\
        awk -v m=0 -v T={params.min_span} '{{ if ($1~/@/) {{C=substr($1,2);D=$2;split($1,Q,"-");S=Q[2]; if (m>0) {{print X;m=0 }} }} else {{ if ($8>m&&($2-$1)>(S*T)) {{m=$8; X=">"C" "D",start="$1",end="$2",score="$8"\\n"$14 }} }} }} END {{print X}}' > {output.monomers}
        #manually print last case since no new line to trigger it

        #awk '{{ if ($1~/@/) {{ C=substr($1,2);n=1;D=$2 }} else {{ if ($1~/[[:digit:]+]/) {{ print ">"C"_"n" "D",start="$1",end="$2",score="$8"\\n"$14; n+=1 }} }}    }}' > {output.monomers}
        '''

rule curate_satellite_repeats:
    input:
       SRF = rules.SRF.output,
       TRF = rules.TRF.output
    output:
        'satellites/{sample}.txt'
    shell:
        '''
        samtools faidx -r <(grep -vFf <(awk '/>/ {{ print substr($1,2) }}' {input.TRF[1]}) <(awk '/>/ {{ print substr($1,2) }}' {input.SRF})) {input.SRF} | cat {input.TRF[1]} - |\
        python -c '
import sys;
def minimal_string(string):
    return sorted([string[i:]+string[:i] for i in range(len(string))])[0];
for line in sys.stdin:
    if line[0] == ">":
        print(line.rstrip());
    else:
        print(minimal_string(line.rstrip()))
        ' > {output}
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
        walltime = '4h'
    shell:
        '''
        minimap2 -c -t {threads} -N1000000 -f1000 -r100,100 <(srfutils.js enlong {input.satellites}) {input.reads} > {output}
        '''
 
rule srfutils:
    input:
        paf = rules.minimap2_align.output,
        coverage = rules.estimate_coverage.output
    output:
        bed = 'satellites/{sample}.abundance.bed',
        abun = 'satellites/{sample}.abundance.txt'
    params:
        size = lambda wildcards, input: int(float(open(input.coverage[0]).read())*2759153975)
    localrule: True
    shell:
        '''
        srfutils.js paf2bed -l 1000 {input.paf} > {output.bed}   # filter and extract non-overlapping regions
        srfutils.js bed2abun -g {params.size} {output.bed} > {output.abun}
        '''

