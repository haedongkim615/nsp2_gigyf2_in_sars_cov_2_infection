#!/usr/local/bin/snakemake

import pandas as pd

####################################################################################################

def detect_quality_offset(infile_path):
    import gzip

    def read_file(infile_path):
        if infile_path.endswith('gz'):
            for line in gzip.open(infile_path, 'rb'):
                line = line.decode('ascii').strip('\n')
                yield line
        else:
            for line in open(infile_path):
                yield line

    n = 0
    for line in read_file(infile_path):
        n += 1
        if n % 4 != 0: continue

        for q in line:
            if ord(q) < ord(';'): return '33' # Q33
            elif ord(q) > ord('J'): raise ValueError('Not Q33') #return '64' # Q64

####################################################################################################

info = pd.read_table('labels.txt').fillna('').astype(str)
for col in ['filepath_r1', 'filepath_r2', 'label', '5_adapter', '3_adapter']:
    if (info[col] == '').any().any(): raise ValueError('Blank info')
    if info[col].str.startswith('#').any(): raise ValueError('Comments should be removed', col)
for col in ['filepath_r1', 'filepath_r2', 'label']:
    if info[col].duplicated().any(): raise ValueError('Duplicated labels', col)
info = info[['filepath_r1', 'filepath_r2', 'label', '5_adapter', '3_adapter']]
info = info.set_index('label')
LABELS = info.T.to_dict()

for index, row in pd.read_table('global_variables.txt').iterrows():
    HOST = row['host']
    VIRUS = row['virus']
    STRAND_WITH_STRANDNESS = row['strand_with_strandness']
    READ_LENGTH = row['read_length']
    MIN_READ_LEN = row['min_read_len']

####################################################################################################

REF_DIR = '/casa/haedong/ref'

GENOME = {
    'grch38': '{REF_DIR}/grch38/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_genomic.fna.gz'.format(REF_DIR=REF_DIR), \
    'sars_cov_2': '{REF_DIR}/sars_cov_2/GCF_009858895.2/GCF_009858895.2_ASM985889v3_genomic.fna.gz'.format(REF_DIR=REF_DIR), \
}

ANNOTATION = {
    'grch38': '{REF_DIR}/grch38/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz'.format(REF_DIR=REF_DIR), \
    'sars_cov_2': '{REF_DIR}/sars_cov_2/GCF_009858895.2/GCF_009858895.2_ASM985889v3_genomic.only_cds.no_orf10.gtf.gz'.format(REF_DIR=REF_DIR), \
}

VIRUS_CHR = {
    'sars_cov_2': 'NC_045512.2', \
}[VIRUS]

####################################################################################################

rule all:
    input: expand('fx_files/{label}.{r}.fastq.gz', label=LABELS, r=['r1', 'r2']), \
           expand('bam_files/{label}.Aligned.sortedByCoord.out.unique.bam', label=LABELS), \
           expand('bam_files/{label}.Aligned.sortedByCoord.out.unique.bam.bai', label=LABELS), \
           expand('coverage_files/{label}.{VIRUS}.coverage.txt', label=LABELS, VIRUS=VIRUS), \
           expand('anno_files/{label}.bed.no_rrna_trna.strand_specific.gz', label=LABELS), \
           expand('piranha_files/{label}.window_{window}.bed.gz', label=[label for label in LABELS if 'antigigyf2' in label], window=[50]), \
           expand('piranha_files/{label}.window_{window}.bed.count.txt', label=[label for label in LABELS if 'antigigyf2' in label], window=[50])

####################################################################################################

rule uncompress_raw_fastq_file:
    input: r1 = lambda wildcards: LABELS[wildcards.label]['filepath_r1'], \
           r2 = lambda wildcards: LABELS[wildcards.label]['filepath_r2']
    output: r1 = temp('raw_files/{label}.r1.fq'), \
            r2 = temp('raw_files/{label}.r2.fq')
    threads: 16
    priority: 50
    run:
        shell('pigz -d -k -c -p {threads} {input.r1} > {output.r1}')
        shell('pigz -d -k -c -p {threads} {input.r2} > {output.r2}')

rule remove_5_adapter_and_trim_5_random_sequence:
    input: r1 = 'raw_files/{label}.r1.fq', \
           r2 = 'raw_files/{label}.r2.fq'
    output: r1 = temp('fx_files/{label}.r1.5_adapter.5_random.fq'), \
            r2 = temp('fx_files/{label}.r2.5_adapter.5_random.fq')
    params: r1 = 'fx_files/{label}.r1.5_adapter.fq', \
            r2 = 'fx_files/{label}.r2.5_adapter.fq'
    threads: 4
    run:
        from Bio.Seq import Seq

        label = wildcards.label
        phred = detect_quality_offset(input[0])
        if STRAND_WITH_STRANDNESS == 'forward':
            adapter_r1 = LABELS[label]['5_adapter']
            adapter_r2 = LABELS[label]['3_adapter']
            adapter_r2 = str(Seq(adapter_r2).reverse_complement())
        elif STRAND_WITH_STRANDNESS == 'reverse':
            adapter_r1 = LABELS[label]['3_adapter']
            adapter_r1 = str(Seq(adapter_r1).reverse_complement())
            adapter_r2 = LABELS[label]['5_adapter']
        min_read_len = MIN_READ_LEN + 4

        shell('cutadapt -g {adapter_r1} -G {adapter_r2} -m {min_read_len} -o {params.r1} -p {params.r2} {input.r1} {input.r2}')
        shell('cat {params.r1} | fastx_trimmer -Q {phred} -f 5 > {output.r1}')
        shell('cat {params.r2} | fastx_trimmer -Q {phred} -f 5 > {output.r2}')
        shell('rm -f {params.r1} {params.r2}')

rule remove_3_adapter_and_trim_3_random_sequence_from_r1:
    input: r1 = 'fx_files/{label}.r1.5_adapter.5_random.fq', \
           r2 = 'fx_files/{label}.r2.5_adapter.5_random.fq'
    output: r1 = temp('fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.fq'), \
            r2 = temp('fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.fq')
    params: r1 = 'fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.fq'
    threads: 4
    run:
        from Bio.Seq import Seq

        label = wildcards.label
        phred = detect_quality_offset(input[0])
        if STRAND_WITH_STRANDNESS == 'forward':
            adapter_r1 = LABELS[label]['3_adapter']
        elif STRAND_WITH_STRANDNESS == 'reverse':
            adapter_r1 = LABELS[label]['5_adapter']
            adapter_r1 = str(Seq(adapter_r1).reverse_complement())

        min_read_len = MIN_READ_LEN + 4

        shell('cutadapt -a {adapter_r1} -m {min_read_len} --discard-untrimmed -o {params.r1} -p {output.r2} {input.r1} {input.r2}')
        shell('cat {params.r1} | fastx_trimmer -Q {phred} -t 4 > {output.r1}')
        shell('rm -f {params.r1}')

rule extract_reads_which_dont_have_3_adapter_from_r1:
    input: r1 = 'fx_files/{label}.r1.5_adapter.5_random.fq', \
           r2 = 'fx_files/{label}.r2.5_adapter.5_random.fq'
    output: r1 = temp('fx_files/{label}.r1.5_adapter.5_random.no_3_adapter_of_r1.fq'), \
            r2 = temp('fx_files/{label}.r2.5_adapter.5_random.no_3_adapter_of_r1.fq')
    threads: 4
    run:
        from Bio.Seq import Seq

        label = wildcards.label
        phred = detect_quality_offset(input[0])
        if STRAND_WITH_STRANDNESS == 'forward':
            adapter_r1 = LABELS[label]['3_adapter']
        elif STRAND_WITH_STRANDNESS == 'reverse':
            adapter_r1 = LABELS[label]['5_adapter']
            adapter_r1 = str(Seq(adapter_r1).reverse_complement())

        min_read_len = MIN_READ_LEN + 4

        shell('cutadapt -a {adapter_r1} -m {min_read_len} --discard-trimmed -o {output.r1} -p {output.r2} {input.r1} {input.r2}')

rule merge_3_adapter_clipped_files_from_r1:
    input: r1_3_adapter = 'fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.fq', \
           r2_3_adapter = 'fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.fq', \
           r1_no_3_adapter = 'fx_files/{label}.r1.5_adapter.5_random.no_3_adapter_of_r1.fq', \
           r2_no_3_adapter = 'fx_files/{label}.r2.5_adapter.5_random.no_3_adapter_of_r1.fq'
    output: r1 = temp('fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.merged.fq'), \
            r2 = temp('fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.merged.fq')
    threads: 4
    run:
        shell('cat {input.r1_3_adapter} {input.r1_no_3_adapter} > {output.r1}')
        shell('cat {input.r2_3_adapter} {input.r2_no_3_adapter} > {output.r2}')

# cutadapt does not discard adapter-containing reads from read 2
# thus, read 2 should be considered as read 1
rule remove_3_adapter_and_trim_3_random_sequence_from_r2:
    input: r1 = 'fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.merged.fq', \
           r2 = 'fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.merged.fq'
    output: r1 = temp('fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.merged.3_adapter_of_r2.3_random.fq'), \
            r2 = temp('fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.merged.3_adapter_of_r2.3_random.fq')
    params: r2 = 'fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.merged.3_adapter_of_r2.fq'
    threads: 4
    run:
        from Bio.Seq import Seq

        label = wildcards.label
        phred = detect_quality_offset(input[0])
        if STRAND_WITH_STRANDNESS == 'forward':
            adapter_r2 = LABELS[label]['5_adapter']
            adapter_r2 = str(Seq(adapter_r2).reverse_complement())
        elif STRAND_WITH_STRANDNESS == 'reverse':
            adapter_r2 = LABELS[label]['3_adapter']

        min_read_len = MIN_READ_LEN + 4

        shell('cutadapt -a {adapter_r2} -m {min_read_len} --discard-untrimmed -o {params.r2} -p {output.r1} {input.r2} {input.r1}')
        shell('cat {params.r2} | fastx_trimmer -Q {phred} -t 4 > {output.r2}')
        shell('rm -f {params.r2}')

rule extract_reads_which_dont_have_3_adapter_from_r2:
    input: r1 = 'fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.merged.fq', \
           r2 = 'fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.merged.fq'
    output: r1 = temp('fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.merged.no_3_adapter_of_r2.fq'), \
            r2 = temp('fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.merged.no_3_adapter_of_r2.fq')
    threads: 4
    run:
        from Bio.Seq import Seq

        label = wildcards.label
        phred = detect_quality_offset(input[0])
        if STRAND_WITH_STRANDNESS == 'forward':
            adapter_r2 = LABELS[label]['5_adapter']
            adapter_r2 = str(Seq(adapter_r2).reverse_complement())
        elif STRAND_WITH_STRANDNESS == 'reverse':
            adapter_r2 = LABELS[label]['3_adapter']

        min_read_len = MIN_READ_LEN + 4

        shell('cutadapt -a {adapter_r2} -m {min_read_len} --discard-trimmed -o {output.r2} -p {output.r1} {input.r2} {input.r1}')

rule merge_3_adapter_clipped_files_from_r1_and_r2:
    input: r1_3_adapter = 'fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.merged.3_adapter_of_r2.3_random.fq', \
           r2_3_adapter = 'fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.merged.3_adapter_of_r2.3_random.fq', \
           r1_no_3_adapter = 'fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.merged.no_3_adapter_of_r2.fq', \
           r2_no_3_adapter = 'fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.merged.no_3_adapter_of_r2.fq'
    output: r1 = temp('fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.merged.3_adapter_of_r2.3_random.merged.fq.gz'), \
            r2 = temp('fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.merged.3_adapter_of_r2.3_random.merged.fq.gz')
    threads: 4
    run:
        shell('cat {input.r1_3_adapter} {input.r1_no_3_adapter} | pigz -p {threads} > {output.r1}')
        shell('cat {input.r2_3_adapter} {input.r2_no_3_adapter} | pigz -p {threads} > {output.r2}')

####################################################################################################

rule remove_adapter_multimer:
    input: r1 = 'fx_files/{label}.r1.5_adapter.5_random.3_adapter_of_r1.3_random.merged.3_adapter_of_r2.3_random.merged.fq.gz', \
           r2 = 'fx_files/{label}.r2.5_adapter.5_random.3_adapter_of_r1.3_random.merged.3_adapter_of_r2.3_random.merged.fq.gz', \
           idx = '{REF_DIR}/truseq_adapter/aligner/bowtie2/genome.1.bt2'.format(REF_DIR=REF_DIR)
    output: temp('adapter_multimer_files/{label}.sam'), \
            'fx_files/{label}.r1.fastq.gz', \
            'fx_files/{label}.r2.fastq.gz'
    threads: 32
    params: idx_dir = '{REF_DIR}/truseq_adapter/aligner/bowtie2/genome'.format(REF_DIR=REF_DIR), \
            out_pfx = 'adapter_multimer_files/{label}'
    run:
        shell('bowtie2 -x {params.idx_dir} -1 {input.r1} -2 {input.r2} -S {params.out_pfx}.sam --threads {threads} --local --un-conc fx_files/{wildcards.label}')
        shell('mv fx_files/{wildcards.label}.1 fx_files/{wildcards.label}.r1.fastq')
        shell('mv fx_files/{wildcards.label}.2 fx_files/{wildcards.label}.r2.fastq')
        shell('pigz -p {threads} fx_files/{wildcards.label}.r1.fastq')
        shell('pigz -p {threads} fx_files/{wildcards.label}.r2.fastq')

rule extract_adapter_multimer:
    input: 'adapter_multimer_files/{label}.sam'
    output: 'adapter_multimer_files/{label}.bam'
    threads: 32
    shell: 'samtools view -bS -F 4 -@ {threads} {input} > {output}'

####################################################################################################

rule merge_genome:
    input: host = GENOME[HOST],
           virus = GENOME[VIRUS]
    output: '{REF_DIR}/{HOST}_{VIRUS}/genomic.fna.gz'
    threads: 24
    shell: 'zcat {input.host} {input.virus} | pigz -p {threads} > {output}'

rule merge_annotation:
    input: host = ANNOTATION[HOST],
           virus = ANNOTATION[VIRUS]
    output: '{REF_DIR}/{HOST}_{VIRUS}/genomic.gtf.gz'
    threads: 8
    shell: 'cat <(zcat {input.host} | tail -n+6 | head -n -1) <(zcat {input.virus} | tail -n+4) | \
            pigz -p {threads} > {output}'

####################################################################################################

rule uncompress_genome_file:
    input: '{filepath}.fna.gz'
    output: temp('{filepath}.fna')
    threads: 16
    shell: 'pigz -d -k -p {threads} {input}'

rule uncompress_annotation_file:
    input: '{filepath}.gtf.gz'
    output: temp('{filepath}.gtf')
    threads: 16
    shell: 'pigz -d -k -p {threads} {input}'

####################################################################################################

rule build_STAR_index_with_sjdbGTFfile:
    input: fa = '{REF_DIR}/{HOST}_{VIRUS}/genomic.fna', \
           gtf = '{REF_DIR}/{HOST}_{VIRUS}/genomic.gtf'
    output: '{{REF_DIR}}/{{HOST}}_{{VIRUS}}/aligner/star_sjdboverhang{LEN}/Genome'.format(LEN=READ_LENGTH-1)
    threads: 32
    params: idx_dir = '{{REF_DIR}}/{{HOST}}_{{VIRUS}}/aligner/star_sjdboverhang{LEN}/'.format(LEN=READ_LENGTH-1)
    run:
        sjdboverhang = READ_LENGTH - 1
        shell('STAR --runMode genomeGenerate --genomeDir {params.idx_dir} --runThreadN {threads} \
               --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang {sjdboverhang}')

rule align_using_STAR_with_hyeshik_options:
    input: r1 = 'fx_files/{label}.r1.fastq.gz', \
           r2 = 'fx_files/{label}.r2.fastq.gz', \
           idx = '{REF_DIR}/{HOST}_{VIRUS}/aligner/star_sjdboverhang{LEN}/Genome'.format(REF_DIR=REF_DIR, HOST=HOST, VIRUS=VIRUS, LEN=READ_LENGTH-1)
    output: 'bam_files/{label}.Aligned.sortedByCoord.out.bam', \
            temp('bam_files/{label}.Aligned.toTranscriptome.out.bam'), \
            'bam_files/{label}.Unmapped.out.mate1.gz', \
            'bam_files/{label}.Unmapped.out.mate2.gz'
    threads: 32
    params: idx_dir = '{REF_DIR}/{HOST}_{VIRUS}/aligner/star_sjdboverhang{LEN}/'.format(REF_DIR=REF_DIR, HOST=HOST, VIRUS=VIRUS, LEN=READ_LENGTH-1), \
            out_pfx = 'bam_files/{label}.'
    run:
        shell('STAR --genomeDir {params.idx_dir} \
                    --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat \
                    --runThreadN {threads} --outFileNamePrefix {params.out_pfx} \
                    --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx \
                    --outSAMunmapped Within KeepPairs --quantMode TranscriptomeSAM GeneCounts \
                    --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
                    --alignSJDBoverhangMin 1 --outSJfilterOverhangMin 12 12 12 12 \
                    --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 \
                    --outSJfilterDistToOtherSJmin 0 0 0 0 --outFilterMismatchNmax 999 \
                    --outFilterMismatchNoverReadLmax 0.04 --scoreGapNoncan -4 --scoreGapATAC -4 \
                    --chimOutType Junctions WithinBAM HardClip --chimScoreJunctionNonGTAG 0 \
                    --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignIntronMin 20 \
                    --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20')
        shell('pigz -p {threads} {params.out_pfx}Unmapped.out.mate1')
        shell('pigz -p {threads} {params.out_pfx}Unmapped.out.mate2')

rule extract_unique_mapping:
    input: 'bam_files/{filename}.bam'
    output: 'bam_files/{filename}.unique.bam'
    shell: 'samtools view -@ {threads} -b -q 4 -o {output} {input}'

####################################################################################################

rule index_bam:
    input: '{filepath}.bam'
    output: '{filepath}.bam.bai'
    threads: 32
    shell: 'samtools index -b -@ {threads} {input}'

rule extract_viral_mappings:
    input: bam = 'bam_files/{label}.Aligned.sortedByCoord.out.unique.bam', \
           bai = 'bam_files/{label}.Aligned.sortedByCoord.out.unique.bam.bai'
    output: temp('bam_files/{label}.Aligned.sortedByCoord.out.unique.{VIRUS}.bam')
    threads: 32
    shell: 'samtools view -@ {threads} -b -q 4 -o {output} {input.bam} {VIRUS_CHR}'

rule calculate_coverage:
    input: 'bam_files/{label}.Aligned.sortedByCoord.out.unique.{VIRUS}.bam'
    output: 'coverage_files/{label}.{VIRUS}.coverage.txt'
    shell: 'bedtools genomecov -ibam {input} -dz -split > {output}'

####################################################################################################

rule make_annotation_of_transcript_unit:
    input: host = ANNOTATION[HOST],
           virus = GENOME[VIRUS]
    output: '{REF_DIR}/{HOST}_{VIRUS}/genomic.transcript_unit.gtf.gz'.format(REF_DIR=REF_DIR, HOST=HOST, VIRUS=VIRUS)
    params: '{REF_DIR}/{HOST}_{VIRUS}/genomic.transcript_unit.gtf'.format(REF_DIR=REF_DIR, HOST=HOST, VIRUS=VIRUS)
    run:
        host = pd.read_table(input.host, header=None, skiprows=5, skipfooter=1)
        host = host.loc[host.iloc[:, 2] == 'gene'].drop_duplicates()

        import gzip
        from Bio import SeqIO

        with gzip.open(input.virus, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                chrom = record.id
                genome_length = len(record.seq)
        virus_pos = pd.DataFrame([[chrom, 'RefSeq', 'gene', 1, genome_length, '.', '+', '.', 'gene_id "SARS-CoV-2 +"; transcript_id ""; gene_biotype "protein_coding";']])
        virus_neg = pd.DataFrame([[chrom, 'RefSeq', 'gene', 1, genome_length, '.', '-', '.', 'gene_id "SARS-CoV-2 -"; transcript_id ""; gene_biotype "lncRNA";']])

        anno = pd.concat([host, virus_pos, virus_neg])
        import csv
        anno.to_csv(params[0], sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)
        shell('pigz {params}')

rule bedtools_intersect_with_same_strandness:
    input: bam = 'bam_files/{label}.Aligned.sortedByCoord.out.unique.bam', \
           anno = '{REF_DIR}/{HOST}_{VIRUS}/genomic.transcript_unit.gtf.gz'.format(REF_DIR=REF_DIR, HOST=HOST, VIRUS=VIRUS)
    output: 'anno_files/{label}.bed.gz'
    threads: 4
    shell: 'bedtools intersect -abam {input.bam} -b {input.anno} -wa -wb -f 0.50 -s -bed | \
            uniq | pigz -p {threads} > {output}'

rule remove_rrna_trna_and_select_stranded_information_from_bed:
    input: 'anno_files/{label}.bed.gz'
    output: 'anno_files/{label}.bed.no_rrna_trna.strand_specific.gz'
    params: 'anno_files/{label}.bed.no_rrna_trna.strand_specific'
    threads: 4
    run:
        if STRAND_WITH_STRANDNESS == 'forward': read_end = '/1'
        else: read_end = '/2'

        read_ids = []
        chunksize = 10 ** 6
        chunks = []
        for chunk in pd.read_table(input[0], header=None, chunksize=chunksize):
            chunk = chunk.loc[chunk.iloc[:, 3].str.endswith(read_end)]
            chunk.loc[:, 'gene_biotype'] = chunk.iloc[:, 20].str.split('gene_biotype').str.get(1).astype(str).str.split('"').str.get(1).fillna('')
            chunk.loc[:, 'description'] = chunk.iloc[:, 20].str.split('description').str.get(1).astype(str).str.split('"').str.get(1).fillna('')

            is_rrna_trna = chunk['gene_biotype'].isin(('rRNA', 'tRNA'))
            is_rrna_pseudogene = ((chunk['gene_biotype'] == 'pseudogene') & chunk['description'].str.contains('ribosomal pseudogene'))
            is_trna_pseudogene = ((chunk['gene_biotype'] == 'pseudogene') & chunk['description'].str.startswith('tRNA-') & (~chunk['description'].str.contains('protein')))

            read_ids += chunk.loc[is_rrna_trna | is_rrna_pseudogene | is_trna_pseudogene].iloc[:, 3].unique().tolist()
            read_ids = list(set(read_ids))

            chunk = chunk.loc[~(is_rrna_trna | is_rrna_pseudogene | is_trna_pseudogene)]
            chunks.append(chunk)

        bed = pd.concat(chunks)
        bed = bed.loc[~bed.iloc[:, 3].isin(read_ids)]
        bed = bed.loc[:, bed.columns.drop(['gene_biotype', 'description'])]
        import csv
        bed.to_csv(params[0], sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)

        shell('pigz -p {threads} {params}')

#####################################################################################################

rule divide_annotation_into_utr_cds:
    input: host = ANNOTATION[HOST], \
           virus_genome = GENOME[VIRUS], \
           virus_annotation = ANNOTATION[VIRUS]
    output: '{REF_DIR}/{HOST}_{VIRUS}/genomic.transcript_unit.non_redundant.utr_cds_divided.gtf.gz'.format(REF_DIR=REF_DIR, HOST=HOST, VIRUS=VIRUS)
    params: '{REF_DIR}/{HOST}_{VIRUS}/genomic.transcript_unit.non_redundant.utr_cds_divided.gtf'.format(REF_DIR=REF_DIR, HOST=HOST, VIRUS=VIRUS)
    run:
        host = pd.read_table(input.host, header=None, skiprows=5, skipfooter=1)
        host.loc[:, 'name'] = host.iloc[:, 8].str.split('gene_id').str.get(1).str.split('"').str.get(1)
        host.loc[:, 'transcript_id'] = host.iloc[:, 8].str.split('transcript_id').str.get(1).str.split('"').str.get(1)

        def find_maximum_length_transcript_id(transcript_ids, transcript2length):
            max_val = 0
            for transcript_id in transcript_ids:
                if transcript2length[transcript_id] > max_val:
                    max_key = transcript_id
                    max_val = transcript2length[transcript_id]
            return max_key

        from collections import defaultdict

        name2transcripts = defaultdict(list)
        name = ''
        for index, row in host.iterrows():
            if row[2] == 'gene':
                if (name != '') and (not was_gene):
                    if gene_biotype == 'protein_coding':
                        candidates_nm, candidates_xm = [], []
                        for transcript_id, has_cds in transcript2has_cds.items():
                            if has_cds and transcript_id.startswith('NM'): candidates_nm.append(transcript_id)
                            elif has_cds and transcript_id.startswith('XM'): candidates_xm.append(transcript_id)

                        transcript_ids_nm = []
                        for transcript_id in candidates_nm:
                            if (left in transcript2left[transcript_id]) and (right in transcript2right[transcript_id]): transcript_ids_nm.append(transcript_id)


                        transcript_ids_xm = []
                        for transcript_id in candidates_xm:
                            if (left in transcript2left[transcript_id]) and (right in transcript2right[transcript_id]): transcript_ids_xm.append(transcript_id)

                        if len(transcript_ids_nm) == 1: name2transcripts[name].append(transcript_ids_nm[0])
                        elif len(transcript_ids_nm) > 1: name2transcripts[name].append(find_maximum_length_transcript_id(transcript_ids_nm, transcript2length))
                        elif len(transcript_ids_xm) == 1: name2transcripts[name].append(transcript_ids_xm[0])
                        elif len(transcript_ids_xm) > 1: name2transcripts[name].append(find_maximum_length_transcript_id(transcript_ids_xm, transcript2length))
                        else:
                            if len(candidates_nm) >= 1: name2transcripts[name].append(find_maximum_length_transcript_id(candidates_nm, transcript2length))
                            elif len(candidates_xm) >= 1: name2transcripts[name].append(find_maximum_length_transcript_id(candidates_xm, transcript2length))
                            else: raise ValueError()

                    else:
                        transcript_ids = []
                        for transcript_id in transcript2length:
                            if (left in transcript2left[transcript_id]) and (right in transcript2right[transcript_id]): transcript_ids.append(transcript_id)

                        if len(transcript_ids) == 1: name2transcripts[name].append(transcript_ids[0])
                        elif len(transcript_ids) > 1: name2transcripts[name].append(find_maximum_length_transcript_id(transcript_ids, transcript2length))
                        else: name2transcripts[name].append(max(transcript2length, key=transcript2length.get))

                chrom, source, feature, left, right, score, strand, frame, attribute, name, transcript_id = row
                gene_biotype = attribute.split('gene_biotype')[1].split('"')[1]
                transcript2length = defaultdict(int)
                transcript2left = defaultdict(list)
                transcript2right = defaultdict(list)
                transcript2has_cds = {}
                was_gene = True

            if row[2] == 'exon':
                transcript2length[row['transcript_id']] += (row[4] - row[3])
                transcript2left[row['transcript_id']].append(row[3])
                transcript2right[row['transcript_id']].append(row[4])
                transcript2has_cds[row['transcript_id']] = False
                was_gene = False

            elif (row[2] == 'start_codon') or (row[2] == 'stop_codon'): transcript2has_cds[row['transcript_id']] = True

        for key, val in name2transcripts.items():
            if len(val) > 1: raise ValueError()


        outfile = open(params[0], 'w')
        name = ''
        for index, row in host.iterrows():
            if row[2] == 'gene':
                if name != '':
                    if was_gene: print(*[chrom, source, 'gene', left, right, score, strand, frame, 'gene_id "{}"; transcript_id ""; gene_biotype "{}";'.format(name, gene_biotype)], sep='\t', file=outfile)
                    else:
                        positions = []
                        for exon_left, exon_right in exons:
                            positions += [exon_left, exon_right]
                        if gene_biotype == 'protein_coding': positions += [cds_left, cds_right]
                        if strand == '+': positions = sorted(positions)
                        elif strand == '-': positions = sorted(positions, reverse=True)

                        if gene_biotype != 'protein_coding':
                            for i in range(len(positions) - 1):
                                left, right = sorted([positions[i], positions[i + 1]])
                                is_exon = any([(exon_left <= left <= right <= exon_right) for exon_left, exon_right in exons])
                                if is_exon: print(*[chrom, source, 'exon', left, right, score, strand, frame, 'gene_id "{}"; transcript_id "{}"; gene_biotype "{}";'.format(name, name2transcripts[name][0], gene_biotype)], sep='\t', file=outfile)
                                else: print(*[chrom, source, 'intron', left + 1, right - 1, score, strand, frame, 'gene_id "{}"; transcript_id "{}"; gene_biotype "{}";'.format(name, name2transcripts[name][0], gene_biotype)], sep='\t', file=outfile)

                        else:
                            if strand == '+': utr_left, utr_right = '5UTR', '3UTR'
                            else: utr_left, utr_right = '3UTR', '5UTR'

                            for i in range(len(positions) - 1):
                                left, right = sorted([positions[i], positions[i + 1]])
                                is_exon = any([(exon_left <= left <= right <= exon_right) for exon_left, exon_right in exons])
                                if is_exon:
                                    if right < cds_left: print(*[chrom, source, utr_left, left, right, score, strand, frame, 'gene_id "{}"; transcript_id "{}"; gene_biotype "{}";'.format(name, name2transcripts[name][0], gene_biotype)], sep='\t', file=outfile)
                                    elif right == cds_left: print(*[chrom, source, utr_left, left, right - 1, score, strand, frame, 'gene_id "{}"; transcript_id "{}"; gene_biotype "{}";'.format(name, name2transcripts[name][0], gene_biotype)], sep='\t', file=outfile)
                                    elif cds_left <= left <= right <= cds_right: print(*[chrom, source, 'CDS', left, right, score, strand, frame, 'gene_id "{}"; transcript_id "{}"; gene_biotype "{}";'.format(name, name2transcripts[name][0], gene_biotype)], sep='\t', file=outfile)
                                    elif cds_right == left: print(*[chrom, source, utr_right, left + 1, right, score, strand, frame, 'gene_id "{}"; transcript_id "{}"; gene_biotype "{}";'.format(name, name2transcripts[name][0], gene_biotype)], sep='\t', file=outfile)
                                    elif cds_right < left: print(*[chrom, source, utr_right, left, right, score, strand, frame, 'gene_id "{}"; transcript_id "{}"; gene_biotype "{}";'.format(name, name2transcripts[name][0], gene_biotype)], sep='\t', file=outfile)
                                    else: raise ValueError()
                                else: print(*[chrom, source, 'intron', left + 1, right - 1, score, strand, frame, 'gene_id "{}"; transcript_id "{}"; gene_biotype "{}";'.format(name, name2transcripts[name][0], gene_biotype)], sep='\t', file=outfile)
                chrom, source, feature, left, right, score, strand, frame, attribute, name, transcript_id = row
                gene_biotype = attribute.split('gene_biotype')[1].split('"')[1]
                exons = []
                was_gene = True

            else:
                if row['name'] != name: raise ValueError()
                if row['transcript_id'] not in name2transcripts[row['name']]: continue

                if row[2] == 'exon': exons.append([row[3], row[4]])
                elif row[2] == 'start_codon':
                    if row[6] == '+': cds_left = row[3]
                    elif row[6] == '-': cds_right = row[4]
                    else: raise ValueError()
                elif row[2] == 'stop_codon':
                    if row[6] == '+': cds_right = row[4]
                    elif row[6] == '-': cds_left = row[3]
                    else: raise ValueError()

                was_gene = False


        import gzip
        from Bio import SeqIO

        with gzip.open(input.virus_genome, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                genome_length = len(record.seq)

        virus = pd.read_table(input.virus_annotation, header=None, skiprows=3, skipfooter=1)
        virus = virus.loc[virus.iloc[:, 2] == 'exon']

        lines = []
        prev_right = 0
        for index, row in virus.iterrows():
            chrom, source, feature, left, right, score, strand, frame, attribute = row
            if strand == '+': strand_opposite = '-'
            else: strand_opposite = '+'
            if prev_right < left:
                if prev_right == 0:
                    print(*[chrom, source, '5UTR', prev_right + 1, left - 1, score, strand, frame, 'gene_id "5UTR"; transcript_id ""; gene_biotype "protein_coding";'], sep='\t', file=outfile)
                    lines.append([chrom, source, 'gene', prev_right + 1, left - 1, score, strand_opposite, frame, 'gene_id "5UTR -"; transcript_id ""; gene_biotype "lncRNA";'])
                else:
                    print(*[chrom, source, 'UTR', prev_right + 1, left - 1, score, strand, frame, 'gene_id "UTR"; transcript_id ""; gene_biotype "protein_coding";'], sep='\t', file=outfile)
                    lines.append([chrom, source, 'gene', prev_right + 1, left - 1, score, strand_opposite, frame, 'gene_id "UTR -"; transcript_id ""; gene_biotype "lncRNA";'])

            print(*[chrom, source, 'CDS', left, right, score, strand, frame, 'gene_id "{}"; transcript_id ""; gene_biotype "protein_coding";'.format(attribute.split('gene_id')[1].split('"')[1])], sep='\t', file=outfile)
            lines.append([chrom, source, 'gene', left, right, score, strand_opposite, frame, 'gene_id "{} -"; transcript_id ""; gene_biotype "lncRNA";'.format(attribute.split('gene_id')[1].split('"')[1])])

            prev_right = right

        print(*[chrom, source, '3UTR', right + 1, genome_length, score, strand, frame, 'gene_id "3UTR"; transcript_id ""; gene_biotype "protein_coding";'], sep='\t', file=outfile)
        lines.append([chrom, source, 'gene', right + 1, genome_length, score, strand_opposite, frame, 'gene_id "3UTR -"; transcript_id ""; gene_biotype "lncRNA";'])

        for line in lines[::-1]:
            print(*line, sep='\t', file=outfile)

        outfile.close()
        shell('pigz {params}')

rule peak_enrichment_analysis:
    input: sample = 'anno_files/{label}.bed.no_rrna_trna.strand_specific.gz', \
           control = lambda wildcards: 'anno_files/{label_ctrl}.bed.no_rrna_trna.strand_specific.gz'.format(label_ctrl=wildcards.label.replace('antigigyf2', 'igg'))
    output: temp('piranha_files/{label}.window_{window}.unannotated.bed')
    shell: 'Piranha -s -z {wildcards.window} -n -l -o {output} <(zcat {input.sample}) <(zcat {input.control})'

rule bedtools_intersect_to_get_gene_annotation_of_the_peak:
    input: bed = 'piranha_files/{label}.window_{window}.unannotated.bed', \
           anno = '{REF_DIR}/{HOST}_{VIRUS}/genomic.transcript_unit.non_redundant.utr_cds_divided.gtf.gz'.format(REF_DIR=REF_DIR, HOST=HOST, VIRUS=VIRUS)
    output: temp('piranha_files/{label}.window_{window}.annotated.bed')
    shell: 'bedtools intersect -a {input.bed} -b {input.anno} -wa -wb -f 0.50 -s | uniq > {output}'

rule sort_by_pval:
    input: 'piranha_files/{label}.window_{window}.annotated.bed'
    output: 'piranha_files/{label}.window_{window}.bed.gz'
    threads: 4
    params: 'piranha_files/{label}.window_{window}.bed'
    run:
        bed = pd.read_table(input[0], header=None)
        bed = bed.sort_values([6, 0, 1])
        bed.loc[:, 'type'] = [row[16].split('gene_biotype')[1].split('"')[1] if row[0] != 'NC_045512.2' else 'SARS-CoV-2' for index, row in bed.iterrows()]
        bed.loc[:, 'name'] = bed[16].str.split('gene_id').str.get(1).str.split('"').str.get(1)
        bed = bed.loc[:, [0, 1, 2, 3, 4, 5, 6, 7, 'type', 'name', 10]]
        bed.columns = ['chrom', 'start', 'end', 'x', 'score_of_bed', 'strand', 'adj_pval', 'covariate', 'type', 'name', 'region']
        bed.to_csv(params[0], sep='\t', header=True, index=False)

        shell('pigz -p {threads} {params}')

rule bedtools_intersect_to_get_read_count_of_the_peak:
    input: bed1 = 'piranha_files/{label1}.window_{window}.unannotated.bed', \
           bed2 = 'anno_files/{label2}.bed.no_rrna_trna.strand_specific.gz'
    output: temp('piranha_files/{label1}.window_{window}.unannotated.{label2}.bed.no_rrna_trna.strand_specific.bed')
    threads: 8
    shell: 'bedtools intersect -a {input.bed1} -b {input.bed2} -wa -wb -s | uniq > {output}'

rule get_read_count_of_the_peak:
    input: expand('piranha_files/{{label1}}.window_{{window}}.unannotated.{label2}.bed.no_rrna_trna.strand_specific.bed', label2=LABELS), \
           'piranha_files/{label1}.window_{window}.bed.gz'
    output: 'piranha_files/{label1}.window_{window}.bed.count.txt'
    threads: 8
    run:
        from collections import defaultdict

        dataframes = []
        for filepath in input[:-1]:
            label = filepath.split('.')[-5]
            bed = pd.read_table(filepath, header=None)

            # if this takes long, use 'chunk'
            peak2count = {}
            for key, grp in bed.groupby([0, 1, 2, 3, 5, 28]):
                peak2count[key] = len(grp.iloc[:, 11].unique())

            lines = []
            for key, val in peak2count.items():
                lines.append(list(key) + [val])
            df = pd.DataFrame(lines)
            df.columns = ['chrom', 'start', 'end', 'x', 'strand', 'description', label]

            df.loc[:, 'type'] = [row['description'].split('gene_biotype')[1].split('"')[1] if row['chrom'] != 'NC_045512.2' else 'SARS-CoV-2' for index, row in df.iterrows()]
            df.loc[:, 'name'] = df['description'].str.split('gene_id').str.get(1).str.split('"').str.get(1)
            df = df[df.columns.drop('description')]
            df = df.set_index(['chrom', 'start', 'end', 'x', 'strand', 'type', 'name'])
            dataframes.append(df)

        count = pd.concat(dataframes, axis=1).fillna(0).astype(int)
        count = count.reset_index()

        anno = pd.read_table(input[-1])
        count_host = pd.merge(count, anno, how='inner', on=['chrom', 'start', 'end', 'x', 'strand', 'type', 'name'])
        count_host = count_host[['chrom', 'start', 'end', 'x', 'score_of_bed', 'strand', 'adj_pval', 'covariate', 'type', 'name', 'region'] + \
                                count_host.columns.drop(['chrom', 'start', 'end', 'x', 'score_of_bed', 'strand', 'adj_pval', 'covariate', 'type', 'name', 'region']).tolist()]
        count_virus = pd.merge(count.loc[count['chrom'] == 'NC_045512.2'], anno.loc[anno['chrom'] == 'NC_045512.2'], how='inner', on=['chrom', 'start', 'end', 'x', 'strand', 'type'])
        count_virus = count_virus.loc[:, ~count_virus.columns.str.endswith('_x')]
        count_virus.columns = [col[:-2] if col.endswith('_y') else col for col in count_virus.columns]
        count_virus = count_virus[['chrom', 'start', 'end', 'x', 'score_of_bed', 'strand', 'adj_pval', 'covariate', 'type', 'name', 'region'] + \
                                  count_virus.columns.drop(['chrom', 'start', 'end', 'x', 'score_of_bed', 'strand', 'adj_pval', 'covariate', 'type', 'name', 'region']).tolist()]

        count = pd.concat([count_host, count_virus])
        count = count.sort_values(['adj_pval', 'chrom', 'start'])
        count.to_csv(output[0], sep='\t', header=True, index=False)
