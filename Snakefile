import re
# pull file names and paths
SAMPLE_FASTQ_TABLE = open(config['sample_fastq_table']) # tab separated text file (no column names) with sample name and fastq file 
SAMPLE_DICT = dict()
for line in SAMPLE_FASTQ_TABLE:
	sample = line.split()[0]
	fastq = line.split()[1]
	if config['forward_suffix'] in fastq and sample not in SAMPLE_DICT:
		SAMPLE_DICT[sample] = {'Forward' : [fastq], 'Reverse' : []}
	elif config['reverse_suffix'] in fastq and sample not in SAMPLE_DICT:
		SAMPLE_DICT[sample] = {'Forward' : [], 'Reverse': [fastq]}
	elif config['forward_suffix'] in fastq and sample in SAMPLE_DICT:
		SAMPLE_DICT[sample]['Forward'].append(fastq)
	elif config['reverse_suffix'] in fastq and sample in SAMPLE_DICT:
		SAMPLE_DICT[sample]['Reverse'].append(fastq)
	else:
		print("parsing error?")


def fastq_by_sample(sample, read):
	if read == 'Forward':
		out = ['fastq_trimmed/' + i for i in SAMPLE_DICT[sample]['Forward']]
	elif read == 'Reverse':
		out = ['fastq_trimmed/' + i for i in SAMPLE_DICT[sample]['Reverse']]
	else:
		print("missing " + sample + ' ' + read)
	return(out)

def pe_fastq_by_lane(lane_sample, read):
	if read == 'Forward':
		out = config['fastq_relative_dir'] + '/' + lane_sample + config['forward_suffix']+ '.fastq.gz'
	elif read == 'Reverse':
		out = config['fastq_relative_dir'] + '/' + lane_sample + config['reverse_suffix'] + '.fastq.gz'
	return(out)
	


# STAR bam
STAR_BAM_OUTPUT = ['STAR_align/' + sample + '/' + 'Aligned.sortedByCoord.out.bam'  \
	for sample in SAMPLE_DICT.keys()] 

# salmon quantification
SALMON_QUANT_OUTPUT = ['salmon_quant/' + sample + '/quant.sf'  \
	for sample in SAMPLE_DICT.keys()] 
	
# fastqc of each sample pair
FASTQC_OUTPUT = ['fastqc/' + sample \
    for sample in SAMPLE_DICT.keys()]

#wildcard_constraints:
#	flowcell_lane_info = "^H\w.*_L\d{3}_"

localrules: all, multiqc


rule all:
	input:
		STAR_BAM_OUTPUT,
		SALMON_QUANT_OUTPUT,
		'fastqc/multiqc_report',
		'salmon_quant/multiqc_report'


rule download_references:
	output:
		#basic_gtf = config['annotation_path'] + 'gencode.v' + config['gencode_version'] + '.basic.annotation.gtf.gz',
		#full_gtf = config['annotation_path'] + 'gencode.v' + config['gencode_version'] + '.annotation.gtf.gz',
		#full_tx = config['annotation_path'] + 'gencode.v' + config['gencode_version'] + '.transcripts.fa.gz',
		pc_tx = config['annotation_path'] + 'gencode.v' + config['gencode_version'] + '.pc_transcripts.fa.gz',
		genome_fa = config['annotation_path'] + config['genome']
	params: 
		fasta_pc_transcripts = 'gencode.v' + config['gencode_version'] + '.pc_transcripts.fa.gz',
		genome_fa = config['genome'],
		ftp = config['ftp'] + config['gencode_version'] + '/'
	shell:
		"""
		rsync -av {params.ftp}{params.fasta_pc_transcripts} {config[annotation_path]}
		rsync -av {params.ftp}{params.genome_fa} {config[annotation_path]}
		"""

rule trim_fastq:
	input:
		r1 = lambda wildcards: pe_fastq_by_lane(wildcards.lane_sample, 'Forward'),
		r2 = lambda wildcards: pe_fastq_by_lane(wildcards.lane_sample, 'Reverse')
	output:
		r1 = 'fastq_trimmed/{lane_sample}' + config['forward_suffix'] + '.fastq.gz',
		r2 = 'fastq_trimmed/{lane_sample}' + config['reverse_suffix'] + '.fastq.gz'
	threads: 8
	conda: 'OGVFB_RNAseq.yml'
	shell:
		"""
		fastp --thread {threads} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2}
		"""

rule fastqc:
	input:
		r1 = lambda wildcards: fastq_by_sample(wildcards.sample, 'Forward'),
		r2 = lambda wildcards: fastq_by_sample(wildcards.sample, 'Reverse'),
	output:
		directory('fastqc/{sample}')
	threads: 4
	conda: 'OGVFB_RNAseq.yml'
	shell:
		"""
		mkdir {output}
		fastqc -f fastq -o {output} -t {threads} {input.r1} {input.r2}
		"""

rule salmon_index:
	input:
		tx = config['annotation_path'] + config['tx'],
		genome = config['annotation_path'] + config['genome']
	output:
		directory(config['annotation_path'] + 'salmon_index_gencode.' + config['gencode_version'] + '.pc_transcripts_salmon130')
	threads: 16
	conda: 'OGVFB_RNAseq.yml'
	shell:
		"""
		grep "^>" <( gunzip -c {input.genome} ) | cut -d " " -f 1 > decoys.txt
		sed -i.bak -e 's/>//g' decoys.txt
		cat {input} > gentrome.fa.gz
		salmon index -t gentrome.fa.gz \
			-d decoys.txt \
			-p {threads} -i {output} -k 31
		"""

rule STAR_index:
	input:
		config['annotation_path'] + config['genome']	
	output:
		directory(config['annotation_path'] + config['genome'] + '_STAR275')
	threads: 16
	conda: 'OGVFB_RNAseq.yml'
	shell:
		"""
		mkdir -p {output}
		zcat {input} > tmp_fasta
		STAR --runThreadN {threads} \
			--runMode genomeGenerate \
			--genomeDir {output} \
			--genomeFastaFiles tmp_fasta \
			--limitGenomeGenerateRAM=140090479317
		rm tmp_fasta 
		"""

rule salmon_quant:
	input:
		r1 = lambda wildcards: fastq_by_sample(wildcards.sample, 'Forward'),
		r2 = lambda wildcards: fastq_by_sample(wildcards.sample, 'Reverse'),
		index = config['annotation_path'] + 'salmon_index_gencode.' + config['gencode_version'] + '.pc_transcripts_salmon130'
	output:
		'salmon_quant/{sample}/quant.sf'
	params:
		'salmon_quant/{sample}'
	threads: 16
	conda: 'OGVFB_RNAseq.yml'
	shell:
		"""
		salmon quant --index {input.index} --libType A --seqBias --gcBias \
			--validateMappings \
			-p {threads} \
			-1 {input.r1} \
			-2 {input.r2} \
			-o {params}
		""" 

rule STAR_align:
	input:
		index = config['annotation_path'] + config['genome'] + '_STAR275',
		r1 = lambda wildcards: fastq_by_sample(wildcards.sample, 'Forward'),
		r2 = lambda wildcards: fastq_by_sample(wildcards.sample, 'Reverse')
	output:
		'STAR_align/{sample}/Aligned.sortedByCoord.out.bam'
	conda: 'OGVFB_RNAseq.yml'
	params:
		scratch = '/lscratch/$SLURM_JOB_ID/STAR_align__{sample}/',
		out = 'STAR_align/{sample}/',
		fq_r1 = lambda wildcards, input: ','.join(input.r1),
		fq_r2 = lambda wildcards, input: ','.join(input.r2)
	threads: 1 
	shell:
		"""
		mkdir -p {params.out}
		STAR --runThreadN {threads} \
			--genomeDir {input.index} \
			--readFilesIn {params.fq_r1} {params.fq_r2} \
			--outSAMtype BAM Unsorted \
			--readFilesCommand zcat \
			--outFileNamePrefix {params.out} \
			--outTmpDir {params.scratch} \
			--limitBAMsortRAM 32000000000 
		mkdir -p {params.scratch}SAMSORT
		samtools sort -T {params.scratch}SAMSORT -o {output} {params.out}Aligned.out.bam 
		rm {params.out}Aligned.out.bam
		samtools index {output}
		"""

rule multiqc:
	input:
		FASTQC_OUTPUT,
		SALMON_QUANT_OUTPUT,
		STAR_BAM_OUTPUT
	output:
		fastqc=directory('fastqc/multiqc_report'),
		salmon=directory('salmon_quant/multiqc_report'),
		star=directory('STAR_align/multiqc_report')
	conda: 'OGVFB_RNAseq.yml'
	shell:
		"""
		multiqc -f -o {output.fastqc} fastqc/
		multiqc -f -o {output.salmon} salmon_quant/
		multiqc -f -o {output.star} STAR_align/
		"""
