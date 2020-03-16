# pull file names and paths
SAMPLE_FASTQ_TABLE = open(config['sample_fastq_table']) # tab separated text file (no column names) with sample name and fastq file 
SAMPLE = []
FASTQ = []
for line in SAMPLE_FASTQ_TABLE:
	SAMPLE.append(line.split()[0])
	FASTQ.append(line.split()[1])
# salmon quantification
SALMON_QUANT_OUTPUT = ['salmon_quant/' + f.split('_R1')[0] + '/quant.sf'  \
	for f in FASTQ if 'R1' in f] 
	
# fastqc of each fastq file pair
FASTQC_OUTPUT = ['fastqc/' + f.split('_R1')[0] \
    for f in FASTQ if 'R1' in f]

#wildcard_constraints:
#	flowcell_lane_info = "^H\w.*_L\d{3}_"

localrules: all, download_references


rule all:
	input:
		SALMON_QUANT_OUTPUT,
		'fastqc/multiqc_report',
		'salmon_quant/multiqc_report'


rule download_references:
	output:
		basic_gtf = config['annotation_path'] + 'gencode.v33.basic.annotation.gtf.gz',
		full_gtf = config['annotation_path'] + 'gencode.v33.annotation.gtf.gz',
		full_tx = config['annotation_path'] + 'gencode.v33.transcripts.fa.gz',
		pc_tx = config['annotation_path'] + 'gencode.v33.pc_transcripts.fa.gz',
		genome_fa = config['annotation_path'] + 'GRCh38.p13.genome.fa.gz'
	params: 
		gtf_basic = 'gencode.v33.basic.annotation.gtf.gz',
		gtf = 'gencode.v33.annotation.gtf.gz',
		fasta_transcripts = 'gencode.v33.transcripts.fa.gz',
		fasta_pc_transcripts = 'gencode.v33.pc_transcripts.fa.gz',
		genome_fa = 'GRCh38.p13.genome.fa.gz',
		ftp = 'rsync://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/'
	shell:
		"""
		rsync -av {params.ftp}{params.gtf_basic} {config[annotation_path]}
		rsync -av {params.ftp}{params.gtf} {config[annotation_path]}
		rsync -av {params.ftp}{params.fasta_transcripts} {config[annotation_path]}
		rsync -av {params.ftp}{params.fasta_pc_transcripts} {config[annotation_path]}
		rsync -av {params.ftp}{params.genome_fa} {config[annotation_path]}
		"""

rule fastqc:
	input:
		r1 = 'fastq_files/{flowcell_lane_info}_R1_001.fastq.gz',
		r2 = 'fastq_files/{flowcell_lane_info}_R2_001.fastq.gz'
	output:
		directory('fastqc/{flowcell_lane_info}')
	threads: 4
	conda: 'OGVFB_RNAseq.yml'
	shell:
		"""
		mkdir {output}
		fastqc -f fastq -o {output} -t {threads} {input.r1} {input.r2}
		"""

rule salmon_index:
	input:
		tx = config['annotation_path'] + 'gencode.v33.pc_transcripts.fa.gz',
		genome = config['annotation_path'] + 'GRCh38.p13.genome.fa.gz'
	output:
		directory(config['annotation_path'] + 'salmon_index_gencode.v33.pc_transcripts_salmon110')
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

rule salmon_quant:
	input:
		r1 = 'fastq_files/{flowcell_lane_info}_R1_001.fastq.gz',
		r2 = 'fastq_files/{flowcell_lane_info}_R2_001.fastq.gz',
		index = config['annotation_path'] + 'salmon_index_gencode.v33.pc_transcripts_salmon110'
	output:
		'salmon_quant/{flowcell_lane_info}/quant.sf'
	params:
		'salmon_quant/{flowcell_lane_info}'
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

rule multiqc:
	input:
		FASTQC_OUTPUT,
		SALMON_QUANT_OUTPUT
	output:
		fastqc=directory('fastqc/multiqc_report'),
		salmon=directory('salmon_quant/multiqc_report')
	conda: 'OGVFB_RNAseq.yml'
	shell:
		"""
		multiqc -f -o {output.fastqc} fastqc/
		multiqc -f -o {output.salmon} salmon_quant/
		"""
