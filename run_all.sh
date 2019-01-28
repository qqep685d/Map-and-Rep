#!/bin/bash
#=====================================
#===== Workflow of this pipline ======
#=====================================
# STEP1: BLAST / Filtering non-chloroplast contigs
# STEP2: Cutting non-homologous region in each contigs
# STEP3: Extracting no-overlapped contigs
# STEP4: Mapping no-overlapped contigs to reference sequence
#=====================================



#=====================================
#============= Settings ==============
#=====================================
#=== INPUT ===
## PATH to each files
REF_FASTA="" # fasta file of Complete chloroplast genome (single fasta)
IN_FASTA="" # fasta file of Assembly contigs
OUT_NAME="MySample"     # using this as the prefix of output files

#=== OUTPUT ===
DIR01="01_blast_output" # output directory of STEP1
DIR02="02_homolog_output" # output directory of STEP2
DIR03="03_no_overlap_output" # output directory of STEP3
# DIR04="04_chromosomer_output" # output directory of STEP4
DIR04="04_repseq_output" # output directory of STEP4

#=== BLAST settings of STEP1 ===
threshold_evalue="1e-10"     # Threshold of e-value
threshold_perc_identity="80" # Threshold of Percent of Identity

#=== settings of STEP3 ===
threshold_aln_length="0"
## If alignment length is less than given values, the contig is not used.

#=====================================

BLS_RESULT="${DIR01}/${OUT_NAME}.blast.txt"
BLS_FASTA="${DIR01}/${OUT_NAME}.blast.fasta"
HLG_RESULT="${DIR02}/${OUT_NAME}.homolog.txt"
HLG_FASTA="${DIR02}/${OUT_NAME}.homolog.fasta"
NOL_RESULT="${DIR03}/${OUT_NAME}.no_overlap.txt"
NOL_FASTA="${DIR03}/${OUT_NAME}.no_overlap.fasta"
SCF_NAME="${DIR04}/${OUT_NAME}"

#=====================================

# #=== [STEP1] ==
function run_step1(){
	## Make BLAST Database
	makeblastdb -in ${REF_FASTA} -out ${REF_FASTA} -dbtype nucl -parse_seqids

	# ## Run BLAST
	mkdir -p ${DIR01}
	blastn -num_threads 1 -db ${REF_FASTA} -query ${IN_FASTA} -evalue ${threshold_evalue} -perc_identity ${threshold_perc_identity} -outfmt "6" > ${BLS_RESULT}

	python program/extract_seqs.py --fasta ${IN_FASTA} --list ${BLS_RESULT} > ${BLS_FASTA}
}

#=== [STEP2] ===
function run_step2(){
	mkdir -p ${DIR02}

	## homologous position list
	python program/cutting_position.py --list ${BLS_RESULT}
	cat 'step2_temp' > ${HLG_RESULT}
	rm -f 'step2_temp'

	## extract homologous sequences in each contigs
	python program/cut_seqs.py --fasta ${BLS_FASTA} --list ${HLG_RESULT} > ${HLG_FASTA}
}

#=== [STEP3] ===
function run_step3(){
	mkdir -p ${DIR03}

	## list of no overlapped contigs
	python program/no_overlap_contigs.py --list ${HLG_RESULT} --cutoff ${threshold_aln_length}
	cat 'step3_temp' > ${NOL_RESULT}
	rm -f 'step3_temp'

	## extract contigs
	python program/extract_seqs.py --fasta ${HLG_FASTA} --list ${NOL_RESULT} > ${NOL_FASTA}
}

#=== [STEP4] ===
# ## Use e-RGA
# function run_step4_e_rga(){
# 	mkdir -p ${DIR04}
#
# 	perl program/e-RGA.pl -ref ${REF_FASTA} -de-novo ${HLG_FASTA} -out 'step3_temp'
#
# 	cat 'step3_temp' > ${SCF_FASTA}
# 	cat 'mergeListstep3_temp' > ${SCF_FASTA}.mergeList
# 	rm -f 'step3_temp' 'mergeListstep3_temp'
# }

# ## Use Chromosomer
# function run_step4_chromosomer(){
# 	mkdir -p ${DIR04}
#
# 	#makeblastdb -in ${REF_FASTA} -out ${REF_FASTA} -dbtype nucl -parse_seqids
#
# 	blastn -num_threads 1 -db ${REF_FASTA} -query ${NOL_FASTA} -evalue ${threshold_evalue} -perc_identity ${threshold_perc_identity} -outfmt "6" > "${SCF_NAME}.blastn.txt"
#
# 	chromosomer fastalength ${NOL_FASTA} "${SCF_NAME}.fasta_length.txt"
#
# 	chromosomer fragmentmap "${SCF_NAME}.blastn.txt" 500 "${SCF_NAME}.fasta_length.txt" "${SCF_NAME}.fragment_map.txt"
#
# 	chromosomer fragmentmapstat "${SCF_NAME}.fragment_map.txt" "${SCF_NAME}.fragment_map_stat.txt"
#
# 	chromosomer fragmentmapbed "${SCF_NAME}.fragment_map.txt" "${SCF_NAME}.fragment_map.bed"
#
# 	chromosomer assemble "${SCF_NAME}.fragment_map.txt" ${NOL_FASTA} "${SCF_NAME}.chromosomer.fasta"
# }

## Use in-house programs
function run_step4_make_seqs(){
	mkdir -p ${DIR04}

	python program/replace_map.py --ref ${REF_FASTA} --list ${NOL_RESULT}
	cat 'step4_temp' > "${SCF_NAME}.map.txt"
	rm -f 'step4_temp'

	python program/replace_seqs.py --ref ${REF_FASTA} --contig ${NOL_FASTA} --list "${SCF_NAME}.map.txt"
	cat 'step4_tempN' > "${SCF_NAME}.gap_filled_NNN.fasta"
	cat 'step4_tempR' > "${SCF_NAME}.gap_filled_Ref.fasta"
	rm -f 'step4_tempN' 'step4_tempR'
}



#=====================================
#=============== Main ================
#=====================================
# run_step1
# run_step2
# run_step3
run_step4_make_seqs

#=====================================
#======= Related information =========
#=====================================
#=== Softwares ===
## STEP1/STEP2
##
## Software: BLAST version 2.5.0
## Homepage: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
##
## Python scripts (in-house scripts):
## - cut_seqs.py : cut sequences in each contigs by using list of cutting position
## - cutting_position.py : make list about cutting position
## - extract_seqs.py : extract contigs in list
## - get_seq_length.py : get length of sequences
## - no_overlap_contigs.py : extract contigs with no-overlapped region
## - replace_map.py : make list about cutting position
## - replace_seqs.py : extract contigs with sequences replaced by reference nucleotides or 'N's
