
DATADIR=/nfs/vendata/oursu/oana/GemPaper_2015-12-07

RSEM=${DATADIR}/bin/rsem/RSEM-1.2.25/rsem
index=${DATADIR}/data/RNAseq/referenceData/encodeHg19Male/rsem_index/encodeHg19Male.hg19.gencode.v19.annotation.protein_coding

for data in G3;
do
 echo ${data}
 bam=${DATADIR}/data/RNAseq/alignments/${data}_ucscgenome/accepted_hits.bam
 ls -lh  ${bam}
 OUT=${DATADIR}/results/data_processing/RNAseq/RSEM/${data}
 mkdir -p ${OUT}
 ${RSEM}-calculate-expression --bam ${bam} ${index} ${OUT}/${data}
 #fastq=/nfs/antdata/analysis/110929_F/barcodes/110929_V2_ACACCT.._yoonsing_Fraenkel_L8_1_sequence.txt
 #${RSEM}-calculate-expression --bowtie2 --bowtie2-path ${DATADIR}/bin/bowtie2/bowtie2-2.2.6/ ${fastq} ${index} ${OUT}/${data}
done

#would be really cool because out it would give me the values to put into DEseq



RSEM=${DATADIR}/bin/rsem/RSEM-1.2.25/rsem
index=${DATADIR}/data/RNAseq/referenceData/encodeHg19Male/rsem_index/encodeHg19Male.hg19.gencode.v19.annotation.protein_coding
OUT=${DATADIR}/results/data_processing/RNAseq/RSEM/${data}
${RSEM}-calculate-expression -p 8 --bam --estimate-rspd --append-names --output-genome-bam ${bam} ${index}
${RSEM}-calculate-expression --bam ${bam} ${index} ${data}