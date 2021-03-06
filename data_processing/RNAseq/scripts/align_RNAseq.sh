
DATADIR=/nfs/vendata/oursu/oana/GemPaper_2015-12-07
TH=/nfs/vendata/oursu/oana/GemPaper_2015-12-07/bin/tophat/tophat-2.1.0.Linux_x86_64/tophat2 

v2=/nfs/antdata/analysis/110929_F/barcodes/110929_V2_ACACCT.._yoonsing_Fraenkel_L8_1_sequence.txt
v3=/nfs/antdata/analysis/110929_F/barcodes/110929_V3_CTCTCT.._yoonsing_Fraenkel_L8_1_sequence.txt
g2=/nfs/antdata/analysis/110929_F/barcodes/110929_G2_CTCTTG.._yoonsing_Fraenkel_L8_1_sequence.txt
g3=/nfs/antdata/analysis/110929_F/barcodes/110929_G3_GAGATG.._yoonsing_Fraenkel_L8_1_sequence.txt
gtf=${DATADIR}/data/RNAseq/referenceData/gencode.v19.annotation.protein_coding.gtf
bw2_index=${DATADIR}/data/RNAseq/referenceData/hg19/hg19

for sp in ${v2} ${v3} ${g2} ${g3};
do
 sample=$(echo $(basename ${sp}) | sed 's/_/\t/g' | cut -f2)
 echo ${sample}
 OUT=${DATADIR}/data/RNAseq/alignments/${sample}
 s=${OUT}/align_${sample}_script.sh
 mkdir -p ${OUT}
 echo "export PATH="'$'"PATH:/nfs/vendata/oursu/oana/GemPaper_2015-12-07/bin/bowtie2/bowtie2-2.2.6/" > ${s}
 echo "${TH} -G ${gtf} --no-novel-juncs --max-insertion-length 0 --max-deletion-length 0 --library-type fr-unstranded  --solexa1.3-quals --keep-tmp -o ${OUT} ${bw2_index} ${sp}" >> ${s}
 chmod 755 ${s}
done

