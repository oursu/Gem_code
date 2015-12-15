DATADIR=$1

#run cufflinks
#=============
#bam files
G2=${DATADIR}/data/RNAseq/alignments/G2/accepted_hits.bam
G3=${DATADIR}/data/RNAseq/alignments/G3/accepted_hits.bam
V2=${DATADIR}/data/RNAseq/alignments/V2/accepted_hits.bam
V3=${DATADIR}/data/RNAseq/alignments/V3/accepted_hits.bam
#run cufflinks here
label=VEH23,GEM23
OUTDIR=${DATADIR}/results/data_processing/RNAseq/cufflinks/
mkdir -p ${OUTDIR}
maskfile=${DATADIR}/data/RNAseq/referenceData/hg19/maskfiles/rRNA.tRNA.mitoRNA_2013-03-10.gtf
numthreads=3
CUFFDIR=${DATADIR}/bin/cufflinks/cufflinks-2.2.1.Linux_x86_64
gtf=/nfs/vendata/oursu/oana/paper_analysis/mRNAseq/seq_analysis/bowtie_indices/gencodeV15/gencode.v15.annotation.gtf.protein_coding.gtf
CMD="${CUFFDIR}/cuffdiff --labels $label -o $OUTDIR --mask-file $maskfile --raw-mapped-norm --num-threads $numthreads  $gtf $V2,$V3 $G2,$G3" 
echo $CMD
