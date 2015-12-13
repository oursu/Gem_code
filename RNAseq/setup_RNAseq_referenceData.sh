
DATADIR=/nfs/vendata/oursu/oana/GemPaper_2015-12-07

gencode_gtf=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
cd ${DATADIR}/data/RNAseq/referenceData/
wget ${gencode_gtf}
#keep only protein coding
zcat -f gencode.v19.annotation.gtf.gz | grep protein_coding | gzip > gencode.v19.annotation.protein_coding.gtf.gz
gunzip gencode.v19.annotation.protein_coding.gtf.gz

#get encode genome
mkdir -p ${DATADIR}/data/RNAseq/referenceData/encodeHg19Male
cd ${DATADIR}/data/RNAseq/referenceData/encodeHg19Male
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.chrom.sizes
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz
for chromo in {1..22};
do
 wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/femaleByChrom/chr${chromo}.fa.gz
done
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/femaleByChrom/chrX.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/femaleByChrom/chrM.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/maleByChrom/chrY.fa.gz

#make bowtie2 index
mkdir -p ${DATADIR}/data/RNAseq/referenceData/encodeHg19Male/bowtie2_index
FASTA=${DATADIR}/data/RNAseq/referenceData/encodeHg19Male/male.hg19.fa
zcat ${FASTA}.gz > ${FASTA}
bowtiepath=/nfs/vendata/oursu/oana/GemPaper_2015-12-07/bin/bowtie2/bowtie2-2.2.6/bowtie2
${bowtiepath}-build ${FASTA} ${DATADIR}/data/RNAseq/referenceData/encodeHg19Male/bowtie2_index/encodeHg19Male.hg19

#make rsem reference files
rsempath=${DATADIR}/bin/rsem/RSEM-1.2.25/rsem
pcgtf=${DATADIR}/data/RNAseq/referenceData/gencode.v19.annotation.protein_coding.gtf
mkdir -p ${DATADIR}/data/RNAseq/referenceData/encodeHg19Male/rsem_index/
${rsempath}-prepare-reference --bowtie2 --bowtie2-path /nfs/vendata/oursu/oana/GemPaper_2015-12-07/bin/bowtie2/bowtie2-2.2.6/ --gtf ${pcgtf} ${DATADIR}/data/RNAseq/referenceData/encodeHg19Male/male.hg19.fa ${DATADIR}/data/RNAseq/referenceData/encodeHg19Male/rsem_index/encodeHg19Male.hg19.gencode.v19.annotation.protein_coding

