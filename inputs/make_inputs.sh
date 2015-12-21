
#General
DATADIR=/nfs/vendata/oursu/oana/GemPaper_2015-12-07
RMYPATH=${DATADIR}/bin/R-3.2.2/bin
export R_LIBS=$R_LIBS:${DATADIR}/bin/R_libraries/

#Inputs to network
PPI=${DATADIR}/data/interactome/9606.mitab.01192011.uniq_miscore-localirefindex3-20110831.digraphno_UBC,EP300_symbol.pkl
EXPRESSED=/nfs/vendata/oursu/oana/GemPaper_2014-12-06/results/DE_analysis/Expressed_above_0.1 #### change this                                                            
#Data processing - genetic screen
#================================
valiHit_cap=3
pvalue=0.05
SI_dev=0.3
veh_min_suvival=0.5
GENETIC_RESULTS=${DATADIR}/results/data_processing/genetic/
GENETIC_OUT=${GENETIC_RESULTS}/DATA_geneticScreen.CELL_PANC1.TREAT_screenForGemResistance.GENOME_hg19.PARAMS_p${pvalue}-vehSurvival${veh_min_suvival}-SIdev${SI_dev}-withValidated
${DATADIR}/src/Gem_code/data_processing/genetic/DataProcessing_genetic.sh ${RMYPATH} ${valiHit_cap} ${pvalue} ${SI_dev} ${veh_min_suvival} ${GENETIC_OUT} ${DATADIR} ${PPI} ${EXPRESSED}

#Data processing - RNAseq 
#========================
#alignments (tophat2)
#${DATADIR}/src/Gem_code/data_processing/RNAseq/align_RNAseq.sh
${DATADIR}/src/Gem_code/data_processing/RNAseq/DataProcessing_RNAseq.sh 

# Kinase analysis and input creation
#=====================================
DEFILE=${DATADIR}/data/RNAseq/Table.S3.1.Cuffdiff_protein_coding_V23_vs_G23_gene_exp.diff
${DATADIR}/src/Gem_code/data_processing/kinase/DataProcessing_kinase.sh ${RMYPATH} ${DATADIR} ${DEFILE} ${DATADIR}/data/kinase/KinaseNameMapping.txt ${EXPRESSED} ${PPI}

#DNase data analysis
#=====================
#get encode data (peaks)
cd ${DATADIR}/data/TFgene/DNASEseq/ENCODE
#wget https://www.encodeproject.org/files/ENCFF001EEF/@@download/ENCFF001EEF.fastq.gz
#wget https://www.encodeproject.org/files/ENCFF001EEG/@@download/ENCFF001EEG.fastq.gz
#Peaks for rep1 and rep2 from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM736519 in ${DATADIR}/data/TFgene/DNASEseq/ENCODE/peaks/



