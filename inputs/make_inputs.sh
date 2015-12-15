
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
${DATADIR}/src/Gem_code/data_processing/RNAseq/DataProcessing_RNAseq.sh ${DATADIR} 


# Kinase analysis and input creation
#=====================================
KINASE_INPUT=${DATADIR}/results/inputs/KINASE_INPUT_${KINASE_ACTIVITYREDUCTION_THRESHOLD}.phen
KINASE_ACTIVITYREDUCTION_THRESHOLD=50

#analysis of the kinase screen 
${RMYPATH}/Rscript ${DATADIR}/src/kinase/RUN_kinase_screenAnalysis.R ${DATADIR}
#analysis of the kinase targets
${RMYPATH}/Rscript ${DATADIR}/src/kinase/RUN_kinase_targetAnalysis.R ${DATADIR}
#converting the kinase results to SAMNet inputs
KINASE_TEMP=${DATADIR}/results/data_processing/kinase/KinaseTargets_for_phen_thresholdActivityReduction${KINASE_ACTIVITYREDUCTION_THRESHOLD}.phen
python ${DATADIR}/src/inputs/prepare_input_for_SAMNet.py --input ${KINASE_TEMP} --input_type phen --PPI ${PPI} --expressed ${EXPRESSED} --out ${KINASE_INPUT}

#DNase data analysis
#=====================
#get encode data (peaks)
cd ${DATADIR}/data/TFgene/DNASEseq/ENCODE
#wget https://www.encodeproject.org/files/ENCFF001EEF/@@download/ENCFF001EEF.fastq.gz
#wget https://www.encodeproject.org/files/ENCFF001EEG/@@download/ENCFF001EEG.fastq.gz
#Peaks for rep1 and rep2 from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM736519 in ${DATADIR}/data/TFgene/DNASEseq/ENCODE/peaks/



