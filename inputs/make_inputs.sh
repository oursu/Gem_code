
#General
RMYPATH=/nfs/vendata/oursu/oana/GemPaper_2015-12-07/src/R/R-3.2.2/bin
export R_LIBS=$R_LIBS:/nfs/vendata/oursu/oana/GemPaper_2015-12-07/src/R/my_libraries/
DATADIR=/nfs/vendata/oursu/oana/GemPaper_2015-12-07
#Params
KINASE_ACTIVITYREDUCTION_THRESHOLD=50
#Inputs to network
PPI=${DATADIR}/data/interactome/9606.mitab.01192011.uniq_miscore-localirefindex3-20110831.digraphno_UBC,EP300_symbol.pkl
KINASE_INPUT=${DATADIR}/results/inputs/KINASE_INPUT_${KINASE_ACTIVITYREDUCTION_THRESHOLD}.phen
EXPRESSED=/nfs/vendata/oursu/oana/GemPaper_2014-12-06/results/DE_analysis/Expressed_above_0.1 #### change this

# Kinase analysis and input creation
#=====================================
#analysis of the kinase screen 
${RMYPATH}/Rscript ${DATADIR}/src/kinase/RUN_kinase_screenAnalysis.R ${DATADIR}
#analysis of the kinase targets
${RMYPATH}/Rscript ${DATADIR}/src/kinase/RUN_kinase_targetAnalysis.R ${DATADIR}
#converting the kinase results to SAMNet inputs
KINASE_TEMP=${DATADIR}/results/data_processing/kinase/KinaseTargets_for_phen_thresholdActivityReduction${KINASE_ACTIVITYREDUCTION_THRESHOLD}.phen
python ${DATADIR}/src/inputs/prepare_input_for_SAMNet.py --input ${KINASE_TEMP} --input_type phen --PPI ${PPI} --expressed ${EXPRESSED} --out ${KINASE_INPUT}

#Genetic input creation
#======================