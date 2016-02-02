RMYPATH=$1
DATADIR=$2
DEFILE=$3
TARGET_NAME_MAPPING=$4
EXPRESSED=$5
PPI=$6

RMYPATH=/nfs/vendata/oursu/oana/GemPaper_2015-12-07/bin/R-3.2.2/bin
DATADIR=/nfs/vendata/oursu/oana/GemPaper_2015-12-07
export R_LIBS=${DATADIR}/bin/R_libraries/
DEFILE=${DATADIR}/data/RNAseq/Table.S3.1.Cuffdiff_protein_coding_V23_vs_G23_gene_exp.diff
TARGET_NAME_MAPPING=${DATADIR}/data/kinase/KinaseNameMapping.txt
EXPRESSED=/nfs/vendata/oursu/oana/GemPaper_2014-12-06/results/DE_analysis/Expressed_above_0.1
PPI=${DATADIR}/data/interactome/9606.mitab.01192011.uniq_miscore-localirefindex3-20110831.digraphno_UBC,EP300_symbol.pkl

#analysis of the kinase screen                                                   
kinase_scripts=${DATADIR}/src/Gem_code/data_processing/kinase/scripts
kinase_results=${DATADIR}/results/data_processing/kinase
${RMYPATH}/Rscript ${DATADIR}/src/Gem_code/general/run_r_rmd.R ${kinase_scripts}/kinase_screenAnalysis.Rmd ${kinase_results}

#analysis of the kinase targets                                         
THRESHOLD=50
FPKMMIN=0.1
intro=${kinase_results}/DATA_kinaseScreen.CELL_PANC1.GENOME_NA.Params_
for kinasehits in ${intro}Hits-Ratio0.3.txt ${intro}Hits-Ratio0.2.txt ${intro}Hits-ttest0.05-Ratio0.3.txt ${intro}Hits-BlissVsLevel0.95.txt;
do                                                                                        
 echo ${kinasehits}=================
 OUT=${kinase_results}/$(basename ${kinasehits} | sed 's/[.]txt//g')
 ${RMYPATH}/Rscript ${kinase_scripts}/kinase_targetAnalysis.R ${DATADIR} ${kinasehits} ${DEFILE} ${THRESHOLD} ${FPKMMIN} ${OUT}/ ${TARGET_NAME_MAPPING}
 KINASE_TEMP=${OUT}/$(basename ${kinasehits} | sed 's/[.]txt//g')".phen-KinaseAffectedDown"
 KINASE_INPUT=${DATADIR}/results/inputs/$(basename ${KINASE_TEMP})-inPPI-expressed.phen
 python ${DATADIR}/src/Gem_code/inputs/prepare_input_for_SAMNet.py --input ${KINASE_TEMP} --input_type phen --PPI ${PPI} --expressed ${EXPRESSED} --out ${KINASE_INPUT}
 KINASE_TEMP=${OUT}/$(basename ${kinasehits} | sed 's/[.]txt//g')".phen-KinaseAffectedUpAndDown"
 KINASE_INPUT=${DATADIR}/results/inputs/$(basename ${KINASE_TEMP})-inPPI-expressed.phen
 python ${DATADIR}/src/Gem_code/inputs/prepare_input_for_SAMNet.py --input ${KINASE_TEMP} --input_type phen --PPI ${PPI} --expressed ${EXPRESSED} --out ${KINASE_INPUT}
done
