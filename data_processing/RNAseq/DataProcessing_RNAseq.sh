
SIG=0.05
FPKMmin=0.1
RMYPATH=/nfs/vendata/oursu/oana/GemPaper_2015-12-07/bin/R-3.2.2/bin
DATADIR=/nfs/vendata/oursu/oana/GemPaper_2015-12-07
export R_LIBS=${DATADIR}/bin/R_libraries/
DEFILE=${DATADIR}/data/RNAseq/Table.S3.1.Cuffdiff_protein_coding_V23_vs_G23_gene_exp.diff
PPI=${DATADIR}/data/interactome/9606.mitab.01192011.uniq_miscore-localirefindex3-20110831.digraphno_UBC,EP300_symbol.pkl
EXPRESSED=${DATADIR}/results/data_processing/RNAseq/cufflinks_DE/Expressed${FPKMmin}

#Assumes quantification with cufflinks, file=${DATADIR}/data/RNAseq/Table.S3.1.Cuffdiff_protein_coding_V23_vs_G23_gene_exp.diff
#get DE genes
OUT=${DATADIR}//results/data_processing/RNAseq/cufflinks_DE
mkdir -p ${OUT}
intro=DATA_RNAseq.CELL_PANC1.GENOME_hg19.TREAT_wtVehvsGem.PARAMS_Cufflinks-proteincoding-V23vsG23-SIG${SIG}-FPKMmin${FPKMmin}
${RMYPATH}/Rscript ${DATADIR}/src/Gem_code/data_processing/RNAseq/scripts/DEgenes_fromCufflinks.R ${DEFILE} ${SIG} ${FPKMmin} ${OUT}/${intro} ${EXPRESSED}

#plot GO enrichment for these genes
OUT=${DATADIR}/results/data_processing/RNAseq/cufflinks_DE
for gofile in $(ls ${OUT}/*GO.txt);
do
 echo $gofile
 ${RMYPATH}/Rscript ${DATADIR}/src/Gem_code/data_processing/RNAseq/scripts/DEgenes_plotGO.R ${gofile} ${gofile}.pdf ${SIG}
done

#Venn diagram
GENETIC=${DATADIR}/results/data_processing/genetic/DATA_geneticScreen.CELL_PANC1.TREAT_screenForGemResistance.GENOME_hg19.PARAMS_p0.05-vehSurvival0.5-SIdev0.3-withValidated.phen
TFFILE=${DATADIR}/results/inputs/DATA_ENCODEplusOurDNase.CELL_PANC1.TREAT_wtVehAndGem.GENOME_hg19.PARAMS_fromINPUT4TFgene.TFA.tfa
K=${DATADIR}/results/data_processing/kinase/DATA_kinaseScreen.CELL_PANC1.GENOME_NA.Params_
DE_GENES=${DATADIR}//results/data_processing/RNAseq/cufflinks_DE/DATA_RNAseq.CELL_PANC1.GENOME_hg19.TREAT_wtVehvsGem.PARAMS_Cufflinks-proteincoding-V23vsG23-SIG${SIG}-FPKMmin${FPKMmin}.UpAndDownReg.txt
for kinaseparams in Hits-Ratio0.3.txt Hits-Ratio0.2.txt Hits-ttest0.05-Ratio0.3.txt Hits-BlissVsLevel0.95.txt;
do 
 kinase=${K}${kinaseparams}
 for direction in Down UpAndDown;
 do
  KINASE_EXP=$(echo ${kinase} | sed 's/[.]txt//g')/$(basename ${kinase} | sed 's/[.]txt//g')".phen-KinaseAffected"${direction}"-expressed"
  ls -lh   ${KINASE_EXP}
  intro=DATA_RNAseq.CELL_PANC1.GENOME_hg19.TREAT_wtVehvsGem.Params_Cufflinks-proteincoding-V23vsG23-SIG${SIG}-FPKMmin${FPKMmin}-withKinase-${kinaseparams}${direction}
  code=${DATADIR}/src/Gem_code/data_processing/RNAseq/scripts/DEgenes_VennDiagram.R
  ${RMYPATH}/Rscript ${code} ${DATADIR} ${OUT}/${intro} ${DE_GENES} ${GENETIC} ${KINASE_EXP} ${DEFILE} ${FPKMmin} > ${OUT}/${intro}.intersectionsText ${TFFILE}
 done
done

#tf enrichments
intro=DATA_RNAseq.CELL_PANC1.GENOME_hg19.TREAT_wtVehvsGem.Params_Cufflinks-proteincoding-V23vsG23-SIG${SIG}-FPKMmin${FPKMmin}
${RMYPATH}/Rscript ${DATADIR}/src/Gem_code/data_processing/RNAseq/scripts/DEgenes_TFenrichment.R ${DE_GENES} ${DEFILE} ${TFFILE} ${OUT}/${intro}.TFenrichments.pdf ${SIG}


#prepare for samnet
FINAL=${DATADIR}/results/inputs/${intro}-inPPI-expressed.phen
python ${DATADIR}/src/Gem_code/inputs/prepare_input_for_SAMNet.py --input ${DE_GENES} --input_type DE --PPI ${PPI} --expressed ${EXPRESSED} --out ${FINAL} --TFgene ${TFFILE}