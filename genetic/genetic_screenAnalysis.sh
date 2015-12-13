#Settings
DATADIR=/nfs/vendata/oursu/oana/GemPaper_2015-12-07
SCREEN_ANALYSIS=${DATADIR}/results/data_processing/genetic/
validatedHits=${DATADIR}/data/genetic/validatedHits
valiHit_cap=3
EntrezGene2GeneSymbol=${DATADIR}/data/genetic/R_org.Hs.eg.db_2011-08-07_entrez_to_geneSymbols
pvalue=0.05
SI_dev=0.3
veh_min_suvival=0.5
genetic_hits=${SCREEN_ANALYSIS}/GeneticHits_Pval${pvalue}_vehSurvival${veh_min_suvival}_SIdev${SI_dev}.tab
genetic_hits_out=${genetic_hits}_with_validated.phen
README=${SCREEN_ANALYSIS}/README


#Analyze screen
mkdir ${SCREEN_ANALYSIS}
Rscript ${DATADIR}/src/genetic/Genetic_hits_analysis.R ${DATADIR} ${DATADIR}/results/data_processing/genetic/ ${pvalue} ${SI_dev} ${veh_min_suvival} ${genetic_hits}

#Prepare genetic hits for SAMNet
python ${DATADIR}/src/genetic/Filter_genetic_hits_for_network_input.py --genetic_hits ${genetic_hits} --genetic_hits_out ${genetic_hits_out} --eg2gs ${EntrezGene2GeneSymbol} --validated ${validatedHits} --scaling ${valiHit_cap}

echo '==========================' > ${README}
echo I ran ${DATADIR}/src/genetic/genetic_screenAnalysis.sh >> ${README}
echo I have created SAMNet INPUT.Genetic as ${genetic_hits_out} >> ${README}
echo '==========================' >> ${README}

