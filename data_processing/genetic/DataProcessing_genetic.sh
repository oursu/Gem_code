RPATH=$1
#PARAMS                                                                                                                                                                 
valiHit_cap=$2
pvalue=$3
SI_dev=$4
veh_min_suvival=$5
#paths
genetic_hits_out=$6
DATADIR=$7
PPI=$8
EXPRESSED=$9

#Settings
SRC=${DATADIR}/src/Gem_code
SCREEN_DATA=${DATADIR}/data/genetic/Table.S2.1.GeneticHitsProcessedData.txt
validatedHits=${DATADIR}/data/genetic/validatedHits
EntrezGene2GeneSymbol=${DATADIR}/data/genetic/R_org.Hs.eg.db_2011-08-07_entrez_to_geneSymbols
#Output names
genetic_hits=$(echo ${genetic_hits_out} | sed 's/-withValidated//g')
FINAL=${DATADIR}/results/inputs/$(basename ${genetic_hits_out})-inPPI-expressed$(basename ${EXPRESSED}).phen

#Analyze screen
${RPATH}/Rscript ${SRC}/data_processing/genetic/scripts/geneticScreen_analysis.R ${SCREEN_DATA} ${pvalue} ${SI_dev} ${veh_min_suvival} ${genetic_hits}

#add validated hits
python ${SRC}/data_processing/genetic/scripts/Filter_genetic_hits_for_network_input.py --genetic_hits ${genetic_hits}.tab --genetic_hits_out ${genetic_hits_out} --eg2gs ${EntrezGene2GeneSymbol} --validated ${validatedHits} --scaling ${valiHit_cap}

#prepare for samnet
python ${DATADIR}/src/Gem_code/inputs/prepare_input_for_SAMNet.py --input ${genetic_hits_out}.phen --input_type phen --PPI ${PPI} --expressed ${EXPRESSED} --out ${FINAL}

#And translate the whole screen (for cross-check with the other datasets)
python ${DATADIR}/src/Gem_code/data_processing/genetic/scripts/Translate_geneticData_toSymbol.py --genetic_data ${SCREEN_DATA} --genetic_data_out ${SCREEN_DATA}.translated$(basename ${EntrezGene2GeneSymbol}) --eg2gs ${EntrezGene2GeneSymbol} 