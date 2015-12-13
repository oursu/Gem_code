from optparse import OptionParser
import os
import networkx as nx
import sys
'''
Author:Oana Ursu
'''

def main():
    parser=OptionParser()
    parser.add_option('--input',dest='input',help='Input')
    parser.add_option('--input_type',dest='input_type',help='Input type. Can be phen or DE')
    parser.add_option('--PPI',dest='ppi',help='PPI pkl',
                      default='/nfs/vendata/oursu/oana/Gem_paper/data/interactome/9606.mitab.01192011.uniq_miscore-localirefindex3-20110831.digraphno_UBC,EP300_symbol.pkl')
    parser.add_option('--expressed',dest='expressed',help='Expressed genes',default='/nfs/vendata/oursu/oana/paper_analysis/networks/2013-04-30/_FDR_0.05thresh0.3_vehMin_0.5_lfc_0.585_minFPKM_0.1_p_0.05_pseudocounts_exprProtsFPKM_0.1/0.1FPKM.genes')
    parser.add_option('--TFgene',dest='TFgene',help='TFgene network')
    parser.add_option('--out',dest='out',help='Out name')
    opts,args=parser.parse_args()


    #Get PPI nodes
    if opts.input_type=='phen' or opts.input_type=='TFgene':
        ppi=nx.read_gpickle(opts.ppi)
        ppi_nodes=ppi.nodes()

    #Get expressed
    expressed=set()
    for line in open(opts.expressed,'r').readlines():
        expressed.add(line.strip())

    #setup the analysis for TFgene files here
    if opts.input_type=='TFgene':
        print 'processing tfgene'
        out=open(opts.out,'w')
        counte=0
        for line in open(opts.input,'r').readlines():
            items=line.strip().split()
            tf=items[0]
            gene=items[1]
            score=items[2]
            if tf not in ppi_nodes:
                if tf=='RARB':
                    print 'losing '+tf+' because not in PPI'
                continue
            #if gene not in expressed:
            #    continue
            if tf not in expressed:
                if tf=='RARB':
                    print 'losing '+tf+' because not expressed'
                continue
            out.write(tf+'\t'+gene+'\t'+score+'\n')
        out.close()
        sys.exit()
    #Get TFgene
    TFgene_genes=set()
    if opts.input_type=='DE':
        for line in open(opts.TFgene,'r').readlines():
            TFgene_genes.add(line.strip().split('\t')[1])
    
    #Read in input
    input_genes={}
    for line in open(opts.input,'r').readlines():
        items=line.strip().split('\t')
        if items[0]=='NA':
            continue
        input_genes[items[0]]={}
        input_genes[items[0]]['score']=items[1]

    print 'New dataset ---------------------------------'
    #Check if expressed
    for input_gene in input_genes.keys():
        input_gene_split=input_gene.split('_')
        #print input_gene_split
        for gene in input_gene_split:
            #print 'checking '+gene 
            if gene in expressed:
                print gene+' is expressed'
                input_genes[input_gene]['expressed']=True
            else:
                print gene+' is NOT expressed'
    #For phen, check if in interactome
    if opts.input_type=='phen':
        for input_gene in input_genes.keys():
            genes=input_gene.split('_')
            #Find first gene name that is in interactome and keep it
            for gene in genes:
                if gene in ppi_nodes:
                    if 'inNet' not in input_genes[input_gene].keys():
                        input_genes[input_gene]['inNet']=[]
                    #ONLY ADD IT IF EXPRESSED TOO
                    if gene in expressed:
                        input_genes[input_gene]['inNet'].append(gene)

    #For DE genes, check it is in the TFgene network
    if opts.input_type=='DE':
        for input_gene in input_genes.keys():
            if input_gene in TFgene_genes:
                if 'inNet' not in input_genes[input_gene].keys():
                    input_genes[input_gene]['inNet']=[]
                input_genes[input_gene]['inNet'].append(input_gene)

    #If multiple inputs are in the network under the same name, take the one with the highest score
    gene_to_highest_scoring_input={}
    for input_gene in input_genes.keys():
        if 'inNet' not in input_genes[input_gene].keys():
            continue
        else:
            genes=input_genes[input_gene]['inNet']
            for gene in genes:
                if gene not in gene_to_highest_scoring_input.keys():
                    gene_to_highest_scoring_input[gene]=input_gene
                if float(input_genes[input_gene]['score'])>float(input_genes[gene_to_highest_scoring_input[gene]]['score']):
                        gene_to_highest_scoring_input[gene]=input_gene
    out=open(opts.out,'w')

    #Write down the input to keep
    for input_gene in input_genes.keys():
        #print input_gene
        #print input_genes[input_gene]
        written=False
        if 'expressed' in input_genes[input_gene].keys():
            if input_genes[input_gene]['expressed']==True:
                if 'inNet' in input_genes[input_gene].keys():
                    for gene in input_genes[input_gene]['inNet']:
                        if not written:
                            if gene_to_highest_scoring_input[gene]==input_gene:
                                out.write(gene+'\t'+input_genes[input_gene]['score']+'\n')
                                written=True
    out.close()
    




  

if __name__=='__main__':
    main()
