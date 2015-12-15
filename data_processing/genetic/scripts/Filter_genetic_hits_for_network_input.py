from optparse import OptionParser
import os
import re

'''
Author:Oana Ursu
'''

def main():
    parser=OptionParser()
    parser.add_option('--genetic_hits',dest='genetic_hits',help='Genetic hits')
    parser.add_option('--genetic_hits_out',dest='genetic_hits_out',help='Genetic hits out')
    parser.add_option('--eg2gs',dest='eg2gs',help='Entrez to gene symbols')
    parser.add_option('--validated',dest='vali_hits',help='Experimentally validated genetic hits')
    parser.add_option('--scaling',dest='scaling',help='Capacity for validated hits')
    opts,args=parser.parse_args()

    #SET IN STONE INPUTS
    #entrez_to_geneSymbol_file=open('/nfs/vendata/oursu/oana/gene_names/data/R_org.Hs.eg.db_2011-08-07_entrez_to_geneSymbols','r')
    #validated_hits='/nfs/vendata/oursu/oana/ResponseNet/datasets/2012-08-20/top_27_'
    #scaling='3'
    entrez_to_geneSymbol_file=open(opts.eg2gs,'r')
    validated_hits=opts.vali_hits
    scaling=opts.scaling
    intermediary_out=re.sub('.tab','',opts.genetic_hits)+'-withoutValidatedHits.phen'

    out=open(intermediary_out,'w')

    #1. Naming
    #=========
    #create a dictioary to use for screen                                                                                                                             
    entrez_geneSymbol_dict={}
    for line in entrez_to_geneSymbol_file.readlines():
        fields=line.strip('\n').split('\t')
        entrez_geneSymbol_dict[str(fields[0])]=fields[1]
    #Read in the genetic hits, and write them down as the new gene symbols
    screenfile=open(opts.genetic_hits,'r')

    screenlines=screenfile.readlines()
    screenlines=screenlines[1:len(screenlines)]
    for line in screenlines:
        fields=line.strip().split('\t')
        #print fields
        entrez=fields[5]
        if entrez=='NA':
            continue
        gene_symbol=fields[4]
        if entrez in entrez_geneSymbol_dict.keys():
            gene_symbol=entrez_geneSymbol_dict[entrez]
        out.write(gene_symbol+'\t'+str(abs(float(fields[20])))+'\n')
    out.close()

    #2. Augment genetic hits with validated dataset                                                                                                                       
    #==============================================                                                                                                                       
    out=open(opts.genetic_hits_out+'.phen','w')
    top_hits=set()                               
    for line in open(validated_hits,'r').readlines():                                                                                                                     
        hit=line.strip()
        #print hit
        top_hits.add(hit) 
    for line in open(intermediary_out,'r').readlines():
        items=line.strip().split('\t')                                                                                                                             
        hit=items[0]                                                                    
        if hit not in top_hits:                                                                                                                                           
            out.write(line)                                                                                                                                               
    for hit in top_hits:                                                                                                                                                 
        out.write(hit+'\t'+scaling+'\n')                                                                                                                                  
    out.close()                                                                                                                                                           
                                                                                                                
    #3. Keep only expressed genes in the interactome - to come. Was not doing this before, so not doing it now, for now.
    
    #4. Print out a stats file



  

if __name__=='__main__':
    main()
