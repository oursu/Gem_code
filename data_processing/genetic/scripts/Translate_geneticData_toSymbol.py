from optparse import OptionParser
import os

'''
Author:Oana Ursu
'''

def main():
    parser=OptionParser()
    parser.add_option('--genetic_data',dest='genetic_data',help='Genetic data to be translated')
    parser.add_option('--genetic_data_out',dest='genetic_data_out',help='Genetic data out')
    parser.add_option('--eg2gs',dest='eg2gs',help='Entrez to gene symbols')
    opts,args=parser.parse_args()

    entrez_to_geneSymbol_file=open(opts.eg2gs,'r')
    out=open(opts.genetic_data_out,'w')

    #1. Naming
    #=========
    #create a dictioary to use for screen                                                                                                                             
    entrez_geneSymbol_dict={}
    print 'making dict'
    counter=1
    for line in entrez_to_geneSymbol_file.readlines():
        counter=counter+1
        print counter
        fields=line.strip('\n').split('\t')
        entrez_geneSymbol_dict[str(fields[0])]=fields[1]
    #Read in the genetic data, and write them down as the new gene symbols
    screenfile=open(opts.genetic_data,'r')

    screenlines_orig=screenfile.readlines()
    out.write('TraslatedGeneSymbol'+'\t'+screenlines_orig[0])
    screenlines=screenlines_orig[1:len(screenlines_orig)]
    print 'writing data'
    counter=1
    for line in screenlines:
        counter=counter+1
        print counter
        fields=line.strip().split('\t')
        #print fields
        entrez=fields[5]
        gene_symbol=fields[4]+'original'
        if entrez in entrez_geneSymbol_dict.keys():
            gene_symbol=entrez_geneSymbol_dict[entrez]
        out.write(gene_symbol+'\t'+line)
    out.close()
  

if __name__=='__main__':
    main()
