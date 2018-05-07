#!/usr/bin/python

import sys
import os
import math
import argparse

import Bio.Ontology
import Bio.Ontology.IO as OntoIO



def read_list(filename):
    
    out = []
    with open(filename, 'r') as file_in:
        line = file_in.readline()  #header
        if not (line[0] == '!' or line[0] == '#'): file_in.seek(0)

        for line in file_in:
            content = line.strip().split('\t') 
            if len(content) <= 1:
                if content[0] != "":
                    if len(content[0].split('_')) < 2:
                        out.append( content[0])
                    else:
                        out.append( "_".join(content[0].split('_')[1:-1]))
            elif content[1] == '1':
                if len(content[0].split('_')) < 2:
                    out.append( content[0])
                else:
                    out.append( "_".join(content[0].split('_')[1:-1]))
            elif content[1] != '0':
                raise Exception("Invalid values in list of genes: second column includes %s instead of 0 or 1"%content[1])
    return out




def run_term(assocs, go_graph, gene_list, corrections):
    from Bio.Ontology import TermForTermEnrichmentFinder

    ef = TermForTermEnrichmentFinder(assocs, go_graph)
    result = ef.find_enrichment(gene_list, corrections)
    return result

def run_parent_child(assocs, go_graph, gene_list, corrections, method):
    from Bio.Ontology import ParentChildEnrichmentFinder

    ef = ParentChildEnrichmentFinder(assocs, go_graph)
    result = ef.find_enrichment(gene_list, corrections, method)
    return result




def check_file(parser, arg, openparam):
    if openparam == 'r':
        if not os.path.exists(arg):
            parser.error("The file %s does not exist!" % arg)
    else:
        try:
            f=open(arg, openparam)
            f.close()
        except:
            parser.error("Cannot create file %s" % arg)




def main():
    main_parser = argparse.ArgumentParser(description='run Gene Ontology')
    subparsers = main_parser.add_subparsers(dest='which', help='type of enrichment analysis')
    subparsers.required = True
   
    parser = argparse.ArgumentParser(add_help=False)
    
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-o', '--out', type=str, required=True, nargs = '+',
                   help='output file')
    
    required.add_argument('-i', '--inp', type=str, required=True,
                   help='input gene list file')
    required.add_argument('-a', '--assoc', type=str, required=True,
                   help='input associations file (.gaf)')
    required.add_argument('-g', '--gograph', type=str, required=True,
                   help='input GO graph file (.obo)')
    
    
    parser.add_argument('-f', '--outputformat', choices=["html","txt", "gml", "png", "tabular"],  nargs = "+",
                   help='output file format', default = ["html"])
    
    parser.add_argument('-c', '--corrections', choices=["bonferroni","bh_fdr", "bonferroni,bh_fdr", "bh_fdr,bonferroni"],  
                   help='multiple hypothesis testing corrections', nargs='+', default=[])
    
    
    parser1 = subparsers.add_parser("term-for-term", parents=[parser])
    parser2 = subparsers.add_parser("parent-child", parents=[parser])
    
    

    #Parent-child
    parser2.add_argument('-m', '--method', choices=["union", "intersection"],  
                   help='method used to compute probabilities', default = "union")
    
    
    #validate args    
    if len(sys.argv) < 2:
        main_parser.print_usage()
        sys.exit(1)
    
    args = main_parser.parse_args()
    if len(args.out) != len(args.outputformat):
        main_parser.error("Number of output files doesn't match numer of formats!")
    check_file(main_parser, args.inp, 'r')
    check_file(main_parser, args.assoc, 'r')
    check_file(main_parser, args.gograph, 'r')
    for f in args.out:
        check_file(main_parser, f, 'w+')
    
    cors = []
    for cor in args.corrections:
        if "," in cor:
            cors += cor.split(",")
        else:
            cors.append(cor)
    args.corrections = list(set(cors))
    
    #Read inputs
    gene_list = read_list(args.inp)
    
    go_graph = OntoIO.read(args.gograph, "obo")
    #assocs = OntoIO.read(args.assoc, "gaf", assoc_format = "in_mem_sql")
    assocs = OntoIO.read(args.assoc, "gaf")
    result=None
    
        
    if args.which == "term-for-term":
        result = run_term(assocs, go_graph, gene_list, args.corrections)
    elif args.which == "parent-child":
        result = run_parent_child(assocs, go_graph, gene_list, args.corrections, args.method)
    else:
        parser.error("Method unimplemented!")

    assert result!= None,  "An error occured while computing result"
    print result
    for outfilename, outputformat in zip(args.out, args.outputformat):
        with open(outfilename, 'w+') as outfile:
            if outputformat == 'html':
                OntoIO.pretty_print(result, go_graph, outfile, outputformat, go_to_url="http://amigo.geneontology.org/amigo/term/")
            else:
                OntoIO.pretty_print(result, go_graph, outfile, outputformat)
                
        
        
    
    
    
if __name__ == "__main__":
   main()












