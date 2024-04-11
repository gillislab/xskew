#!/usr/bin/env python
import argparse
import os
import sys
import pandas as pd
import logging

gitpath = os.path.expanduser("~/git/xskew")
sys.path.append(gitpath)

from Bio import SeqIO
from Bio.Seq import Seq 

from xskew.tools import *
from genome.genome import *



'''
ENSEMBL
https://useast.ensembl.org/index.html
https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/
https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/


NCBI
https://www.ncbi.nlm.nih.gov/genome
https://www.ncbi.nlm.nih.gov/assembly/model/
https://www.ncbi.nlm.nih.gov/assembly/basics/

https://www.ncbi.nlm.nih.gov/genome/215  (Macaca mulatta)
  -> Reference genome = Macaca mulatta Mmul_10

Genbank (GCA) 
Genome:
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Macaca_mulatta/ 
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Macaca_mulatta/latest_assembly_versions/GCA_003339765.3_Mmul_10/

Annotation:
https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9544/103/GCF_003339765.1_Mmul_10/


Refseq (GCF)
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Macaca_mulatta/
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Macaca_mulatta/annotation_releases/103/GCF_003339765.1_Mmul_10/

GRC genome reference consortium (human and mouse genomes)

'''

if __name__ == "__main__":
    FORMAT = '%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--debug',
                        action="store_true",
                        dest='debug',
                        help='debug logging')

    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        dest='verbose',
                        help='verbose logging')


    parser.add_argument('-s', '--source', 
                        metavar='source', 
                        type=str,
                        required=True, 
                        help='genome source [ensembl|refseq|genbank] refseq=GCF genbank=GCA')  

    parser.add_argument('-g', '--genome', 
                        metavar='genome', 
                        type=str,
                        required=True, 
                        help='fasta genome file')

    parser.add_argument('-a', '--annotation', 
                        metavar='annotation', 
                        type=str, 
                        help='GTF annotations. ')

    parser.add_argument('-r', '--report', 
                        metavar='report', 
                        type=str, 
                        help='assembly report file  ')

    parser.add_argument('-o', '--outdir', 
                        metavar='outdir',
                        default=None, 
                        type=str, 
                        help='output directory, normally same as where genome file is')
  

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)


    outdir = args.outdir
    if args.outdir is None:
        outdir = os.path.dirname(os.path.abspath(os.path.expanduser(args.genome)))
        
    logging.debug(f'source = {args.source} genome={args.genome} report={args.report} annot={args.annotation} outdir={outdir} ')
    
    if args.source == 'ensembl':
        prepare_genome_ensembl(args.genome, args.annotation, outdir)
    elif args.source == 'refseq':
        prepare_genome_refseq(args.genome, args.annotation, args.report, outdir)
        
    