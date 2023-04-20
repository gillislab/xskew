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


def prepare_genome_ensembl(genomefile, annotfile, outdir):
    '''
    for ensemble, no assembly report is needed. 
    make <chrnum>label.txt files for each chromosome, with <chrnum> as value inside. 
    make <chrnum>.fa and <chrnum>.fa.fai for each chromosome using index_region
        
     
    '''
    logging.debug(f'creating symlinks')
    if os.path.exists(f'{outdir}/genome.fa'):
        os.remove(f'{outdir}/genome.fa')
    os.symlink(genomefile, f'{outdir}/genome.fa')
    if os.path.exists(f'{outdir}/annotation.gtf'):
        os.remove(f'{outdir}/annotation.gtf')
    os.symlink(annotfile, f'{outdir}/annotation.gtf')

    sfh = open(genomefile)
    seq_dict = SeqIO.to_dict(SeqIO.parse(sfh, "fasta"))
    sfh.close()
    logging.debug(f'read genome file {genomefile} w/ {len(seq_dict)} records. keys={seq_dict.keys()}')
    chrlabels = list(seq_dict.keys())
    for chrlabel in chrlabels:
        fn = f'{outdir}/{chrlabel}label.txt'
        logging.debug(f'writing {fn}...')
        with open(fn, 'w') as f:
            f.write(f'{chrlabel}\n')
    
        make_chr_index( infile=genomefile,
                        genomedir = outdir, 
                        chr = chrlabel, 
                        outfile = f'{outdir}/{chrlabel}.fa')
        logging.debug(f'handled chromsome {chrlabel}')
    
    star_genome(outdir, "6" ,f'{outdir}/annotation.gtf',  f'{outdir}/genome.fa')
    samtools_dict(f'{outdir}/genome.fa',f'{outdir}/genome.dict')
    
    logging.info(f'done.')
            

def parse_assembly_report(reportfile):
    '''
    return list of tuples  ( chrnum,   contigid )
    
    '''
    colnames = ['Sequence-Name','Sequence-Role','Assigned-Molecule',
                'Assigned-Molecule-Location/Type','GenBank-Accn',
                'Relationship','RefSeq-Accn','Assembly-Unit',
                'Sequence-Length','UCSC-style-name']   
    df = pd.read_csv(reportfile, comment="#", sep='\t', header=None)
    df.columns = colnames
    amdf = df[df['Sequence-Role'] == 'assembled-molecule']
    amdf = amdf[amdf['Assembly-Unit'] == 'Primary Assembly' ]
    mdf = amdf[['Assigned-Molecule','RefSeq-Accn']]
    tlist = list(mdf.itertuples(index=False, name=None))
    logging.debug(f'extracted {len(tlist)} maps in {reportfile}')
    return tlist

def prepare_genome_refseq(genomefile, annotfile, reportfile, outdir):
    logging.debug(f'creating symlinks')
    if os.path.exists(f'{outdir}/genome.fa'):
        os.remove(f'{outdir}/genome.fa')
    os.symlink(genomefile, f'{outdir}/genome.fa')
    if os.path.exists(f'{outdir}/annotation.gtf'):
        os.remove(f'{outdir}/annotation.gtf')
    os.symlink(annotfile, f'{outdir}/annotation.gtf')
    
    tlist = parse_assembly_report(reportfile)
    logging.debug(tlist)
    for (chrlabel, contig) in tlist:
        fn = f'{outdir}/{chrlabel}label.txt'
        logging.debug(f'writing {fn}...')
        with open(fn, 'w') as f:
            f.write(f'{contig}\n')

        make_chr_index( infile=genomefile,
                        genomedir = outdir, 
                        chr = chrlabel, 
                        outfile = f'{outdir}/{chrlabel}.fa')
        logging.debug(f'handled chromsome {chrlabel}')
        
    star_genome(outdir,"6" ,f'{outdir}/annotation.gtf',  f'{outdir}/genome.fa')
    samtools_dict(f'{outdir}/genome.fa',f'{outdir}/genome.dict')
    samtools_faidx(f'{outdir}/genome.fa',f'{outdir}/genome.fa.fai')
    logging.info(f'done.')


def prepare_genome_genbank():
    pass








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
        
    