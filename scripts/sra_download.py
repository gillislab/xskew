#!/usr/bin/env python
import argparse
import os
import sys
import pandas as pd
import logging

gitpath = os.path.expanduser("~/git/scqc")
sys.path.append(gitpath)
from scqc.utils import *

def parse_tsv(infile, species):
    rdf = load_df(infile)
    srcurls = rdf.loc[rdf.sciname == species ,['file_url']]
        for i in range(srcurls.shape[0]):
            try:
                download_wget(srcurls.file_url[i], destpath=f'{destdir}/{srcurls.run_id[i]}.sra',rate = '50M')
            except KeyError:
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

    parser.add_argument('-i', '--infile', 
                        metavar='infile', 
                        type=str, 
                        help='a runs.tsv file with file_url and sciname columns ')

    parser.add_argument('-s', '--species', 
                        metavar='species', 
                        type=str, 
                        help='Linnean species in quotes, e.g.  Sus scrofa')    

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)


