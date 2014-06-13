#!/usr/bin/env python
import gem
from gem.pipeline import *
import argparse
import sys
import logging

logging.basicConfig()
log = logging.getLogger('denovo')
#gem.loglevel('info')

def filter_map(args):
    all_mappings = gem.files.open(args.input)
    filtered = gem.filter.unmapped(all_mappings, eval(args.mismatches))

    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    for line in filtered:
        if args.format == 'map':
            print line.to_map()
        if args.format == 'fastq':
            print line.to_fastq()

    #exit(0)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest="input", help='Input file with all the mappings', required=True)
    parser.add_argument('-t', '--threads', dest="threads", default=1, help='Number of threads')
    parser.add_argument('-f', '--format', dest="format", choices=['map', 'fastq'], help='Output file format')
    parser.add_argument('-m', '--mismatches', dest="mismatches", default="5", help='Mismatches to consider a read unmapped (>=)')

    args = parser.parse_args()
    filter_map(args)

