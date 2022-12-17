#!/usr/bin/env python3
import csv, os, sys, re
import io
import gzip
import argparse
import pysam

#from sqlitedict import SqliteDict
#import tempfile

parser = argparse.ArgumentParser(
    description='Get sRNA size distribution and info')

parser.add_argument('-b', '--bamfile', required=True, help='bamfile for input')
parser.add_argument('-o', '--out', required=False, help='outfile base name')
parser.add_argument('-v','--debug',default=False,required=False,help='debugging flag')
parser.add_argument('--maxrecords', required=False,default=-1, type=int,
                    help='maximum number of records to process')
# add min coverag
parser.add_argument('--tempdir', required=False, default='/scratch', help='scratch folder')
parser.add_argument('--maxsize', required=False, default=30, help='max read size to report')
parser.add_argument('--minsize', required=False, default=18, help='min read size to report')
parser.add_argument('--sra2name', required=False, type=argparse.FileType('r'),
                    help='map of sranames to experiment names')
args = parser.parse_args()
outname=""
if args.out == None:
    outname=os.path.splitext(args.bamfile)[0]
else:
    outname=args.out

sra2name = {}
if args.sra2name != None:
    csvin = csv.reader(args.sra2name,delimiter="\t")
    for r in csvin:
        sra2name[r[0]] = r[1]
samfile = pysam.AlignmentFile(args.bamfile, "rb")
iter = samfile.fetch()

n=0
librarynames = {}
perfectmatch = re.compile(r'^(\d+)M$')
sRNA = {}
sRNAInfo = {}
sRNABlackList = set()

for read in iter:
    (srr,readid) = read.qname.split(".")
    if srr not in librarynames:
        librarynames[srr] = 0
    seqlen = len(read.query_alignment_sequence)
    m=perfectmatch.match(read.cigarstring)
    if  seqlen > args.maxsize or seqlen < args.minsize or not m:
        #      print("skipping {} as it isn't a perfect match ({})".format(read.query_alignment_sequence,read.cigarstring))
        continue
    if read.query_alignment_sequence not in sRNA:
        sRNA[read.query_alignment_sequence] = {srr: 0}
    elif srr not in sRNA[read.query_alignment_sequence]:
        sRNA[read.query_alignment_sequence][srr] = 0

    sRNA[read.query_alignment_sequence][srr] += 1
    if read.query_alignment_sequence not in sRNAInfo:
        strand = "+"
        if read.is_reverse:
            strand = '-'
        sRNAInfo[read.query_alignment_sequence] = { 'ref': read.reference_id,
                                                    'start'         : read.reference_start,
                                                    'end'           : read.reference_end,
                                                    'strand'        : strand}
    elif sRNAInfo[read.query_alignment_sequence]['ref'] != read.reference_id or sRNAInfo[read.query_alignment_sequence]['start'] != read.reference_start:
        sRNABlackList.add(read.query_alignment_sequence)

    librarynames[srr] += 1
    if (args.maxrecords > 0) and (n > args.maxrecords):
        break
    n += 1
    if (n % 10**6) == 0: # counter
        print("Processed {} records".format(n))

skipcount = 0
with gzip.open("{}.reads.tsv.gz".format(outname),"w") as outfh:
    outtsv=csv.writer(io.TextIOWrapper(outfh,newline="", write_through=True),delimiter="\t")

    srrorder = sorted(librarynames.keys())
    header   = ['CHROM','START','END','STRAND','SEQ','LENGTH','UNIQUE','TOTAL_COUNT']
    #  replace SRR with names if we gave a tab files of SRR to name
    if len(sra2name) > 0:
        for n in sorted(srrorder,key=lambda x: sra2name[x] if x in sra2name else x):
            if n in sra2name:
                header.append(sra2name[n])
            else:
                header.append(n)
    else:
        header.extend(srrorder)
    outtsv.writerow(header)

    table = []
    for seq in sRNA:
        unique = 1
        if seq in sRNABlackList:
            skipcount +=1
            unique = 0
        row = [samfile.get_reference_name(sRNAInfo[seq]['ref']),
               sRNAInfo[seq]['start'],sRNAInfo[seq]['end'], sRNAInfo[seq]['strand'], seq, len(seq),unique,sum(sRNA[seq].values()) ]
        for n in srrorder:
            if n in sRNA[seq]:
                row.append(sRNA[seq][n])
            else:
                row.append(0)
        table.append(row)

    for row in sorted(table,key=lambda x: x[2],reverse=True):
        outtsv.writerow(row)

print("There were {} reads that did not have a unique mapping location".format(skipcount))
