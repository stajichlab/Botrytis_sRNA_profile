#!/usr/bin/env python3
import csv, os, sys, re, subprocess
import io
import argparse
import pysam

parser = argparse.ArgumentParser(
    description='Get sRNA size distribution and info')

parser.add_argument('-b', '--bamfile', required=True, help='bamfile for input')
parser.add_argument('-o', '--out', required=False, help='outfile base name')
parser.add_argument('--maxrecords', required=False,default=-1, type=int,
                    help='maximum number of records to process')
parser.add_argument('--maxsize', required=False, default=30, help='max read size to report')
parser.add_argument('--minsize', required=False, default=12, help='min read size to report')
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

lengths = {}
n=0
names = {}
for read in iter:
    (srr,readid) = read.qname.split(".")
    if srr not in names:
        names[srr] = 0
#    print(read.qname,read.query_alignment_length,len(read.query_sequence))
#    print(read.query_alignment_start)
#    print(read.query_sequence)
#    print(read.query_alignment_sequence)
    if read.query_alignment_length not in lengths:
        lengths[read.query_alignment_length] = {srr: 0}
    elif srr not in lengths[read.query_alignment_length]:
        lengths[read.query_alignment_length][srr] = 0
    lengths[read.query_alignment_length][srr] += 1
    names[srr] += 1
    if (args.maxrecords > 0) and (n > args.maxrecords):
        break
    n += 1

with open("{}.read_lengths_counts.tsv".format(outname),"w") as outfh, open("{}.read_lengths_percent.tsv".format(outname),"w") as outpfh:
    outtsv=csv.writer(outfh,delimiter="\t")
    outtsvP=csv.writer(outpfh,delimiter="\t")

    srrorder = sorted(names.keys())
    header   = []
    #  replace SRR with names if we gave a tab files of SRR to name
    if len(sra2name) > 0:
        for n in sorted(srrorder,key=lambda x: sra2name[x] if x in sra2name else x):
            if n in sra2name:
                header.append(sra2name[n])
            else:
                header.append(n)
    else:
        header = srrorder
    header.insert(0,"SIZE")
    outtsv.writerow(header)
    outtsvP.writerow(header)
    for size in sorted(lengths):
        if size > args.maxsize or size < args.minsize:
            continue
        row = [size]
        rowP = [size]
        for n in srrorder:
            if n in lengths[size]:
                row.append(lengths[size][n])
                rowP.append("%.2f"%( 100*lengths[size][n] / names[n]))
            else:
                row.append(0)
                rowP.append("%.2f"%(0))
        outtsv.writerow(row)
        outtsvP.writerow(rowP)
