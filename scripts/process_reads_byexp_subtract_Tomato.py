#!/usr/bin/env python3
"""Count Reads overlapping features, stratefy by experiment which is groups of reads."""
import argparse
import csv
import gzip
import io
import logging
import os

from pybedtools import BedTool

parser = argparse.ArgumentParser(
    description='Subtract Tomato perfect matches from Bc perfect matches')

parser.add_argument('-o', '--outdir', default='results', help="output folder for results")
parser.add_argument('-b', '--botrytis', default="results/Botrytis_matchAligned.out.reads.tsv.gz",
                    help='Botrytis mapping reads in TSV format, separated by library - generated by get_BAM_sRNA_unique.py')

parser.add_argument('-t', '--tomato', default="results/Tomato_matchAligned.out.reads.tsv.gz",
                    help='Tomato mapping reads in TSV format, separated by library - generated by get_BAM_sRNA_unique.py')


parser.add_argument('-s', '--skiplow', action=argparse.BooleanOptionalAction, default=2, required=False,
                    help='Skip low count reads')

parser.add_argument('-f', '--featurefile', default="genomes/botrytis/BcinereaB05-10.features.sort.gff3.gz",
                    help='gene and TEfile for read overlap classification')

parser.add_argument('-v', '--debug', default=False, action=argparse.BooleanOptionalAction, help='Print debugging messages')
parser.add_argument('-l', '--log', default="process_reads.log",
                    help='Log file')

args = parser.parse_args()
if args.debug:
    logging.basicConfig(filename=args.log, encoding='utf-8', level=logging.DEBUG)
else:
    logging.basicConfig(filename=args.log, encoding='utf-8', level=logging.INFO)

outsave = os.path.join(args.outdir, "{}.tomato_subtract.tsv.gz".format("Botrytis"))
outsavebedF = os.path.join(args.outdir, "{}.F.tomato_subtract.bed".format("Botrytis"))
outsavebedR = os.path.join(args.outdir, "{}.R.tomato_subtract.bed".format("Botrytis"))

exprnames = set()
expr_count = {}
id2size = {}

with gzip.open(args.botrytis, "r") as bc,  gzip.open(args.tomato, "r") as sl, gzip.open(outsave, 'w') as outB, open(outsavebedF, 'w') as outBedF, open(outsavebedR, 'w') as outBedR:
    slreader = csv.reader(io.TextIOWrapper(sl, newline=""), delimiter="\t")
    slheader = next(slreader)

    slhits = set()

    for row in slreader:
        name = row[4]
        if name in slhits:
            print(f"found a non-unique sRNA read {name}")
        else:
            slhits.add(name)
    logging.debug(f"Done processing {args.tomato}")

    bcreader = csv.reader(io.TextIOWrapper(bc, newline=""), delimiter="\t")
    bcwriter = csv.writer(io.TextIOWrapper(outB, newline="", write_through=True), delimiter="\t")
    bcbedF = csv.writer(outBedF, delimiter="\t")
    bcbedR = csv.writer(outBedR, delimiter="\t")

    bcheader = next(bcreader)
    bcwriter.writerow(bcheader)

    col2expr = {}

    exprcolstart = 8
    s = 8

    for col in bcheader[exprcolstart:]:
        expr = col[0:1]  # first letter of library name is experiment
        col2expr[s] = expr  # match column number to experimental group
        exprnames.add(expr)
        s += 1

    for row in bcreader:
        name = row[4]
        if int(row[7]) <= args.skiplow:
            continue

        if name not in slhits:
            bcwriter.writerow(row)
            expr_count[name] = {}
            for ename in exprnames:
                expr_count[name][ename] = 0

            for colct in range(exprcolstart, len(row)):
                expr = col2expr[colct]
                expr_count[name][expr] += int(row[colct])

            id2size[name] = {'length': int(row[5]),   # this LENGTH col
                             'unique': int(row[6]),   # this is UNIQUE col
                             'count': int(row[7])}    # this TOTAL_COUNT column
            if row[4] == "+":
                bcbedF.writerow([row[0], row[1], row[2], row[3]])
            else:
                bcbedR.writerow([row[0], row[1], row[2], row[3]])

logging.debug(f"Done processing {args.botrytis}")

# forward read matching
genomefeatures = BedTool(args.featurefile)
sRNAF = BedTool(outsavebedF)
sRNAR = BedTool(outsavebedR)

# this would only be reads that overlap we need to also consider those which do not overlap
readsInFeaturesF = sRNAF.intersect(genomefeatures, wao=True)

logging.debug("Processing fwd strand features")

sizeprofile = {}
sizeprofileuniq = {}

sizeprofileexpr = {}
sizeprofileexpruniq = {}

for ename in exprnames:
    sizeprofileexpr[ename] = {}
    sizeprofileexpruniq[ename] = {}

types = set()

# for now I am repeating these routines as I don't quite want to
# make this generic

for i in readsInFeaturesF:
    seq = i[3]
    strand = i[4]
    print(i)
    read_feature_direction = 'sense'
    if strand == "-":
        read_feature_direction = 'antisense'
    len = id2size[seq]['length']
    unique = id2size[seq]['unique']
    ct = id2size[seq]['count']

    type = i[6]
    TypeClass = "None"
    if type == "exon":
        TypeClass = type
    elif type == "tRNA" or type == "intron":
        TypeClass = type
    elif type == 'match':  # this is how we coded RepeatMasker results
        grpcol = i[12]
        grp = {}
        for nm in grpcol.split(";"):
            (key, val) = nm.split("=")
            v = val.split('/')[0]
            grp[key] = f"TE.{v}"
        TypeClass = grp['type']
    elif type == ".":
        TypeClass = "None"
    else:
        continue

    types.add(TypeClass)
    if len not in sizeprofile:
        sizeprofile[len] = {TypeClass: {read_feature_direction: ct}}
    elif TypeClass not in sizeprofile[len]:
        sizeprofile[len][TypeClass] = {read_feature_direction: ct}
    elif read_feature_direction not in sizeprofile[len][TypeClass]:
        sizeprofile[len][TypeClass][read_feature_direction] = ct
    else:
        sizeprofile[len][TypeClass][read_feature_direction] += ct

    if unique == 1:
        if len not in sizeprofileuniq:
            sizeprofileuniq[len] = {TypeClass: {read_feature_direction: ct}}
        elif TypeClass not in sizeprofileuniq[len]:
            sizeprofileuniq[len][TypeClass] = {read_feature_direction: ct}
        elif read_feature_direction not in sizeprofileuniq[len][TypeClass]:
            sizeprofileuniq[len][TypeClass][read_feature_direction] = ct
        else:
            sizeprofileuniq[len][TypeClass][read_feature_direction] += ct

    # now do per-experiment count
    for expr in exprnames:
        ecount = expr_count[name][expr]
        if len not in sizeprofileexpr[expr]:
            sizeprofileexpr[expr][len] = {TypeClass: {read_feature_direction: ecount}}
        elif TypeClass not in sizeprofileexpr[expr][len]:
            sizeprofileexpr[expr][len][TypeClass] = {read_feature_direction: ecount}
        elif read_feature_direction not in sizeprofileexpr[expr][len][TypeClass]:
            sizeprofileexpr[expr][len][TypeClass][read_feature_direction] = ecount
        else:
            sizeprofileexpr[expr][len][TypeClass][read_feature_direction] += ecount

# now process reverse strand features and check if sRNA is on same strand as the feature
readsInFeaturesR = sRNAR.intersect(genomefeatures, wao=True)
logging.debug("Processing rev strand features")

for i in readsInFeaturesR:
    seq = i[3]
    strand = i[4]
    print(i)
    read_feature_direction = 'sense'
    if strand == "+":
        read_feature_direction = 'antisense'
    len = id2size[seq]['length']
    unique = id2size[seq]['unique']
    ct = id2size[seq]['count']

    type = i[6]
    TypeClass = "None"
    if type == "exon" or type == "tRNA" or type == "intron":
        TypeClass = type
    elif type == 'match':  # this is how we coded RepeatMasker results
        grpcol = i[12]
        grp = {}
        for nm in grpcol.split(";"):
            (key, val) = nm.split("=")
            v = val.split('/')[0]
            grp[key] = f"TE.{v}"
        TypeClass = grp['type']
    elif type == ".":
        TypeClass = "None"
    else:
        continue

    types.add(TypeClass)
    if len not in sizeprofile:
        sizeprofile[len] = {TypeClass: {read_feature_direction: ct}}
    elif TypeClass not in sizeprofile[len]:
        sizeprofile[len][TypeClass] = {read_feature_direction: ct}
    elif read_feature_direction not in sizeprofile[len][TypeClass]:
        sizeprofile[len][TypeClass][read_feature_direction] = ct
    else:
        sizeprofile[len][TypeClass][read_feature_direction] += ct

    if unique == 1:
        if len not in sizeprofileuniq:
            sizeprofileuniq[len] = {TypeClass: {read_feature_direction: ct}}
        elif TypeClass not in sizeprofileuniq[len]:
            sizeprofileuniq[len][TypeClass] = {read_feature_direction: ct}
        elif read_feature_direction not in sizeprofileuniq[len][TypeClass]:
            sizeprofileuniq[len][TypeClass][read_feature_direction] = ct
        else:
            sizeprofileuniq[len][TypeClass][read_feature_direction] += ct

    # now do per-experiment count
    for expr in exprnames:
        ecount = expr_count[name][expr]
        if len not in sizeprofileexpr[expr]:
            sizeprofileexpr[expr][len] = {TypeClass: {read_feature_direction: ecount}}
        elif TypeClass not in sizeprofileexpr[expr][len]:
            sizeprofileexpr[expr][len][TypeClass] = {read_feature_direction: ecount}
        elif read_feature_direction not in sizeprofileexpr[expr][len][TypeClass]:
            sizeprofileexpr[expr][len][TypeClass][read_feature_direction] = ecount
        else:
            sizeprofileexpr[expr][len][TypeClass][read_feature_direction] += ecount

for desc in ['sense', 'antisense', 'combined']:
    ctheader = ['SIZE']
    ctheader.extend(sorted(list(types)))
    ctheader.append("TOTAL")

    with open(os.path.join(args.outdir, f"Botrytis_size_{desc}_profile.tsv"), "w") as ofh:
        bcwriter = csv.writer(ofh, delimiter="\t")
        bcwriter.writerow(ctheader)

        for size in sorted(sizeprofile):
            row = [size]
            sum = 0
            for type in sorted(types):
                ct = 0
                if type in sizeprofile[size]:
                    if desc == "sense" or desc == "antisense":
                        ct = sizeprofile[size][type][desc]
                    else:
                        ct = sum(sizeprofile[size][type].values())
                row.append(ct)
                sum += ct
            row.append(sum)
            bcwriter.writerow(row)

    with open(os.path.join(args.outdir, "Botrytis_size_combined_profile_uniq.tsv"), "w") as ofh:
        bcwriter = csv.writer(ofh, delimiter="\t")
        bcwriter.writerow(ctheader)
        for size in sorted(sizeprofileuniq):
            row = [size]
            sum = 0
            for type in sorted(types):
                ct = 0
                if type in sizeprofileuniq[size]:
                    if desc == "sense" or desc == "antisense":
                        ct = sizeprofileuniq[size][type][desc]
                    else:
                        ct = sum(sizeprofileuniq[size][type].values())
                row.append(ct)
                sum += ct
            row.append(sum)
            bcwriter.writerow(row)
