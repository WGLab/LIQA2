#!/usr/bin/python

from __future__ import print_function  # load print function in python3
from collections import defaultdict
import math, sys, os, re, pysam, time
from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()

# set up auto dictionary function
def auto_dict():
    return defaultdict(auto_dict)

###############################################################################
###  ARGUMENT SETTINGS
###############################################################################

# checking whether argument is valid or not
validArgList = ["-bam", "-ref", "-out"]
for argIndex in range(1, len(sys.argv)):
    if sys.argv[argIndex][0] == "-" and sys.argv[argIndex] not in validArgList:
        print("Argument \'"+sys.argv[argIndex]+"\' is invalid!")
        sys.exit()

bamFileExists = 0
refFileExists = 0
outFileExists = 0

for argIndex in range(1, len(sys.argv)):
    if sys.argv[argIndex] == "-bam":  ## load in BAM file
        argIndex += 1
        bamFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        bamTmp = sys.argv[argIndex].split("/")
        bamFile = bamFileAbsPath + "/" + bamTmp[len(bamTmp)-1]
        bamFileExists = 1
    elif sys.argv[argIndex] == "-ref":  ## load in annotation file
        argIndex += 1
        refFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        refTmp = sys.argv[argIndex].split("/")
        refGeneFile = refFileAbsPath + "/" + refTmp[len(refTmp)-1]
        refFileExists = 1
    elif sys.argv[argIndex] == "-out":  ## load in annotation file
        argIndex += 1
        outFileAbsPath = os.path.dirname(os.path.abspath(sys.argv[argIndex]))
        outTmp = sys.argv[argIndex].split("/")
        outFile = outFileAbsPath + "/" + outTmp[len(outTmp)-1]
        outFileExists = 1

if bamFileExists == 0 or refFileExists == 0 or outFileExists == 0:
    print("Please provide arguments:")
    print("-bam\tIndexed bam file")
    print("-ref\tGene annotation file")
    print("-out\tOutput file")
    sys.exit()

# load gene information
geneStructureInformation = auto_dict()
geneLineCount = auto_dict()

with open(refGeneFile, "r") as FP:
    for line in FP:
        line = line.strip("\n")
        tmpinf = line.split("\t")
        gene = tmpinf[0]

        if not bool(geneStructureInformation[gene]):
            geneLineCount[gene] = 0
            geneStructureInformation[gene][geneLineCount[gene]] = line
        else:
            geneLineCount[gene] += 1
            geneStructureInformation[gene][geneLineCount[gene]] = line

#############################################
## Calculate coverage and support thresholds
#############################################

def calculate_coverage_threshold():
    total_exons = 0
    total_reads = 0
    for gene in geneStructureInformation:
        numofExons = geneLineCount[gene]
        total_exons += numofExons
        for i in range(1, numofExons + 1):
            total_reads += count_reads_in_exon(gene, i)
    average_reads_per_exon = total_reads / total_exons
    coverage_threshold = int(average_reads_per_exon * 0.8)  # Adjust as needed
    return coverage_threshold

def calculate_support_threshold():
    all_read_counts = []
    for gene in geneStructureInformation:
        numofExons = geneLineCount[gene]
        for i in range(1, numofExons + 1):
            read_count = count_reads_in_exon(gene, i)
            all_read_counts.append(read_count)
    support_threshold = sorted(all_read_counts)[len(all_read_counts) // 2]  # Median - also adjust as needed
    return support_threshold

def count_reads_in_exon(gene, exon_index):
    exon_info = geneStructureInformation[gene][exon_index].split("\t")
    exon_start = int(exon_info[3])
    exon_end = int(exon_info[4])
    exon_read_count = 0
    for read in bamFilePysam.fetch(geneChr, exon_start, exon_end):
        exon_read_count += 1
    return exon_read_count

# Calculate coverage and support thresholds
coverage_threshold = calculate_coverage_threshold()
support_threshold = calculate_support_threshold()

###############################################################################
###  PROCESSING DATA FOR EACH GENE
###############################################################################

# Open output file
OUT = open(outFile, 'w')

# Write header
OUT.write("GeneName\tIsoformName\tReadPerGene_corrected\tRelativeAbundance\tinfor_ratio\n")

# Iterate over each gene
for gene in geneStructureInformation:
    geneCount += 1
    tmpTime = (time.time() - startTime) / 60.0

    sameReadCount = auto_dict()
    readStart = auto_dict()
    readEnd = auto_dict()
    readCigar = auto_dict()
    noveliso_ct = auto_dict()
    noveliso = auto_dict()

    numofExons = geneLineCount[gene]
    tmpgeneinf = geneStructureInformation[gene][0].split("\t")
    geneChr = tmpgeneinf[1]
    geneStart = int(tmpgeneinf[3])
    geneEnd = int(tmpgeneinf[4])

    # Calculate gene and isoform length information
    tmpisoinf = tmpgeneinf[5].split(";")
    tmpgeneinf[5] = tmpisoinf[0]
    if len(tmpisoinf) == 2:
        tmpisolength = tmpisoinf[1].split(",")
        genelength = int(tmpisolength[0]) - geneStart
        for iii in range(len(tmpisolength) - 1):
            tmpisolength[iii] = (int(tmpisolength[iii]) - geneStart)
            if tmpisolength[iii] == 666:
                tmpisolength[iii] = 0
            else:
                tmpisolength[iii] = tmpisolength[iii] / 2566
    else:
        genelength = 100
        tmpisolength = tmpisoinf[0].split(",")
        for iii in range(len(tmpisolength) - 1):
            tmpisolength[iii] = 0

    ## Load all reads information mapped to the specific gene within this loop using pysam
    for read in bamFilePysam.fetch(geneChr, geneStart, geneEnd):
        line = str(read)
        tmpinf = line.split("\t")
        tmpReadName = tmpinf[0]
        tmpReadChr = geneChr
        tmpReadStart = int(tmpinf[3]) + 1
        tmpReadCigar = ""

        # Adjust to different Pysam Version
        if ")]" in tmpinf[5]:  # vector format
            tmpinf[5] = tmpinf[5].rstrip(")]")
            tmpinf[5] = tmpinf[5].lstrip("[(")
            tmpinfcigar = tmpinf[5].split("), (")
            for cc in tmpinfcigar:
                ttcc = cc.split(", ")
                if ttcc[0] == "3":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "N"
                if ttcc[0] == "2":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "D"
                if ttcc[0] == "1":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "I"
                if ttcc[0] == "0":
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "M"
                if not (ttcc[0] == "3" or ttcc[0] == "2" or ttcc[0] == "1" or ttcc[0] == "0"):
                    tmpReadCigar = tmpReadCigar + ttcc[1] + "X"
        else:  # 100M10N100M format
            tmpReadCigar = tmpinf[5]

        if not bool(sameReadCount[tmpReadName]):
            sameReadCount[tmpReadName] = 1
        else:
            sameReadCount[tmpReadName] += 1

        readStart[tmpReadName][sameReadCount[tmpReadName]] = tmpReadStart
        readCigar[tmpReadName][sameReadCount[tmpReadName]] = tmpReadCigar.replace('=', "M")

    # Load structure information of the specific gene within this loop
    tmpgeneinf[5] = tmpgeneinf[5].rstrip(",")
    isoformNames = tmpgeneinf[5].split(",")
    exonStarts = [None] * numofExons
    exonEnds = [None] * numofExons
    exonIndicators = auto_dict()

    for i in range(1, numofExons + 1):
        tmpinf = geneStructureInformation[gene][i].split("\t")
        exonStarts[i - 1] = int(tmpinf[3]) + 1
        exonEnds[i - 1] = int(tmpinf[4])
        tmpinf[5] = tmpinf[5].rstrip(",")
        tmpExonIndicators = tmpinf[5].split(",")

        for j in range(len(tmpExonIndicators)):
            exonIndicators[isoformNames[j]][i - 1] = int(tmpExonIndicators[j])

    lociIndicators = auto_dict()
    for i in range(len(isoformNames)):
        for j in range(len(exonStarts)):
            if exonIndicators[isoformNames[i]][j] == 1:
                for k in range(exonStarts[j], exonEnds[j] + 1):
                    lociIndicators[isoformNames[i]][k] = 1

    #########################################################################################################################################
    ## START TO ANALYZE EACH READ
    ##################################################################################################################################################

    qualifiedRead = auto_dict()
    prereadCount = 0
    readCount = 0
    readCount1iso = 0
    fragmentStart = auto_dict()
    fragmentEnd = auto_dict()
    CompatibleMatrix = auto_dict()
    tmpCompatibleMatrix = auto_dict()
    qualitycheck = auto_dict()
    readEnd1 = auto_dict()
    readEnd2 = auto_dict()
    Tforkmf = []
    Eforkmf = []

    for readName in sameReadCount:

        # load CIGAR information
        cigarNumberRead1 = auto_dict()
        cigarNumberRead2 = auto_dict()
        cigarMatchRead1 = auto_dict()
        cigarMatchRead2 = auto_dict()
        cigarInfCountRead1 = 0
        cigarInfCountRead2 = 0
        cigarInfCountRead1tmp = 0
        cigarInfCountRead2tmp = 0
        qualitycheck[readName] = 0
        readEnd1[readName] = 0
        readEnd2[readName] = 1

        tmp1 = re.split("([A-Z])", readCigar[readName][1])
        for i in range(len(tmp1) - 1):
            if tmp1[i].isalpha():
                cigarMatchRead1[cigarInfCountRead1] = tmp1[i]
                cigarInfCountRead1 += 1
            else:
                cigarNumberRead1[cigarInfCountRead1] = int(tmp1[i])
                cigarInfCountRead1tmp += 1

        if sameReadCount[readName] == 2:
            tmp2 = re


# Close output file
OUT.close()
