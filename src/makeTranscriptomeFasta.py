'''
Created on Jun 16, 2011

@author: helena_storvall
@Changed 2011-07-05: Start position is converted from zero-based to one-based
'''
import GenomeFetcher
import optparse
import os
import time
import json
import copy
from interval import Interval
import subprocess
import sys

def find_overlaps(chrom, lineList):
    lineList.sort()
    transDict = {}
    while lineList:
        startA, endA, transNameA = lineList[0]
        geneIntervalA = Interval.between(startA, endA)
        allTransId = [transNameA]
        olGenes = [lineList[0]]
        geneStart = copy.deepcopy(startA)
        geneEnd = copy.deepcopy(endA)
        # Iterate over all consecutive lines to the first line
        # If the transcripts come from the same chromosome and strand, and if their gene interval overlap,
        # then
        for tbLine in lineList[1:]:
            startB, endB, transNameB = tbLine
            geneIntervalB = Interval.between(startB, endB)
            if geneIntervalB.overlaps(geneIntervalA):
                if geneEnd < endB: #set geneEnd to the highest end
                    geneEnd = copy.deepcopy(endB)
                    geneIntervalA = Interval.between(geneStart, geneEnd)
                olGenes.append(tbLine)
                allTransId.append(transNameB)
            else:
                break                        
        for gene in olGenes: #Remove all lines corresponding to overlapping transcripts
            lineList.remove(gene)
        for transName in allTransId:
            if transDict.has_key(transName):
                print '1. Warning! Transcript %s occurs multiple times in the file'%transName
                transDict[transName] = max(transDict.get(transName), len(allTransId))
            else:
                transDict[transName] = len(allTransId)
    return transDict

def create_overlap_json(annFile, t, jsonFile):
    af = open(annFile)
    chrDict = {}
    nameDict = {}
    for line in af:
        tabs = line.strip('\n').split('\t')
        transName = tabs[t]
        chrom, strand = tabs[t+1:t+3]
        start, end = map(int, tabs[t+3:t+5])
        if chrom.count('_') == 0:
            chrDict.setdefault(chrom, []).append([start, end, transName])
    af.close()
    allTransDict = {}
    for chrom, lineList in chrDict.items():
        chromTransDict = find_overlaps(chrom, lineList)
        for transName, val in chromTransDict.items():
            if allTransDict.has_key(transName):
                print '2. Warning! Transcript %s occurs multiple times in the file'%transName
                allTransDict[transName] = max(allTransDict.get(transName), val)
            else:
                allTransDict[transName] = val
    jf = open(jsonFile, 'w+')
    json.dump(allTransDict, jf)
    jf.close()
    return allTransDict

def main(ass, annFile, annType, spec, fastaDir, mergedFile, jsonFile):
    if jsonFile is None:
        jsonFile = os.path.join(fastaDir, '%s_transcriptOverlaps.json'%annType)
    if (annType.split('+')[0].lower() == 'refgene' or annType.split('+')[0].lower() == 'ensgene' or annType.split('+')[0].lower() == 'mgcgenes'):
        t = 1
    elif annType.split('+')[0].lower() == 'knowngene':
        t = 0
    if os.path.isfile(jsonFile):
        jf = open(jsonFile)
        transDict = json.load(jf)
        jf.close()
    else:
        transDict = create_overlap_json(annFile, t, jsonFile)
    af = open(annFile)
    gf = GenomeFetcher.GenomeFetcher(species = spec, assembly = ass)
    if fastaDir is None:
        fastaDir = 'transcriptomeFastaFiles'
    if mergedFile is None:
        mergedFile = os.path.join(fastaDir, '%s.%s.transcriptome.merged.fa'%(ass, annType))
    if os.path.isfile(mergedFile):
        subprocess.call('rm %s'%mergedFile, shell = True)
    mf = open(mergedFile, 'a')
    allFileNames = []
    i = 0
    for line in af:
        i+=1
        tabs = line.strip('\n').split('\t')
        #print 'line %i'%i, tabs
        transName = tabs[t]
        chrom, strand = tabs[t+1:t+3]
        if chrom.endswith('random') or chrom.count('_') > 0: continue
        overlaps = transDict.get(transName, None)
        if overlaps == None:
            overlaps = 100
            print 'Warning! The transcript %s was not in the dictionary. Default of 100 overlaps will be used.'%transName
        outDir = os.path.join(fastaDir, 'transcripts_m_%i'%overlaps)
        if not os.path.isdir(outDir):
            os.mkdir(outDir) 
        fileName = os.path.join(outDir, '%s.fa'%transName)
        tf = open(fileName, 'w')
        falseExonStarts = map(int, tabs[t+8].rstrip(',').split(','))
        exonStarts = []
        for fes in falseExonStarts:#Changed 2011-07-05
            exonStarts.append(fes+1)
        exonEnds = map(int, tabs[t+9].rstrip(',').split(','))
        if strand == '-':
            exonStarts.reverse()
            exonEnds.reverse()
        tf.write('>%s\t%s\t%s\t%s\t%s\n'%(transName, chrom, strand,tabs[t+8].rstrip(','),tabs[t+9].rstrip(',')))
        mf.write('>%s\t%s\t%s\t%s\t%s\n'%(transName, chrom, strand,tabs[t+8].rstrip(','),tabs[t+9].rstrip(',')))
        
        for es, ee in zip(exonStarts, exonEnds):
            #print '(chrom: %s, es: %i, ee: %i, strand: %s)'%(chrom, es, ee, strand)
            seq = gf.get_seq_from_to(chrom, es, ee, strand)
            tf.write(seq)
            mf.write(seq+'\n')
        tf.close()
        if i%1000 == 0:
            print 'transcript %i processed %s'%(i, time.asctime())
    mf.close()
    #return mergedFile, fastaDir

        
    
if __name__=='__main__':
    opts = optparse.OptionParser()
    opts.formatter=optparse.TitledHelpFormatter()
    opts.add_option('-a','--assembly',dest='ass')
    opts.add_option('-n','--annotation-file',dest='annFile')
    opts.add_option('-t','--annotation-type', dest='annType')
    opts.add_option('-s','--species',dest='spec')
    opts.add_option('-o','--out-dir',dest='outDir', default = 'transcriptomeFastaFiles')
    opts.add_option('-j','--json-file',dest='jsonFile')
    opts.add_option('-m','--merge-name', dest='mergeName')
    opts.add_option('--print-help', dest = 'help', default = False, action = 'store_true')
    (o, args) = opts.parse_args()
    
    if not o.ass is None and not o.annFile is None and not o.annType is None and not o.spec is None:
        print 'Running '+time.asctime()
        main(o.ass, o.annFile, o.annType, o.spec, o.outDir, o.jsonFile, o.mergeName)
        print 'Done '+time.asctime()
    else:
        o.help = True
    
    if len(sys.argv) < 2 or o.help:
        print 'Non-optional arguments'
        print '-a followed by assembly (eg hg19)'
        print '-s followed by species (eg hsa)'
        print '-n followed by annotation file (eg a refGene or ensGene file)'
        print '-t followed by annotation type (refGene, ensGene, mgcGenes and knownGene formats are supported)'
        print 'Optional arguments'
        print '-o followed by output directory (default = transcriptomeFastaFiles)'
        print '-j followed by an overlapJson file (only needed when jsonFile is already created, and does not have default name or location from transcriptomeFasta2)'
        
