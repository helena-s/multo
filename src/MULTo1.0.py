'''
Created on Apr 12, 2011

@author: helena_storvall
'''
"""
Query Kmers using bowtie
"""

import os, sys, subprocess, time
import threading, Queue
import optparse
import math
from array import array
import mappingConversionSimple
import json
#import make_bisulfiteGenome3
#import transcriptomeFasta2
import random
import fcntl
import shutil
import makeBisulfiteFasta
import makeTranscriptomeFasta


def reverse_seq(seq):
    rcMapDNA = {'A':'T', 'C': 'G', 'G': 'C', 'U': 'A', 'T': 'A', 'N': 'N',
            'a':'t', 'c': 'g', 'g': 'c', 'u': 'a', 't': 'a', 'n': 'n'}
    return ''.join([rcMapDNA.get(b, 'N') for b in reversed(seq)])

def get_seqDict(fa_f, cur_text, kMin, kMax, startPos, myChr, lastBlock):
    #faDict = {}#I write the reads to a file instead of making a dictionary
    no_reads = 0
    if lastBlock:
        lastPos = len(cur_text)-kMin+1
    else:
        lastPos = len(cur_text)-kMax+1
    for i in range(0, lastPos):
        pos = str(startPos + i)
        try:
            seq = cur_text[i:i+kMax]
        except:
            seq = cur_text[i:]
        if seq.count('N') > 0: 
            seq = seq[:seq.find('N')]
        readName = myChr+'.'+pos
        if not seq == '':#No point in printing a readname of a non-existing sequence (non-existing due to Ns in the sequence)
            fa_f.write('%s\t%s\n'%(readName, seq))
        #still count every position that can contain a read, even if the read could not be created due to Ns.
        no_reads += 1
    #print '\tDone making fasta dictionary'+' '+time.asctime() 
    return no_reads

def get_pairedSeqDict(fa_f, cur_text, kMin, kMax, startPos, myChr, insertInfo):
    #faDict = {}
    insertLength, lengthStd, no_randoms = map(int, insertInfo.split(','))
    #l = insertLength #Use this if we want to create the exact same in
    no_reads = 0
    for i in range(0, len(cur_text)-insertLength+1):
        pos = str(startPos +i)
        r_idx = 0
        while r_idx < no_randoms:
            r_idx += 1
            l = int(round(random.gauss(insertLength, lengthStd), 0))
            #print 'Randomized length for transcript %s, insert no %i: '%(myChr, r_idx), l
            if i+l > len(cur_text):
                continue
            seq1 = cur_text[i:i+kMax]
            seq2 = reverse_seq(cur_text[i+l-kMax:i+l])
            if seq1.count('N') > 0: 
                npos = seq1.find('N')
                if npos < kMin: continue
                else: seq1 = seq1[:npos]
            if seq2.count('N') > 0:
                npos = seq2.find('N')
                if (kMax-1-npos) < kMin: continue
                else: seq2 = seq2[npos+1:]
            readName = '%s.%s.r%i'%(myChr, pos, r_idx)
            #faDict[readName] = [seq1, seq2]
            fa_f.write('%s\t%s\t%s\n'%(readName, seq1, seq2))
            no_reads += 1
    return no_reads, len(cur_text)

def files_to_dict(faFile, baFile):
    fa_f = open(faFile)
    ba_f = open(baFile)
    faDict = {}
    binArrayDict = {}
    if PAIREDEND:
        for line in fa_f:
            readName, seq1, seq2 = line.strip('\n').split('\t')
            faDict[readName] = [seq1, seq2]
    else:
        for line in fa_f:
            readName, seq = line.strip('\n').split('\t')
            faDict[readName] = seq
    fa_f.close()
    if COUNTONLY:
        for line in ba_f:
            myTrans, AI, length, kMax = line.strip('\n').split('\t')
            binArrayDict[myTrans] = [int(AI), int(length), [0 for i in range(0,int(kMax))]]
    else:
        for line in ba_f:
            myTrans, transLength = line.strip('\n').split('\t')
            binArray = array('B', [0 for x in range(int(transLength))])
            binArrayDict[myTrans] = binArray
    return faDict, binArrayDict


def split_dict(faDict, bowtieList):
    #print '\t\tSplitting dictionary into unique and non-unique'+' '+time.asctime() 
    inListDict = {}
    restDict = faDict
    for pos in bowtieList:
        inListDict[pos] = restDict.get(pos, None)#restDict[pos]
        if inListDict[pos] == None:
            print 'WARNING! No such id in faDict: %s.\nTerminating program.', pos
            #print 'faDict', faDict
            sys.exit()
        del restDict[pos]
    return inListDict, restDict


def mateMapCheck(bowtieList):
    newBowtieList = []

    bowtieDict = {}
    for read in bowtieList:
        readName, mate = read.split('/')
        bowtieDict.setdefault(readName, []).append(mate)
    if MULTIMAP:
        for readName, mateList in bowtieDict.items():
            if len(mateList) == 2:
                newBowtieList.append(readName)
    else:       #For genomic queries it is unnescesary to check if both mates are in the list, since bowtie will only count mappings as valid if both ends map.
                #Observe though, for Bowtie2 mates will be mapped individually if they are not mapped as a pair, unless this option is disabled!!! 
        newBowtieList = bowtieDict.keys()
    return newBowtieList
        

def write_paired_files(faDict, k,  tmpdir):
    newFasta1 = os.path.join(tmpdir, 'tmp_%s_%i_r1.fa'%(NAME,k))
    newFasta2 = os.path.join(tmpdir, 'tmp_%s_%i_r2.fa'%(NAME,k))
    nf1 = open(newFasta1, 'w')
    nf2 = open(newFasta2, 'w')
    shortReads = []
    for key, val in faDict.iteritems():
        r1, r2 = val
        if len(r1) < k and len(r2) < k:
            shortReads.append(key)
        else:
            seq1 = r1[:min(k, len(r1))]
            seq2 = r2[:min(k, len(r2))]
            nf1.write(">%s\n%s\n" % (key, seq1))
            nf2.write(">%s\n%s\n" % (key, seq2))  
    nf1.close()
    nf2.close()
    newFasta = '%s,%s'%(newFasta1,newFasta2)
    return newFasta, shortReads

def write_fasta(faDict, k, tmpdir):
    #print '\t\tMaking fasta file with readlengths %i bp'%k +' '+time.asctime() 
    if PAIREDEND:
        return write_paired_files(faDict, k, tmpdir)
    newFasta = os.path.join(tmpdir, 'tmp_%s_%i.fa'%(NAME,k))
    nf = open(newFasta, 'w')
    shortReads = []
    for key, val in faDict.iteritems():
        if len(val) >= k:
            seq = val[:k]
            nf.write(">%s\n%s\n" % (key, seq))
        else:
            shortReads.append(key)
    nf.close()
    return newFasta, shortReads
        
def write_toBinArray(binArrayDict, outDict, k, startPos):
    #print '\t\tWriting %i bp unique to array'%k +' '+time.asctime()
    for key in outDict.iterkeys():
        myChr, pos = key.split('.')[:2]
        binArray = binArrayDict[myChr]
        binArray[int(pos)-startPos] = k
    del outDict
    #print '\t\tdone'+' '+time.asctime()

def add_count(binArrayDict, outDict, k, startPos):
    for key in outDict.iterkeys():
        #myChr, pos, r_idx = key.split('.')
        myChr = key.split('.')[0]
        binArrayDict[myChr][2][k-1] += 1
    del outDict
    
def array_toFile(binArrayDict, kMin, kMax, outdir):          
    for myChr, binArray in binArrayDict.items():
        binFile = os.path.join(outdir, '%s.unique%i-%i.btxt'%(myChr,kMin,kMax))
        if os.path.isfile(binFile):
            if OVERRIDE:
                print 'Regenerating the file %s'%binFile
                subprocess.call('rm %s'%binFile, shell = True)
            else:
                print 'The file %s has already been generated. Chose override option if you wish to re-generate it.'%binFile
                sys.exit()
        f=open(binFile,'ab')
        binArray.tofile(f)
        f.close()

def count_toFile(binArrayDict, kMin, kMax, outdir):
    outFile = os.path.join(outdir, '%s.uniqueLengths%i-%i.counts.txt'%(NAME,kMin,kMax))
    outFile2 = os.path.join(outdir, '%s.uniqueLengths%i-%i.readCounts.txt'%(NAME,kMin,kMax))
    if os.path.isfile(outFile):
        if OVERRIDE:
            print 'Regenerating the file %s'%outFile
            subprocess.call('rm %s'%outFile, shell = True)
        else:
            print 'The file %s has already been generated. Chose override option if you wish to re-generate it.'%outFile
            sys.exit()
    f=open(outFile,'w')
    f2 = open(outFile2, 'w')
    for myChr, binArray in binArrayDict.items():
        AI, length, propList = binArray
        if int(AI) == 0: continue
        UIk = 0
        allProp = []
        allCounts = []
        for i in range(kMin-1, kMax):#p in propList:
            p = propList[i]
            UIk +=p
            ULk = round(float(UIk*length)/float(AI), 1)
            allCounts.append(str(UIk))
            allProp.append(str(ULk))
        f.write('%s\t%s\t%i\t%i\n'%(myChr, '\t'.join(allProp), length, AI))
    f.close()

def run_bowtie(faFile, bowtie_idx, tmpdir, logf, insertInfo, multiMap = '-m 1', sup = '2,3,4,5,6,7,8', missmatch = '0', st = ''):
    #print '\t\tMapping %s with bowtie'%faFile+' '+time.asctime()
    file = open(LOCKFILE, 'a')
    fcntl.flock(file.fileno(), fcntl.LOCK_EX)
    if PAIREDEND:
        firstFile, secondFile = faFile.split(',')
        filePrefix = firstFile.split('/')[-1].strip('.fa')
        bowtieFile = os.path.join(tmpdir, '%s.bowtie%s'%(filePrefix,st))
        insertLength, lengthStd, no_randoms = map(int, insertInfo.split(','))
        lowLim = insertLength - lengthStd*3
        highLim = insertLength + lengthStd*3
        cmd = 'bowtie -v %s %s -p %i -f -I %i -X %i --suppress %s %s -1 %s -2 %s > %s'% (missmatch, multiMap, BOWPROC, lowLim, highLim, sup, bowtie_idx, firstFile, secondFile, bowtieFile)
    else:
        filePrefix = faFile.split('/')[-1].strip('.fa')
        bowtieFile = os.path.join(tmpdir, '%s.bowtie%s'%(filePrefix,st))
        cmd = 'bowtie -v %s %s -p %i -f --suppress %s %s %s > %s'% (missmatch, multiMap, BOWPROC, sup, bowtie_idx, faFile, bowtieFile)
    file.write('%s\tStarted command: %s\n'%(time.asctime(), cmd))    
    subprocess.call(cmd , shell=True, stderr=logf)
    file.write('%s\tDone with command: %s\n\n'%(time.asctime(), cmd)) 
    file.close()# unlocks the file, next thread can start run_bowtie
    #print cmd % (missmatch, multiMap, sup, bowtie_idx, faFile, bowtieFile)
    #print 'finished mapping with bowtie'+' '+time.asctime()
    #print '\t\tdone'+' '+time.asctime()
    return bowtieFile

def bowtie_fileToList(bowtieFile):
    #print '\t\tMaking bowtieFile %s into list'%bowtieFile +' '+time.asctime()
    bf = open(bowtieFile)
    bl = []
    for line in bf:
        bPos = line.strip('\n').split('\t')[0]
        bl.append(bPos)
    if not KEEP:
        subprocess.call('rm %s'%bowtieFile, shell=True)
    #print '\t\tdone making bowtieFile'+' '+time.asctime()
    return bl

def remove_genomicHits(bowtieFile, k):
    genBowtieFile = bowtieFile+'.gen'
    transBowtieFile = bowtieFile+'.trans'
    subprocess.call('grep chr %s > %s'%(bowtieFile,genBowtieFile), shell = True )
    subprocess.call('grep -v chr %s > %s'%(bowtieFile,transBowtieFile), shell = True )
    if not KEEP:
        subprocess.call('rm %s'%bowtieFile, shell=True)
    convGenBowtieFile = genBowtieFile+'.conv'
    convTransBowtieFile = transBowtieFile+'.conv'
    TX.convert_transcriptomefile(genBowtieFile, convGenBowtieFile, k)
    TX.convert_transcriptomefile(transBowtieFile, convTransBowtieFile, k)
    if not KEEP:
        subprocess.call('rm %s'%genBowtieFile, shell=True)
        subprocess.call('rm %s'%transBowtieFile, shell=True)
    bl = []
    
    tf = open(convTransBowtieFile)
    tfDict = {}
    for line in tf:
        tabs = line.strip('\n').split('\t')
        tfDict.setdefault(tabs[0], []).append(tabs[1:])
    if not KEEP:
        subprocess.call('rm %s'%convTransBowtieFile, shell=True)
    
    gf = open(convGenBowtieFile)
    gfDict = {}
    for line in gf:
        tabs = line.strip('\n').split('\t')
        gfDict.setdefault(tabs[0], []).append(tabs[1:])
    if not KEEP:
        subprocess.call('rm %s'%convGenBowtieFile, shell=True)
    
    for key, val in tfDict.items():
        if len(val) > 1:continue#If the read is mapped to more than one place in the transcriptome, it is not unique
        else:
            gVal = gfDict.get(key, None)
            if gVal == None:
                bl.append(key)
                continue
            if len(gVal) > 1: continue#If there is more than one genomic mapping, read is not unique                
            else:
                if not val[0][1] == gVal[0][1]: continue #If chromosomes are different
                elif not int(val[0][2]) == int(gVal[0][2]) or not int(val[0][3]) == int(gVal[0][3]): continue #If start or end position is different (between genomic and transcriptomic hit)
                else: 
                    bl.append(key)
    return bl
            

def remove_multimapping(bowtieFile, k):
    #print 'Removing multimapping from %s'%bowtieFile
    convBowtieFile = bowtieFile+'.genomic'
    TX.convert_transcriptomefile(bowtieFile, convBowtieFile, k)
    if not KEEP:
        subprocess.call('rm %s'%bowtieFile, shell=True)
    bf = open(convBowtieFile)
    bowtieDict = {}
    for line in bf:
        tabs = line.strip('\n').split('\t')#[readName, strand, chrom, start, end]
        bowtieDict.setdefault(tabs[0], {}).setdefault(tabs[2], [{},{}])
        bowtieDict[tabs[0]][tabs[2]][0].setdefault(int(tabs[3]), 0)
        bowtieDict[tabs[0]][tabs[2]][1].setdefault(int(tabs[4]), 0)
    if not KEEP:
        subprocess.call('rm %s'%convBowtieFile, shell=True)
    bl = []
    for readName, chrDict in bowtieDict.items():
        if len(chrDict) > 1: continue#If the read maps to several chromosome, it can not be unique
        for chrom, val in chrDict.items():
            startDict = val[0]
            endDict = val[1]
            if len(startDict) == 1 or len(endDict) == 1:
                #If there is only one startPosition and/or one endPosition, 
                #the read maps to overlapping places in different isoforms.
                #This is counted as unique.
                bl.append(readName)
    return bl        

def search_kmere(faDict, k, startPos, bowtie_idx, tmpdir, multiMap, sup, missmatch, logf, insertInfo):
    #print '\tSearch kmere of length %i' %k+' '+time.asctime()
    '''
    print '\t\tSearching %s\tBlock: %i\tk: %i\tlen(faDict) = %i\t%s'%(NAME, blockid, k,len(faDict), time.asctime())
    '''
    if len(faDict) > 0:   
        faFile, shortReads = write_fasta(faDict, k, tmpdir)
        if BISULFITE:
            forward_idx, reverse_idx = bowtie_idx.split(',')
            forwardFile = run_bowtie(faFile, forward_idx, tmpdir, logf, insertInfo, multiMap, sup, missmatch, '.fw')
            reverseFile = run_bowtie(faFile, reverse_idx, tmpdir, logf, insertInfo, '-k 1 --norc', sup, missmatch, '.rc')
            #Remove all reads in reverseFile from forwardFile because they are multimapping
            #If a read exsist only in reverseFile, it means it was already multimapping in the forward strand, and can be ignored
            bowtieFile = forwardFile.replace('.fw', '')
            subprocess.call('grep -Fxvf %s %s > %s'%(reverseFile, forwardFile, bowtieFile), shell = True)
            if not KEEP:
                subprocess.call('rm %s'%forwardFile, shell=True)
                subprocess.call('rm %s'%reverseFile, shell=True)
            bowtieList = bowtie_fileToList(bowtieFile)
        else: 
            bowtieFile = run_bowtie(faFile, bowtie_idx, tmpdir, logf, insertInfo, multiMap, sup, missmatch)
            if MULTIMAP:
                if GENE:
                    bowtieList = remove_multimapping(bowtieFile, k)
                elif TRANS:
                    bowtieList = remove_genomicHits(bowtieFile, k)
            else:
                bowtieList = bowtie_fileToList(bowtieFile)
        if not KEEP:
            if PAIREDEND: 
                fa1,fa2 = faFile.split(',')
                subprocess.call('rm %s'%fa1, shell=True)
                subprocess.call('rm %s'%fa2, shell=True)
            else:
                subprocess.call('rm %s'%faFile, shell=True) #Remove the fasta file when it is mapped with bowtie
        if PAIREDEND:
            bowtieList = mateMapCheck(bowtieList)
        uniqueDict, faDict = split_dict(faDict, bowtieList)
        shortDict, faDict = split_dict(faDict, shortReads)
        return uniqueDict, faDict, shortDict
    else:
        #print 'The dictionary for kmere of %i bp in block %i was empty' %(k, blockid)
        return {}, {}, {}

def multiple_search_kmere(faDict, kLow, kHigh, binArrayDict, startPos, bowtie_idx, tmpdir, multiMap, sup, missmatch, logf, insertInfo):
    if len(faDict) == 0:
        return
    for k in range(kLow, kHigh):
        uniqueDict, faDict, shortDict = search_kmere(faDict, k, startPos, bowtie_idx, tmpdir, multiMap, sup, missmatch, logf, insertInfo)
        del shortDict
        if COUNTONLY:
            add_count(binArrayDict, uniqueDict, k, startPos)
        else:
            write_toBinArray(binArrayDict, uniqueDict, k, startPos)
     #We already know they are unique at the kHigh, so if they are not at kHigh-1 they must be at kHigh.
     #No need to run bowtie again, just write to array.
    restDict = {}
    for key, val in faDict.iteritems():
        if len(val) >= kHigh:
            restDict[key] = val
    del faDict
    if COUNTONLY:
        add_count(binArrayDict, restDict, kHigh, startPos)
    else:
        write_toBinArray(binArrayDict, restDict, kHigh, startPos)

def merge_dict(dictA, dictB):
    dictC = dict(dictA.items() + dictB.items())
    return dictC

def stepwise_kmereSearch(kMin, kMax, faDict, binArrayDict, startPos, bowtie_idx, tmpdir, multiMap, sup, missmatch, logf, insertInfo):
    jumpStep = 25
    if len(faDict) == 0:
        return    
    #print '\tQuerying max length %i bp'%kMax +' '+time.asctime()
    faDict, restDict, shortDict = search_kmere(faDict, kMax, startPos, bowtie_idx, tmpdir, multiMap, sup, missmatch, logf, insertInfo)
    faDict = merge_dict(faDict, shortDict)
    del restDict #We do not need to write those who are not unique at kMax to the array, they should remain 0
    kLow = kMin
    #print '\tQuerying min length: %i bp'%kLow+' '+time.asctime()#25
    uniqueDict, faDict, shortDict = search_kmere(faDict, kLow, startPos, bowtie_idx, tmpdir, multiMap, sup, missmatch, logf, insertInfo)
    if len(shortDict) > 0:
        print 'Warning! ShortDict should never contain items after search with minLength'
    if COUNTONLY:
        add_count(binArrayDict, uniqueDict, kLow, startPos)
    else:
        write_toBinArray(binArrayDict, uniqueDict, kLow, startPos)    
    kHigh = kMin+jumpStep
    while kHigh < kMax:
        #print '\tQuerying: %i bp'%kHigh+' '+time.asctime()#50
        uniqueDict, faDict, shortDict = search_kmere(faDict, kHigh, startPos, bowtie_idx, tmpdir, multiMap, sup, missmatch, logf, insertInfo)
        uniqueDict = merge_dict(uniqueDict, shortDict)
        #print '\tQuerying: %i-%i bp'%(kLow, kHigh)+' '+time.asctime()#25-50
        multiple_search_kmere(uniqueDict, kLow+1, kHigh, binArrayDict, startPos, bowtie_idx, tmpdir, multiMap, sup, missmatch, logf, insertInfo)
        kLow = kHigh
        kHigh += jumpStep
    #print '\tQuerying last block: %i-%i bp'%(kHigh, kMax)+' '+time.asctime()#75-100
    multiple_search_kmere(faDict, kLow+1, kMax, binArrayDict, startPos, bowtie_idx, tmpdir, multiMap, sup, missmatch, logf, insertInfo)


def process_transcripts(inInfo, kMin, kMax, tmpdir, BLOCKSIZE, PAIREDEND, COUNTONLY, insertInfo):
    trans_dir,isoName = inInfo.split(',')
    print '\tprocessing %s... %s' %(trans_dir, time.asctime())
    transFiles = [f for f in os.listdir(trans_dir) if f.endswith(".fa")]
    blockid = 1
    startPos = 0 #Always stays 0, we never interrupt in the middle of a transcript
    blockReads = 0
    blockTmpDir = os.path.join(tmpdir, 'block_%i_%s'%(blockid, isoName))
    if not os.path.isdir(blockTmpDir): os.mkdir(blockTmpDir)
    faFile = '%s/block_%i_%s_faDict.txt'%(blockTmpDir, blockid, isoName)
    baFile = '%s/block_%i_%s_binArrayDict.txt'%(blockTmpDir, blockid, isoName)
    fa_f = open(faFile, 'w')
    ba_f = open(baFile, 'w')
    print '\t\tMaking fasta dictionary for %s, block %i'%(trans_dir, blockid)+' '+time.asctime() 
    for trans_f in transFiles:
        tf = open(os.path.join(trans_dir, trans_f), 'r')
        myTrans = trans_f.split('.')[0]
        header = tf.readline()
        cur_text = tf.read().replace('\n','').upper()

        if PAIREDEND:
            insertLength, lengthStd, no_randoms = map(int, insertInfo.split(','))
            if len(cur_text) < insertLength-lengthStd*3: continue #Skip all short transcripts
            no_reads, length = get_pairedSeqDict(fa_f, cur_text, kMin, kMax, startPos, myTrans, insertInfo)
        else:
            no_reads = get_seqDict(fa_f, cur_text, kMin, kMax, startPos, myTrans, True)
            length = len(cur_text)
            
        if COUNTONLY: 
            ba_f.write('%s\t%i\t%i\t%i\n'%(myTrans, no_reads, length, kMax))
        else:
            ba_f.write('%s\t%i\n'%(myTrans, len(cur_text)))
            
        blockReads += no_reads
        if blockReads >= BLOCKSIZE:
            print '\t\tDone making fasta dictionary for %s, block %i'%(trans_dir, blockid)+' '+time.asctime()
            fa_f.close()
            ba_f.close()
            blockid += 1
            blockReads = 0
            blockTmpDir = os.path.join(tmpdir, 'block_%i_%s'%(blockid, isoName))
            if not os.path.isdir(blockTmpDir): os.mkdir(blockTmpDir)
            faFile = '%s/block_%i_%s_faDict.txt'%(blockTmpDir, blockid, isoName)
            baFile = '%s/block_%i_%s_binArrayDict.txt'%(blockTmpDir, blockid, isoName)
            fa_f = open(faFile, 'w')
            ba_f = open(baFile, 'w')
            print '\t\tMaking fasta dictionary for %s, block %i'%(trans_dir, blockid)+' '+time.asctime() 
            
    if blockReads >= 0:     
        fa_f.close()
        ba_f.close()
        print '\t\tDone making fasta dictionary for %s, block %i'%(trans_dir, blockid)+' '+time.asctime()

          
def process_chromosomes(inInfo, kMin, kMax, tmpdir, BLOCKSIZE):
    
    #chrFiles = [f for f in os.listdir(chrom_dir) if f.endswith(".fa")]
    print 'This is the tmpdir: %s'%tmpdir
    chr_f, chrName = inInfo.split(',')
    print '\tprocessing %s... %s' %(chr_f, time.asctime())
    in_f = chr_f.split('/')[-1]
    myChr = in_f.split('.')[0]   
    blockid = 1
    startPos = 0 #Always stays 0, we never interrupt in the middle of a transcript
    blockReads = 0
    blockTmpDir = os.path.join(tmpdir, 'block_%i_%s'%(blockid, chrName))
    if not os.path.isdir(blockTmpDir): os.mkdir(blockTmpDir)
    faFile = '%s/block_%i_%s_faDict.txt'%(blockTmpDir, blockid, chrName)
    baFile = '%s/block_%i_%s_binArrayDict.txt'%(blockTmpDir, blockid, chrName)
    fa_f = open(faFile, 'w')
    ba_f = open(baFile, 'w')
    chr_fh = open(chr_f, 'r')
    header = chr_fh.readline() # the >chr line
    cur_text = chr_fh.read(BLOCKSIZE)
    cur_text = cur_text.replace('\n','').upper()
    firstBlockSize = len(cur_text)
    print 'firstBlockSize', firstBlockSize
    while cur_text:
        print '\t\tMaking fasta dictionary for %s, block %i'%(chrName, blockid)+' '+time.asctime() 
        ba_f.write('%sb%i\t%i\n'%(chrName, blockid, len(cur_text)-kMax+1))
        no_reads = get_seqDict(fa_f, cur_text, kMin, kMax, 0, '%sb%i'%(chrName, blockid), False)#OBS! Changed from startPos to 0, because we now use position in block rather than chromosomal position
            
        startPos += len(cur_text) - (kMax-1)
        cur_text = cur_text[-(kMax-1):] + chr_fh.read(BLOCKSIZE)
        cur_text = cur_text.replace('\n','').upper()
        blockid += 1
        blockTmpDir = os.path.join(tmpdir, 'block_%i_%s'%(blockid, chrName))
        if not os.path.isdir(blockTmpDir): os.mkdir(blockTmpDir)
        faFile = '%s/block_%i_%s_faDict.txt'%(blockTmpDir, blockid, chrName)
        baFile = '%s/block_%i_%s_binArrayDict.txt'%(blockTmpDir, blockid, chrName)
        fa_f = open(faFile, 'w')
        ba_f = open(baFile, 'w')
        if len(cur_text) < firstBlockSize: 
            ba_f.write('%sb%i\t%i\n'%(chrName, blockid, len(cur_text)))
            no_reads = get_seqDict(fa_f, cur_text, kMin, kMax, 0, '%sb%i'%(chrName, blockid), True)
            break
    print '\tDone processing %s %s' %(chr_f, time.asctime())

def get_tmpDir(outdir, chr_f):
    t1dir = os.path.join(outdir, 'tmp_files/')
    if not os.path.isdir(t1dir):
        try:
            os.mkdir(t1dir)
        except:
            print 'Directory %s was already created by another thread'%t1dir
    logf = open(os.path.join(outdir, 'bowtieLog.txt'), 'w')
    return t1dir, logf
    
def query_chromosome(tmpdir, kMin, kMax, outdir, bowtie_idx, bowtie_proc, ass, override, bs, jsonFile, transLevel, geneLevel, keepFiles, bisulfite, smallChrom, missmatch, countonly, pairedEnd, insertInfo, resumeRun):
    """collect all kmers from chromosome into a dictionary. Map only unique ones using bowtie
    and save the set of globally unique ones."""
    
    global KEEP
    KEEP = keepFiles
    global OVERRIDE
    OVERRIDE = override
    global BLOCKSIZE
    BLOCKSIZE = bs
    
    global TRANS
    TRANS = transLevel
    global GENE
    GENE = geneLevel
    global MULTIMAP
    MULTIMAP = False
    if TRANS or GENE:
        MULTIMAP = True
    global NAME
    NAME = tmpdir.split('/')[-1]
    global BISULFITE
    BISULFITE = bisulfite
    global COUNTONLY
    COUNTONLY = countonly
    global PAIREDEND
    PAIREDEND = pairedEnd
    global BOWPROC
    BOWPROC = bowtie_proc
    global LOCKFILE
    LOCKFILE = os.path.join(outdir, 'bowtie_lockfile.txt')
    
    if not os.path.exists(LOCKFILE):
        # create the counter file if it doesn't exist
        file = open(LOCKFILE, "w")
        file.write("Lockfile for MULTo1.0\n")
        file.close()
    multoDir = os.path.join(outdir, 'MULTo_files')
    if not os.path.isdir(multoDir): 
        try: os.mkdir(multoDir)
        except: print 'Directory %s was already created by another thread'%multoDir

    
    if MULTIMAP:
        global TX
        TX = mappingConversionSimple.tx2genome(jsonFile, ass)
        try:
            isoforms = int(tmpdir.split('_')[-1])
        except:
            print '\t\tFor multimapping transcripts, number of isoforms must be stated last in folder name! Current foldername: ',tmpdir
            sys.exit()
        useM = int(isoforms) +1 #Add one for genomic hits.
        print '\t\tUsing the bowtie flag -m %i for the folder %s'%(useM, tmpdir) #OBS! Remove after debugging!
        multiMap = '-m %i -a'%useM
        sup = '5,6,7,8'
        
    elif BISULFITE:
        strand = tmpdir.split('_')[-1]#will be fw or rv
        if strand == 'rv': #Reverse the order of the bowtie-indexes
            fw_idx,rv_idx = bowtie_idx.split(',')
            bowtie_idx = rv_idx+','+fw_idx
        multoDir = os.path.join(multoDir, strand)
        if not os.path.isdir(multoDir): 
            try: os.mkdir(multoDir)
            except: print 'Directory %s was already created by another thread'%multoDir
        multiMap = '-m 1 --norc'
        sup = '2,3,4,5,6,7,8'
        #process_chromosomes(chr_f, kMin, kMax, outdir, bowtie_idx, tmpdir, logf, multiMap, sup)
    else:
        multiMap = '-m 1'
        sup = '2,3,4,5,6,7,8'
        #process_chromosomes(chr_f, kMin, kMax, outdir, bowtie_idx, tmpdir, logf, multiMap, sup)
    #resumeRun = True
    if resumeRun:
        print 'RESUME RUN MODE: checking... %s'%os.path.join(tmpdir, '%s.done'%NAME)
        if os.path.isfile(os.path.join(tmpdir, '%s.done'%NAME)):
            print 'RESUME RUN MODE: The block %s is fully processed. Skipping.'%NAME
            return
        
    faFile = os.path.join(tmpdir, '%s_faDict.txt'%NAME)
    baFile = os.path.join(tmpdir, '%s_binArrayDict.txt'%NAME)
    if os.path.exists(faFile) and os.path.exists(baFile):
        faDict, binArrayDict = files_to_dict(faFile, baFile)
    else:
        print 'WARNING! The block %s had no dictionary files!. Skipping.'%NAME
        return
    print 'MULTODIR: %s'%multoDir
    print '\t\tPerforming stepwise_kmereSearch on folder %s %s'%(tmpdir, time.asctime()) 
    logf = open(os.path.join(outdir, 'bowtieLog.txt'), 'w')
    stepwise_kmereSearch(kMin, kMax, faDict, binArrayDict, 0, bowtie_idx, tmpdir, multiMap, sup, missmatch, logf, insertInfo)
    print '\t\tWriting binArray to file %s'%time.asctime()
    if COUNTONLY:
        count_toFile(binArrayDict, kMin, kMax, multoDir)
    else:
        array_toFile(binArrayDict, kMin, kMax, multoDir)
    df = open(os.path.join(tmpdir, '%s.done'%NAME), 'w')
    df.write('Block finished at %s'%time.asctime())
    df.close()


def check_index(bowtie_idx):
    allSuffixes = ['.1.ebwt','.2.ebwt','.3.ebwt','.4.ebwt','.rev.1.ebwt','.rev.2.ebwt']
    filesExist = True
    for suf in allSuffixes:
        if not os.path.isfile(bowtie_idx+suf):
            print 'The bowtie-index file %s does not exist. Make sure that indexes has been generated'%(bowtie_idx+suf)
            filesExist = False
    if not filesExist:
        print 'Program needs bowtie-indexes to run. Terminating.'
        sys.exit()

def merge_binary_files(multoDir, keepFiles, kMin, kMax):
    binFiles = [f for f in os.listdir(multoDir) if f.endswith('btxt')]
    print binFiles
    blockDict = {}
    for f in binFiles:
        chromID = f.split('.')[0]
        try:
            myChr, blockID = chromID.split('b')[:2]
        except:
            print 'f: %s, chromID: %s'%(f, chromID)
            sys.exit()
        if blockID.isdigit():
            blockDict.setdefault(myChr, []).append([int(blockID), os.path.join(multoDir, f)])
    print blockDict
    for myChr, fileList in blockDict.items():
        outFile = os.path.join(multoDir, '%s.unique%i-%i.btxt'%(myChr,kMin,kMax))
        fileList.sort()
        print fileList
        fout = file(outFile, 'wb')
        for f in fileList:
            fin = file(f[1], 'rb')
            shutil.copyfileobj(fin, fout, 65536)
            fin.close()
            if not keepFiles:
                subprocess.call('rm %s'%f[1], shell=True)
        fout.close()
        
    
def query_all_chrs_threaded(genome_d, kMin, kMax, outdir, bowtie_idx, ass, bowtie_proc, override, bs, jsonFile, transLevel, geneLevel, keepFiles, bisulfite, bisRun, smallChrom, missmatch, countonly, pairedEnd, insertInfo, dictFilesExist):
    #print genome_d
    #print os.listdir(genome_d)
    if not os.path.isdir(outdir):
        create_allSubdirs(outdir)
    print outdir
    extra_opts = ''
    if override: extra_opts += '-O '
    if keepFiles: extra_opts += '-K '
    if smallChrom: extra_opts += '-S '
    if countonly: extra_opts += '-C '
    tmpdir = os.path.join(outdir, 'tmp_files/')
    logf = open(os.path.join(outdir, 'bowtieLog.txt'), 'w')
    if not os.path.isdir(tmpdir):
        try:
            os.mkdir(tmpdir)
        except:
            print 'Directory %s was already created by another thread'%tmpdir
    
    if pairedEnd: 
        extra_opts += '-P '
        if not len(insertInfo.split(',')) == 3:
            print 'WARNING! You need to specify insert length, standard deviation and number of randomized reads with -l flag for paired end reads.\nTerminating program.'
            sys.exit()
        extra_opts += '-l %s '%insertInfo
        
        
    if transLevel or geneLevel:
        #check_index(bowtie_idx)
        if bisulfite:
            print 'WARNING! Bisulfite queries are not supported at transLevel and geneLevel, only genomic.\nTerminating program.'
            sys.exit()
        print 'Creating transcriptome uniqueness files'
        
        jf = '-j %s '%jsonFile
        extra_opts += '-j %s '%jsonFile
        if transLevel: extra_opts += '-T '
        elif geneLevel: extra_opts += '-G '
        #Creating fasta dictionary files       

        if not dictFilesExist:
            infiles = [os.path.join(genome_d, f) for f in os.listdir(genome_d) if os.path.isdir(os.path.join(genome_d, f))]
            
            print '%s\tCreating dictionaries with fasta sequence, and writing them to files.'%time.asctime()
            for in_fp in infiles:
                isoName = in_fp.split('/')[-1]
                cmd_queue.put("python src/MULTo1.0.py -c %s -k %i -m %i -o %s -z %i %s" % (in_fp+','+isoName, kMin, kMax, tmpdir, bs, extra_opts))
            cmd_queue.join()
        
        
        
    elif bisulfite:        
        if transLevel or geneLevel:
            print 'WARNING! Bisulfite queries are not supported at transLevel and geneLevel, only genomic.\nTerminating program.'
            sys.exit()
        print 'Creating bisulfite genome uniqueness files'
        extra_opts += '-B '
        if not len(bowtie_idx.split(',')) == 2:
            print 'WARNING! Two bowtie indexes must be provided for bisulfite queries.(2)'
            print 'Terminating program.'
            sys.exit()
        else:
            fw_idx,rv_idx = bowtie_idx.split(',')
            check_index(fw_idx)
            check_index(rv_idx)
            print 'OBS! Make sure that forward index is entered first, and reverse index follows'
            print 'Forward index: %s ?'%fw_idx
            print 'Reverse index: %s ?'%rv_idx
        if not dictFilesExist:
            if bisRun not in [0,1,2]:
                print 'The choice of bisulfite run mode must be 0 (only forward), 1 (only reverse) or 2 (both)\nTerminating program'
                sys.exit()
            fwDir = os.path.join(genome_d, 'forward')
            rvDir = os.path.join(genome_d, 'reverse')
            fwFiles = [os.path.join(fwDir, f) for f in os.listdir(fwDir) if f.endswith(".fa")]
            rvFiles = [os.path.join(rvDir, f) for f in os.listdir(rvDir) if f.endswith(".fa")]
            if bisRun == 0 or bisRun == 2:
                for fw_f in fwFiles:
                    if fw_f.endswith('merged.fa'):continue
                    chrName = fw_f.split('/')[-1].split('.')[0]+'_fw'
                    cmd_queue.put("python src/MULTo1.0.py -c %s -k %i -m %i -o %s -z %i %s" % (fw_f+','+chrName, kMin, kMax, tmpdir, bs, extra_opts))
            if bisRun == 1 or bisRun == 2:
                for rv_f in rvFiles:
                    if rv_f.endswith('merged.fa'):continue
                    chrName = rv_f.split('/')[-1].split('.')[0]+'_rv'
                    cmd_queue.put("python src/MULTo1.0.py -c %s -k %i -m %i -o %s -z %i %s" % (rv_f+','+chrName, kMin, kMax, tmpdir, bs, extra_opts))   
            cmd_queue.join()
            
        
    elif smallChrom:
        extra_opts += '-S'
        check_index(bowtie_idx)
        if not dictFilesExist:
            infiles = [os.path.join(genome_d, f) for f in os.listdir(genome_d) if os.path.isdir(os.path.join(genome_d, f))]
            print '%s\tCreating dictionaries with fasta sequence, and writing them to files.'%time.asctime()
            for in_fp in infiles:
                isoName = in_fp.split('/')[-1]
                cmd_queue.put("python src/MULTo1.0.py -c %s -k %i -m %i -o %s -z %i %s" % (in_fp+','+isoName, kMin, kMax, tmpdir, bs, extra_opts))
            cmd_queue.join()
    else:
        print 'Initiating threads'
        check_index(bowtie_idx)
        if not dictFilesExist:
            infiles = [os.path.join(genome_d, f) for f in os.listdir(genome_d) if f.endswith(".fa")]
            for in_fp in infiles:
                if in_fp.endswith('merged.fa'):continue
                chrName = in_fp.split('/')[-1].split('.')[0]
                cmd_queue.put("python src/MULTo1.0.py -c %s -k %i -m %i -o %s -z %i %s" % (in_fp+','+chrName, kMin, kMax, tmpdir, bs, extra_opts))
            cmd_queue.join()
    
    #Performing uniqueness queries
    multoDir = os.path.join(outdir, 'MULTo_files')
    if not os.path.isdir(multoDir): 
        try: os.mkdir(multoDir)
        except: print 'Directory %s was already created by another thread'%multoDir
    print '%s\tProcessing blocks in parallell.'%time.asctime()
    allBlockDirs = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if os.path.isdir(os.path.join(tmpdir, f))]
    for blockDir in allBlockDirs:
        cmd_queue.put("python src/MULTo1.0.py -e %s -k %i -m %i -o %s -b %s -a %s -z %i -v %s -p %i %s" % (blockDir, kMin, kMax, outdir, bowtie_idx, ass, bs, missmatch, bowtie_proc, extra_opts))    
    cmd_queue.join()
    print '%s\tDone processing blocks in parallell.'%time.asctime()
    
    
    if countonly:
        print '%s\tMerging count files.'%time.asctime()
        mergedCountFile = '%s/all.uniqueLengths%i-%i.counts.txt'%(multoDir,kMin,kMax)
        f=open(mergedCountFile,'w')
        kList = ['%i'%a for a in range(kMin, kMax+1)]
        f.write('#transcriptName\t%s\tTotalLength\tSimulatedInserts\n'%'\t'.join(kList))
        f.close()
        subprocess.call('less %s/block*.uniqueLengths%i-%i.counts.txt >> %s'%(multoDir,kMin,kMax,mergedCountFile), shell = True)
        print '%s\tDone merging count files.'%time.asctime()
    
    if not transLevel or geneLevel or smallChrom:
        #print 'IMPLEMENT ME!!! THE BINARY FILES FROM ALL BLOCKS IN A CHROMOSOME MUST BE MERGED!!!'
        if bisulfite:
            merge_binary_files(os.path.join(multoDir, 'fw'), keepFiles, kMin, kMax)
            merge_binary_files(os.path.join(multoDir, 'rv'), keepFiles, kMin, kMax)
        else:
            merge_binary_files(multoDir, keepFiles, kMin, kMax)

def convert_annFile(annFile, annType, convAnnFile):
    print 'Converting the annotation file to make sure all transcript names are unique. %s'%time.asctime()
    if (annType.split('+')[0].lower() == 'refgene' or annType.split('+')[0].lower() == 'ensgene' or annType.split('+')[0].lower() == 'mgcgenes'):
        t = 1
    elif annType.split('+')[0].lower() == 'knowngene':
        t = 0
    af = open(annFile)
    nameDict = {}
    for line in af:
         tabs = line.strip('\n').split('\t')
         transName = tabs[t]
         chrom, strand = tabs[t+1:t+3]
         start, end = map(int, tabs[t+3:t+5])
         nameDict.setdefault(transName, []).append(tabs)
    af.close()
    caf = open(convAnnFile, 'w')
    for tName, transList in nameDict.iteritems():
        if len(transList) > 1:
            i = 0
            for trans in transList:
                i += 1
                transName = tName+'_c%i'%i
                trans[t] = transName
                caf.write('\t'.join(trans)+'\n')
        else:
            trans = transList[0]
            caf.write('\t'.join(trans)+'\n')
    caf.close()
    print 'Done converting annotation file. %s'%time.asctime()
    
def merge_genomic(genome_d, ass):
    genomeMerged = os.path.join(genome_d, '%s_merged.fa'%ass)
    if os.path.isfile(genomeMerged):
        subprocess.call('rm %s'%genomeMerged, shell = True)
    print 'Merging genomic fasta files to %s %s'%(genomeMerged, time.asctime())
    inFiles = [os.path.join(genome_d, f) for f in os.listdir(genome_d) if f.endswith(".fa")]
    mf = open(genomeMerged, 'w')
    for f in inFiles:
        tf = open(f)
        for line in tf:
            mf.write(line)
        tf.close()
    mf.close()
    print 'Done merging genomic fasta files to %s'%time.asctime()
    return genomeMerged
    
    
def build_index(bowtie_idx, genome_d, mergedFasta, bisulfite, transcriptome, ass):
    if bisulfite:
        print 'Creating bowtie-indexes for bisulfite genome forward and reverse %s'%time.asctime()
        if not len(bowtie_idx.split(',')) == 2:
            print 'WARNING! Names of two comma separated bowtie indexes (forward and reverse) must be provided to create bisulfite bowtie-index.\nTerminating program.'
            sys.exit()
        elif not len(mergedFasta.split(',')) == 2:
            print 'WARNING! Two comma separated fasta-files (one forward and one reverse) must be provided to create bisulfite bowtie-index.\nTerminating program.'
            sys.exit()
        forward_idx, reverse_idx = bowtie_idx.split(',')
        forwardFasta, reverseFasta = mergedFasta.split(',')
        print 'OBS! Make sure that forward fa-file and index-file are entered first, and reverse index follows'
        print 'Forward index: %s , Forward fasta-file: %s'%(forward_idx, forwardFasta)
        print 'Reverse index: %s , Reverse fasta-file: %s'%(reverse_idx, reverseFasta)
        if not os.path.isfile(forwardFasta):
            print 'WARNING! The forwardFasta-file "%s" does not exist!\nTerminating program.'%forwardFasta
            sys.exit()
            #merge_files(forwardFasta, os.path.join(genome_d, 'forward'), False)
        if not os.path.isfile(reverseFasta):
            print 'WARNING! The reverseFasta-file "%s" does not exist!\nTerminating program.'%reverseFasta
            sys.exit()
            #merge_files(reverseFasta, os.path.join(genome_d, 'reverse'), False)
        cmd_queue.put('bowtie-build %s %s'%(forwardFasta, forward_idx))
        cmd_queue.put('bowtie-build %s %s'%(reverseFasta, reverse_idx))
        
    elif transcriptome:
        print 'Creating bowtie-indexes for transcriptome+genome %s'%time.asctime()
        if not len(mergedFasta.split(',')) == 2:
            print 'WARNING! Two comma separated fasta-files (first genome, second transcriptome) must be provided to create transcriptome bowtie-index.'
            print 'Terminating program.'
            sys.exit()
        genomeFasta, transcriptomeFasta = mergedFasta.split(',')
        if not os.path.isfile(genomeFasta):
            print 'WARNING! The genomeFasta-file "%s" does not exist!\nTerminating program.'%genomeFasta
            sys.exit()
        if not os.path.isfile(genomeFasta):
            print 'WARNING! The transcriptomeFasta-file "%s" does not exist!\nTerminating program.'%transcriptomeFasta
            sys.exit()
        cmd_queue.put('bowtie-build %s,%s %s'%(genomeFasta, transcriptomeFasta, bowtie_idx))
    else:
        print 'Creating bowtie-indexes for genome %s'%time.asctime()
        if not os.path.isfile(mergedFasta):
            mergedFasta = merge_genomic(genome_d, ass)
        cmd_queue.put('bowtie-build %s %s'%(mergedFasta, bowtie_idx))
    cmd_queue.join()
    print 'Done creating bowtie-indexes %s'%time.asctime()
        #subprocess.call('bowtie-build %s %s'%(forwardFasta, forward_idx), shell = True)
    #bowtie-build transcriptFasta/mm9/refGene/all_mm9_refGene.fa,/mnt/crick/sandberglab/bowtieindex/mm9/mm9_norandom_0.12.5.fa /mnt/crick/sandberglab/bowtieindex/mm9_refGene/mm9_norandom_refGene

def create_allSubdirs(wholeDir):
    subDirs = wholeDir.split('/')
    for i in range(1,len(subDirs)+1):
        currDir = '/'.join(subDirs[0:i])
        if not os.path.isdir(currDir):
            os.mkdir(currDir)


def download_chromosomeFasta(spec, ass):
    chromDir = 'files/%s/%s/fastaFiles/genomeFasta/noRandomChrom/'%(spec, ass)
    genomeFastaDir ='files/%s/%s/fastaFiles/genomeFasta/'%(spec, ass)
    create_allSubdirs(chromDir)
    subprocess.call('wget --timestamping -P %s ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/chromosomes/*'%(chromDir, ass), shell = True)
    randomDir = 'files/%s/%s/fastaFiles/genomeFasta/randomChrom/'%(spec, ass)
    if not os.path.isdir(randomDir):
            os.mkdir(randomDir)
    for f in os.listdir(chromDir):
        #if len(f.split('.')[0]) > 5:
        if f.count('_') > 0:
            subprocess.call('mv %s %s'%(os.path.join(chromDir, f), randomDir), shell = True)
        if not f.endswith('gz'):
            subprocess.call('mv %s %s'%(os.path.join(chromDir, f), genomeFastaDir), shell = True)
    for f in os.listdir(chromDir):
        print 'To gzip: %s'%os.path.join(chromDir, f)
        subprocess.call('gzip -d %s'%os.path.join(chromDir, f), shell = True)
        
if __name__=='__main__':
    opts = optparse.OptionParser()
    opts.formatter=optparse.TitledHelpFormatter()
    opts.add_option('-g','--genome-dir',dest='genome_d')
    opts.add_option('-c', '--chromosome', dest= 'chr_f')
    opts.add_option('-e', '--blockDir', dest= 'blockDir')
    opts.add_option('-b','--bowtie-idx',dest='bowtie_idx')
    opts.add_option('-s','--species',dest='spec')
    opts.add_option('-a','--assembly', dest='ass')
    opts.add_option('-o','--out-dir', dest='outdir')
    opts.add_option('-v','--missmatches', dest='missmatches', default = '0')
    opts.add_option('-p','--processors', dest='multi', default = '25')
    opts.add_option('-k','--low-kmer', dest='low', default = '20')
    opts.add_option('-m', '--high-kmer', dest= 'high', default = '255')
    opts.add_option('-l', '--insert-info', dest= 'insertInfo', default = '250,25,5')#Give insert length, standard deviation and number of randomizations
    opts.add_option('-z', '--blocksize', dest = 'bs', default = 10000000)
    opts.add_option('-O', '--override', dest = 'override', default = False, action = 'store_true')
    opts.add_option('-S', '--small-chrom', dest = 'smallChrom', default = False, action = 'store_true')
    opts.add_option('-C', '--count-only', dest = 'countonly', default = False, action = 'store_true')
    opts.add_option('-P', '--paired-end', dest = 'pairedEnd', default = False, action = 'store_true')
    opts.add_option('-D', '--dict-files-exist', dest = 'dictFilesExist', default = False, action = 'store_true')
    #Following inputs are used for bisulfite genome runs
    opts.add_option('-r','--bisulfite-run-mode',dest='bisRun', default = '2')
    opts.add_option('-B','--bisulfite',dest='bisulfite', default = False, action = 'store_true')
    #Following inputs are used for transcriptome uniqueness runs
    opts.add_option('-j','--json-file', dest='jsonFile')
    opts.add_option('-T', '--transcript-level', dest = 'transLevel', default = False, action = 'store_true')
    opts.add_option('-G', '--gene-level', dest = 'geneLevel', default = False, action = 'store_true')
    #Following inputs are used if you want to create fasta-files and/or bowtie-indexes using this script
    opts.add_option('-I','--create-index',dest='createIndex', default = False, action = 'store_true')
    opts.add_option('-F','--create-fasta',dest='createFasta', default = False, action = 'store_true')
    opts.add_option('-J','--create-json',dest='createJson', default = False, action = 'store_true')
    opts.add_option('-A','--create-all',dest='createAll', default = False, action = 'store_true')
    
    #Following inputs are used when creating new transcriptomeFasta-files
    opts.add_option('-n','--annotation-file',dest='annFile')
    opts.add_option('-t','--annotation-type', dest='annType')
    #keep-files should only be used for debugging
    opts.add_option('-K','--keep-files',dest='keepFiles', default = False, action = 'store_true')
    opts.add_option('-R','--resume-run',dest='resumeRun', default = False, action = 'store_true')
    opts.add_option('--print-help', dest = 'help', default = False, action = 'store_true')
    
    
    (o, args) = opts.parse_args()
    #print 'Running MULTo1.0, %s'%time.asctime()
    print ' '.join(sys.argv)
    bowtie_proc = min(int(o.multi)-1, 15)
    #no_threads = max(1, o.multi-bowtie_proc)
    
    if not o.chr_f is None:
        if o.transLevel or o.geneLevel or o.smallChrom:
            process_transcripts(o.chr_f, int(o.low), int(o.high), o.outdir, int(o.bs), o.pairedEnd, o.countonly, o.insertInfo)
        else:
            process_chromosomes(o.chr_f, int(o.low), int(o.high), o.outdir, int(o.bs))
    elif not o.blockDir is None:
        if o.bisulfite:
            if not len(o.bowtie_idx.split(',')) == 2:
                print 'WARNING! Two bowtie indexes must be provided for bisulfite queries. (1)'
                print 'Terminating program.'
                sys.exit()
            else:
                #If reverse index is entered first the results will be faulty
                print 'OBS! Make sure that forward index is entered first, and reverse index follows.'
                print 'Forward index: %s ?'%o.bowtie_idx.split(',')[0]
                print 'Reverse index: %s ?'%o.bowtie_idx.split(',')[1]
        query_chromosome(o.blockDir, int(o.low), int(o.high), o.outdir, o.bowtie_idx, int(o.multi), o.ass, o.override, int(o.bs), o.jsonFile, o.transLevel, o.geneLevel, o.keepFiles, o.bisulfite, o.smallChrom, o.missmatches, o.countonly, o.pairedEnd, o.insertInfo, o.resumeRun)
        #print 'Finished with file: %s '%o.chr_f+time.asctime()
    elif not o.ass is None and not o.spec is None:
        print 'Running unique_range with multiple chromosomes'+time.asctime()
        # Setup multi-threading support
        jsonFile = None
        if o.genome_d is None:
            genomeDir = 'files/%s/%s/fastaFiles/genomeFasta/noRandomChrom/'%(o.spec, o.ass)
        else:
            genomeDir = o.genome_d
        if not os.path.isdir(genomeDir):
            download_chromosomeFasta(o.spec, o.ass)
        genomeMerged = os.path.join(genomeDir, '%s_genomeMerge.fa'%o.ass)
        if not os.path.isfile(genomeMerged):
            genomeMerged = merge_genomic(genomeDir, o.ass)
            
        cmd_queue = Queue.Queue()
        num_threads = int(o.multi)-bowtie_proc
        def process_command(i, q):
            while True:
                cmd = q.get()
                print "[%i] %s" % (i, cmd)
                subprocess.call(cmd, shell=True)
                q.task_done()                    
        # set up some threads
        for i in range(num_threads):
            print 'starting thread %i'%i+time.asctime()
            worker = threading.Thread(target=process_command, args=(i, cmd_queue,))
            worker.setDaemon(True)
            worker.start()
        
        #BISULFITE MULTo RUN
        if o.bisulfite: 
            print 'Running in bisulfite mode'
            outDir = 'files/%s/%s/MULfiles/%s_bisulfite_%s-%s'%(o.spec, o.ass, o.ass, o.low, o.high)
            bisulfiteDir = 'files/%s/%s/fastaFiles/bisulfiteFasta/'%(o.spec, o.ass)
            inFastaDir = bisulfiteDir
            fwDir = os.path.join(bisulfiteDir, 'forward')
            reDir = os.path.join(bisulfiteDir, 'reverse')
            fwMerged = os.path.join(bisulfiteDir, '%s_forward_bisulfiteMerge.fa'%(o.ass))
            reMerged = os.path.join(bisulfiteDir, '%s_reverse_bisulfiteMerge.fa'%(o.ass))
            #if o.bisRun == 2 or o.bisRun == 0 and (not os.path.isdir(fwDir) or len(os.listdir(fwDir)) == 0):
            if not os.path.isdir(fwDir) or len(os.listdir(fwDir)) == 0 or not os.path.isdir(reDir) or len(os.listdir(reDir)) == 0 or o.createFasta:
                print 'Creating bisulfite fasta files %s'%time.asctime()
                create_allSubdirs(bisulfiteDir)
                if not os.path.isdir(genomeDir):
                    download_chromosomeFasta(o.spec, o.ass, genomeDir)
                makeBisulfiteFasta.transform_dir(genomeDir, bisulfiteDir, True, o.ass, o.multi, fwMerged+','+reMerged)
                print 'Done creating bisulfite fasta files %s'%time.asctime()
                
            bisIdx = []
            if o.bowtie_idx is None:
                print 'debug1', o.bisRun
                if int(o.bisRun) == 0 or int(o.bisRun) == 2:
                    print 'debug2'
                    fwIdxDir = 'files/%s/%s/bowtie_indexes/mm9_bisulfite_fw/'%(o.spec, o.ass)
                    fwIdx = os.path.join(fwIdxDir, 'mm9_bisulfite_fw')
                    if not os.path.isdir(fwIdxDir) or len(os.listdir(fwIdxDir))==0 or o.createIndex:
                        create_allSubdirs(fwIdxDir)
                        cmd_queue.put('bowtie-build %s %s'%(fwMerged, fwIdx))
                    bisIdx.append(fwIdx)
                if int(o.bisRun) == 1 or int(o.bisRun) == 2:    
                    print 'debug3'
                    reIdxDir = 'files/%s/%s/bowtie_indexes/mm9_bisulfite_re/'%(o.spec, o.ass)
                    reIdx = os.path.join(reIdxDir, 'mm9_bisulfite_re')
                    if not os.path.isdir(reIdxDir) or len(os.listdir(reIdxDir))<6 or o.createIndex:
                        create_allSubdirs(reIdxDir)
                        cmd_queue.put('bowtie-build %s %s'%(reMerged, reIdx))
                    bisIdx.append(reIdx)  
                bowtieIdx = ','.join(bisIdx)
                cmd_queue.join()
            else:
                if o.bisRun == 2 and not o.bowtie_idx.count(',')==1:
                    print 'WARNING! Two bowtie-indexes needs to be provided for bisulfite runs\n Terminating program'%(o.bowtie_idx+suff)
                    sys.exit()
                bisIdx = o.bowtie_idx.split(',')
                for bIdx in bisIdx:
                    check_index(bIdx)
                bowtieIdx = o.bowtie_idx
            print 'bowtieIdx', bowtieIdx
        
        #TRANSCRIPTOME MULTo RUN
        elif o.transLevel or o.geneLevel: 
            transcriptome = True
            if o.geneLevel:
                runType='geneLevel'
            else: 
                runType = 'transLevel'
            if o.annType is None:
                print 'WARNING! You need to enter an annotation type (ensGene, refGene, knownGene) for transcriptome runs.\nTerminating program.'
                sys.exit()
            if o.annFile is None:
                annFile = 'files/%s/%s/annotationFiles/%s.txt'%(o.spec, o.ass, o.annType)
            else:
                annFile = o.annFile
            convAnnFile = annFile+'.conv'
            transcriptomeDir = 'files/%s/%s/fastaFiles/transcriptomeFasta/%s/'%(o.spec, o.ass, o.annType)
            if o.pairedEnd:
                ii = o.insertInfo.replace(',', '-')
                outDir = 'files/%s/%s/MULfiles/%s_pairedEnd%s_%s_%s_%s-%s'%(o.spec, o.ass, o.ass, ii, o.annType, runType, o.low, o.high)
            else:
                outDir = 'files/%s/%s/MULfiles/%s_%s_%s_%s-%s'%(o.spec, o.ass, o.ass, o.annType, runType, o.low, o.high)
            transcriptomeMerged = os.path.join(transcriptomeDir, '%s_%sMerge.fa'%(o.ass, o.annType))
            if not os.path.isdir(transcriptomeDir) or len(os.listdir(transcriptomeDir)) == 0 or o.createFasta:
                print 'Creating transcriptome fasta files %s'%time.asctime()
                if not os.path.isfile(annFile):
                    annDir = 'files/%s/%s/annotationFiles/'%(o.spec, o.ass)
                    linkPath = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/%s.txt.gz'%(o.ass, o.annType)
                    try:
                        subprocess.call('wget --timestamping -P %s.gz %s'%(annFile, linkPath))
                        subprocess.call('gzip -d %s.gz'%annFile)
                    except:
                        print 'WARNING! Faulty link %s. Check if annotation type is correct, or input annotation file manually. \nTerminating program.'%linkPath
                create_allSubdirs(transcriptomeDir)
                if not os.path.isdir(genomeDir):
                    download_chromosomeFasta(o.spec, o.ass, genomeDir)
                if not os.path.isfile(convAnnFile):
                    convert_annFile(annFile, o.annType, convAnnFile)
                makeTranscriptomeFasta.main(o.ass, convAnnFile, o.annType, o.spec, transcriptomeDir, transcriptomeMerged, None)
                print 'Done creating transcriptome fasta files %s'%time.asctime()
            inFastaDir = transcriptomeDir
            mergedFasta = genomeMerged+','+transcriptomeMerged
            #Bowtie index
            if o.bowtie_idx is None:
                transIdxDir = 'files/%s/%s/bowtie_indexes/%s_%s_transcriptome/'%(o.spec, o.ass, o.ass, o.annType)
                transIdx = os.path.join(transIdxDir, '%s_%s_transcriptome'%(o.ass, o.annType))
                if not os.path.isdir(transIdxDir) or len(os.listdir(transIdxDir))<6 or o.createIndex:
                    create_allSubdirs(transIdxDir)
                    cmd_queue.put('bowtie-build %s %s'%(mergedFasta, transIdx))
                    cmd_queue.join()
                bowtieIdx = transIdx
            else:
                check_index(o.bowtie_idx)
                bowtieIdx = o.bowtie_idx
            #Annotation json file
            if o.jsonFile is None:
                jsonFile = annFile.strip('.txt')+'.json'
            else:
                jsonFile = o.jsonFile

            if not os.path.isfile(jsonFile):
                if not os.path.isfile(annFile):
                    linkPath = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/%s.txt.gz'%(o.ass, o.annType)
                    try:
                        subprocess.call('wget --timestamping -P %s.gz %s'%(annFile, linkPath))
                        subprocess.call('gzip -d %s.gz'%annFile)
                    except:
                        print 'WARNING! Faulty link %s. Check if annotation type is correct, or input annotation file manually. \nTerminating program.'%linkPath
                if not os.path.isfile(convAnnFile):
                    convert_annFile(annFile, o.annType, convAnnFile)
                mappingConversionSimple.refseq2json(convAnnFile, jsonFile)

        #GENOME MULTo RUN
        else: 
            
            inFastaDir = genomeDir
            mergedFasta = genomeMerged
            outDir = 'files/%s/%s/MULfiles/%s_%s-%s'%(o.spec, o.ass, o.ass, o.low, o.high)
            if o.bowtie_idx is None:
                genomeIdxDir = 'files/%s/%s/bowtie_indexes/%s_genome/'%(o.spec, o.ass, o.ass)
                if not os.path.isdir(genomeIdxDir):
                    create_allSubdirs(genomeIdxDir)
                bowtieIdx = os.path.join(genomeIdxDir, '%s_no_random'%o.ass)
            else:
                bowtieIdx = o.bowtie_idx
            print 'bowtie-index file: %s'%bowtieIdx
            if len(os.listdir(genomeIdxDir))<6 or o.createIndex:
                print 'Creating bowtie index'
                cmd_queue.put('bowtie-build %s %s'%(mergedFasta, bowtieIdx))
                cmd_queue.join()
                print 'Done creating bowtie index'
            else:
                check_index(bowtieIdx)
                                
          
        print 'Starting to create uniqueness files %s'%time.asctime()
        query_all_chrs_threaded(inFastaDir, int(o.low), int(o.high), outDir, bowtieIdx, o.ass, bowtie_proc, o.override, int(o.bs), jsonFile, o.transLevel, o.geneLevel, o.keepFiles, o.bisulfite, int(o.bisRun), o.smallChrom, o.missmatches, o.countonly, o.pairedEnd, o.insertInfo, o.dictFilesExist)
        print 'Finished all chr in %s %s'%(o.genome_d, time.asctime())
    
    else:
        o.help = True

    
    if len(sys.argv) < 2 or o.help:
        print 'Non-optional arguments for uniqueness query from file:'
        print '-a followed by assembly (eg hg19)'
        print '-s followed by species (eg hsa)'
        print '-t followed by annotation type, if the run is on a transcriptome (refGene, ensGene, mgcGenes and knownGene formats are supported)'
        
        
        print '\nInputs to change the type of run'
        print '-G to activate gene level run'
        print '-T to activate transcript level run'
        print '-S to bin small chromosomes/transcripts and run together to increase efficiency (default for transcriptome runs)'
        print '-B to activate bisulfite run'
        print '-r followed by bisulfite run mode. 0 will create forward strand uniqueness files, 1 will create reverse strand, and 2 (default) creates both'
        print '-P to activate paired end run'
        print '-C to activate countonly mode (the output will only be number of unique positions in each transcript, instead of positional info. Default for paired end runs.'
        
        print '\nInput for changing properties of the run'
        print '-v followed by number of allowed missmatches. Default = 0.'
        print '-p followed by number of processors to use. Default = 25.'
        print '-k followed by lowest k-mere to query. Default = 20.'
        print '-m followed by highest k-mere to query. Default = 255.'
        print '-l followed by insert length, standard deviation form insert length, and number of randomized insert to be created from each position. \n   Comma separated. Used only for paired end runs. Default = 250,25,5.'
        print '-z followed by blocksize, ie the number of "reads" to be ran together with bowtie.'
                
        print '\nIf input or output files are not in the default directories, they can be specified by following options'
        print '-g followed by genome directory (input fasta files)'
        print '-o followed by output directory'
        print '-b followed by bowtie index. Observe that two indexes (comma separated) is needed for bowtie runs'
        print '-n followed by annotation file (eg a refGene or ensGene file)'
        print '-j followed by a json file created by mappingConversionSimple.py. Contains exon-position info for transcripts'
        print '-c followed by chromosome fasta file'
        print '-e followed by the directory for the block to be processed. (Should never need manual input.)'
        
        print '\nOptional arguments for overwriting files and resuming runs'
        print '-O to overwrite existing uniqueness files in the given directory'
        print '-D to indicate that the fasta-dictionary and output-dictionary files are already created (will save runtime if applicable)'
        print '-R to indicate that a crashed run is resumed to avoid having to restart from beginning'
        print '-I to create new bowtie indexes'
        print '-F to create new fasta files (for bisulfite genome or transcripts)'
        print '-A to create both bowtie index and fasta files'

        print '\nHelp and debugging'
        print '-K to keep all intermediate files. ONLY USE FOR DEBUGGING!'
        print '--print-help to see this help message'
        print 'Observe that when creating fasta files for transcripts or bisulfite genome, the flag -g is used to input the original chromosomal fasta files.'


