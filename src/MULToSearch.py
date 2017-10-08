import os, sys, array, argparse
import ConfigParser
import optparse
import copy

class QueryUniqueome():
    def __init__(self, spec, ass, MULToDir, bisulfite):        
        self.binary_folder= MULToDir
        self.chromFileDict = {}
        self.reverseChromFileDict = {}
        self.nonSupportedChr = {}
        self.bisulfite = bisulfite
        self.reverse_binary_folder = None
        #if bisulfite:
        if bisulfite:
            print 'Bisulfite run'
            self.bisulfite = True
            originalFolder = copy.deepcopy(self.binary_folder)
            if not os.path.isdir(os.path.join(originalFolder, 'fw')):
                print 'The forward directory %s for bisulfite query does not exist. It is required for the script to run.\n Terminating program.'%os.path.join(originalFolder, 'forward')
                sys.exit()
            else:
                fw_files = [ f for f in os.listdir(os.path.join(originalFolder, 'fw')) if f.endswith('btxt')]
                for f in fw_files:
                    chrom = f.split('_')[0]
                    self.chromFileDict[chrom] = f
                self.binary_folder = os.path.join(originalFolder, 'fw')
            if not os.path.isdir(os.path.join(originalFolder, 'rv')):
                print 'WARNING! The reverse directory %s for bisulfite query does not exist. Reverse querys will be approximated from forward strand uniqueness.'%os.path.join(originalFolder, 'reverse')
            else:
                rv_files = [ f for f in os.listdir(os.path.join(originalFolder, 'rv')) if f.endswith('btxt')]
                if not len(fw_files) == len(rv_files):
                    print 'WARNING! The reverse and forward directory does not contain same number of chromosomes.'
                    print 'Forward: %s'%','.join(fw_files)
                    print 'Reverse: %s'%','.join(rv_files)
                    print 'Reverse querys will be approximated from forward strand uniqueness.'
                else:
                    for f in rv_files:
                        chrom = f.split('_')[0]
                        self.reverseChromFileDict[chrom] = f
                    self.reverse_binary_folder = os.path.join(originalFolder, 'rv')
        else:
            all_files = [ f for f in os.listdir(self.binary_folder) if f.endswith('btxt')]
            for f in all_files:
                chrom = f.split('.')[0]
                self.chromFileDict[chrom] = f
            
    
    def get_chrFile(self, chr, strand):
        if self.bisulfite:
            if strand == '-' and not self.reverse_binary_folder is None:
                fName = os.path.join(self.reverse_binary_folder, self.reverseChromFileDict.get(chr, 'no_file'))
            else:
                fName = os.path.join(self.binary_folder, self.chromFileDict.get(chr, 'no_file'))
        else:
            fName = os.path.join(self.binary_folder, self.chromFileDict.get(chr, 'no_file'))
        if not os.path.isfile(fName):
            print fName
            if not self.nonSupportedChr.has_key(chr):
                self.nonSupportedChr[chr] = 0 #Save all non-supported so they are not printed more than once.
                print 'The chromosome "%s" is not supported'%chr
                print 'Supported chromosomes: %s'%(', '.join(self.chromFileDict.keys()))
            return None, None
        fileSize = os.path.getsize(fName)
        return open(fName, 'rb'), fileSize
    
    def get_minus_data(self, infh, start, end, strand, kLow, kHigh, fileSize):
        #print 'Converting from plus to minus strand uniqueness.'
        plusStart = max(0, start-kHigh+1)#If plusStart is smaller than zero, it will be set to zero.
        noPositions = end-plusStart+1
        infh.seek(plusStart, 0)
        plusData = map(ord, infh.read(noPositions))
        plusData.reverse()
        minusData = []
        mp = 0
        while mp < end-start+1:
            foundLim = False
            outside = False
            for k in range(kLow, kHigh+1):
                plusPos = mp+k-1
                if plusPos >= len(plusData):#This happens if start is close to the beginning of the chromosome. Then plusPos can be outside the chromosome.
                    outside = True
                    break
                uLengthPlus = plusData[plusPos]
                if uLengthPlus <= k and uLengthPlus != 0:
                    minusData.append(k)
                    foundLim = True
                    break
            mp +=1
            if outside:
                plusPos = len(plusData)-1
                k = plusPos-mp+1
                uLengthPlus = plusData[plusPos]
                #print k, plusPos, mp
                if uLengthPlus <= k and uLengthPlus != 0:
                    minusData.append(kLow)
                    foundLim = True
            if not foundLim:
                minusData.append(0)
        return minusData
    
    def positions_unique(self, data, kmer_threshold):
        'finding unique positions at readlength %i'%kmer_threshold
        uniqueList = []
        for v in data:
            if v <= kmer_threshold and v != 0:
                uniqueList.append(kmer_threshold)
            else:
                uniqueList.append(0)
        return uniqueList       
        
    def read_block(self, infh, start, end, strand, kLow, kHigh, fileSize, chr):
        if start < 0 or end > fileSize-1:
            print 'The query is out of range. Chromosome (%s) coordinates: %i-%i. Query coordinates: %i-%i'%(chr, 0,fileSize-1, start, end)
            return None
        if strand == '-' and self.reverse_binary_folder is None:
            data = self.get_minus_data(infh, start, end, strand, kLow, kHigh, fileSize)
        else:
            infh.seek(start, 0)
            data = map(ord, infh.read(end-start+1))
            if kLow == kHigh:
                data = self.positions_unique(data, kLow)
        return data

    def unique_within_block(self, data, kLow = 20, kHigh = 255):
        #print 'running positions_unique_range'
        reducedData = []
        noUnique = {}
        for k in range(kLow, kHigh+1):
            maxPos = len(data)-k
            if maxPos+1 > 0:
                noUnique.setdefault(k, [0, maxPos+1])
        for i in range(0, len(data)):
            v = data[i]
            if v == 0: 
                reducedData.append(v)
                continue
            for k in range(v, kHigh+1):
                maxPos = len(data)-k
                if i <= maxPos:
                    #maxPos+1 = number of positions the read of length k can be placed in. (allPos)
                    noUnique[k][0]+=1
            if len(data)-i >= v and v <= kHigh:
                reducedData.append(v)
            else:
                reducedData.append(0)
                #The position will be unique for all reads longer than v
                #To only query reads within the region, startPosition of the read can never 
                #be higher than regionEnd-readLength
        for val in noUnique.values():
            val.append(self.fraction_unique(val[0], val[1]))
        #OBS! If kLow > len(data) the dictionary will be empty
        return noUnique, reducedData #{kMere:[numberOfUnique, totalPositions], ...} 
    
    def unique_full_block(self, data, kLow, kHigh):
        noUnique = {}
        for k in range(kLow, kHigh+1):
            noUnique.setdefault(k, [0, len(data)])
        for v in data:
            if v == 0: continue
            for k in range(max(v, kLow), kHigh+1):
                noUnique[k][0]+=1
        for val in noUnique.values():
            val.append(self.fraction_unique(val[0], val[1]))
        return noUnique


    def number_unique(self, uniqueList):
        return [sum(uniqueList), len(uniqueList)]

    def fraction_unique(self, uniquePos, allPos):
        return float(uniquePos)/float(allPos)

    
    def single_kmer_unique(self, chr, start, end, strand, kmer_threshold, wholeBlock):
        #print 'Searching for uniqueness at readlength %i'%kmer_threshold
        return self.kmer_range_unique(chr, start, end, strand, kmer_threshold, kmer_threshold, wholeBlock)
  
    def kmer_range_unique(self, chr, start, end, strand, kLow, kHigh, wholeBlock):   
        if start > end:
            print 'Warning! Start coordinate (%i) was higher than end coordinate (%i). Switching coordinates.'%(start, end)
            start, end = end, start
        start = start -1
        end = end-1
        fh, fileSize = self.get_chrFile(chr, strand)
        if fh is None:
            return None, None
        if not self.reverse_binary_folder is None and strand == '-':
            tmp_start = copy.deepcopy(start)
            tmp_end = copy.deepcopy(end)
            start = fileSize - tmp_end-1
            end = fileSize - tmp_start-1
        #print 'start: %i, end: %i'%(start,end)
        data = self.read_block(fh, start, end, strand, kLow, kHigh, fileSize, chr)
        if data == None:
            return None, None
        if wholeBlock:
            noUnique = self.unique_full_block(data, kLow, kHigh)
        else:
            noUnique, data = self.unique_within_block(data, kLow, kHigh)
        return noUnique, data
    
    def bed2coord(self, bedline):
        p=bedline.strip().split("\t")
        chrom = p[0]
        start = int(p[1])+1#convert to one-based coordinate
        end = int(p[2])
        if len(p) >= 6:
            strand = p[5]
            name = p[3]
        else:
            strand = '+'
            name = '%s:%i-%i'%(chrom, start, end)
        return [chrom,strand,start,end,name]
    
    def gff2coord(self, gffline):
        p=gffline.strip().split("\t")
        chrom = p[0]
        start = int(p[3]) #convert to zero-based coordinates
        end = int(p[4])
        if len(p) >= 7:
            strand = p[6]
        else:
            strand = '+'
        if len(p) >= 9:
            name = p[8].split(';')[0]
        else:
            name = '%s:%i-%i'%(chrom, start, end)
        return [chrom,strand,start,end,name]
    
    
    def iterate_over_file(self, inFile, outFile, kmer_threshold, positions, counter, wholeBlock, fileType):
        """Used to iterate over bed or gff files, reporting either the fraction of unique positions, the location of unique positions, or both to a bed-file"""
        of = open(outFile, 'w')
        f = open(inFile)
        of.write('#chrom\tstart\tend\tregionName\tscore\tstrand')
        if counter:
            of.write('\treadLength\tunique_pos\tall_pos\tfraction_unique')
        if positions:
            of.write('\tpositionalUniqeness')
        of.write('\n')
        for line in f:
            if line[0] == '#': continue
            if inFile.endswith("bed") or fileType == 'bed':
                chrom,strand,start,end,name = self.bed2coord(line)
            elif inFile.endswith("gff") or inFile.endswith("gtf") or fileType == 'gff':
                chrom,strand,start,end,name = self.gff2coord(line)
            else:
                print 'The program only supports bed or gff input files. Make sure you have entered fileType flag if the file ending is not .bed, .gff or .gtf'
                sys.exit()
            noUnique, data = self.single_kmer_unique(chrom, start, end, strand, kmer_threshold, wholeBlock)
            if noUnique == None: continue
            outLine = '%s\t%i\t%i\t%s\t.\t%s'%(chrom, start-1, end, name, strand)
            if counter:
                for key, val in noUnique.iteritems():
                    outLine += '\t%i\t%i\t%i\t%.2f'%(key, val[0], val[1], val[2])
            if positions:
                positionLine = ''
                for u in data:
                    if u == 0:
                        positionLine += '0'
                    else:
                        positionLine += '1'
                outLine += '\t'+positionLine
            outLine += '\n'
            of.write(outLine)
        of.close()
        f.close()


if '__main__' == __name__:

    # Parse command line options                                                            
    opts = optparse.OptionParser()
    opts.formatter=optparse.TitledHelpFormatter()
    opts.add_option('-i','--infile', dest = 'inFile')
    opts.add_option('-o','--outfile', dest = 'outFile')
    opts.add_option('-g','--file-type', dest ='fileType')
    opts.add_option('-u','--multo-dir', dest ='MULToDir')
    #opts.add_option('-f','--configurationfile', dest = 'confFile', default = 'src/uniqueome.conf')
    opts.add_option('-s','--startcoord', dest = 'start', type=int)
    opts.add_option('-e','--endcoord', dest = 'end', type=int)
    opts.add_option('-t','--strand', dest = 'strand')
    opts.add_option('-c','--chromosome', dest = 'chrom')
    opts.add_option('-a','--assembly', dest = 'ass')
    opts.add_option('-p','--species', dest = 'spec')
    opts.add_option('-k','--kmer-length', dest = 'kmerThresh')#, type = int)
    opts.add_option('-m','--kmer-min', dest ='kMin', default = '20')#, type = int)
    opts.add_option('-x','--kmer-max', dest ='kMax', default = '255')#, type = int)
    opts.add_option('-W', '--whole-block', dest = 'wholeBlock', default=False, action="store_true")
    opts.add_option('-B','--bisulfite',dest='bisulfite', default = False, action = 'store_true')
    opts.add_option('--counter', dest= 'counter', default=False, action="store_true")
    opts.add_option('--positions', dest = 'positions', default=False, action="store_true")
    opts.add_option('--all', dest = 'all', default=False, action="store_true")
    opts.add_option('--print-help', dest = 'help', default = False, action = 'store_true')
    #args = opts.parse_args()
    
    (o, args) = opts.parse_args()

    # read in part of chromosome or full
    
    if o.MULToDir is None:
        if o.bisulfite:
            md = 'files/%s/%s/MULfiles/%s_bisulfite_%s-%s/MULTo_files/'%(o.spec, o.ass, o.ass, o.kMin, o.kMax)
        else:
            md = 'files/%s/%s/MULfiles/%s_%s-%s/MULTo_files/'%(o.spec, o.ass, o.ass, o.kMin, o.kMax)
    else:
        md = o.MULToDir
    if not os.path.isdir(md) or os.listdir(md)==0:
        print 'ERROR: The uniqueness files for the selected assembly (%s) could not be found.\n '%o.ass
        print 'Non-existing or empty folder: %s'%md
    else:
        uq = QueryUniqueome(o.spec, o.ass, md, o.bisulfite)
        
    if not o.start is None and not o.end is None and not o.chrom is None and not o.ass is None and not o.strand is None:
        if o.strand is None:
            print 'WARNING! No strand information given. Forward strand will be used by default'
            o.strand = '+'
        if (o.counter == False and o.positions == False) or o.all == True:
            o.counter = True
            o.positions = True
        if not o.kmerThresh is None:
            print 'Querying region with kmer threshold = %s'%o.kmerThresh
            noUnique, data = uq.single_kmer_unique(o.chrom, o.start, o.end, o.strand, int(o.kmerThresh), o.wholeBlock)
        else:
            print 'Querying region with kMin = %s and kMax = %s'%(o.kMin, o.kMax)
            noUnique, data = uq.kmer_range_unique(o.chrom, o.start, o.end, o.strand, int(o.kMin), int(o.kMax), o.wholeBlock)

            
        print 'Limit\tunique_pos\tall_pos\tfraction_unique'
        if o.counter:
            for key in sorted(noUnique.iterkeys()):
                val = noUnique[key]
                print '%i\t%i\t%i\t%.2f'%(key, val[0], val[1], val[2])
        if o.positions:
            if not o.kmerThresh is None:
                positionLine = ''
                for u in data:
                    if u == 0:
                        positionLine += '0'
                    else:
                        positionLine += '1'
                print positionLine
            else:
                print data
    elif not o.inFile is None and not o.outFile is None:
        if o.counter == False and o.positions == False and o.all == False:
            o.counter = True
        elif o.all == True:
            o.counter = True
            o.positions = True
        if o.kmerThresh is None:
            o.kmerThresh = 50
            print 'Querying regions in file with default threshold k = 50'
        else:
            print 'Querying regions in file with kmer threshold = %s'%o.kmerThresh
        uq.iterate_over_file(o.inFile, o.outFile, int(o.kmerThresh), o.positions, o.counter, o.wholeBlock, o.fileType)
    else:
        o.help = True
    
    if len(sys.argv) < 2 or o.help:
        print 'Non-optional arguments for uniqueness query from file:'
        print '-o followed by output file'
        print '-i followed by input file (bed, gff or gtf format)'
        print 'Non-optional arguments for uniqueness query in a single region:'
        print '-s followed by start-coordinate (one-based)'
        print '-e followed by end-coordinate (one-based)'
        print '-c followed by chromosome'
        print '-t followed by strand ("+" or "-")'
        print '-a followed by assembly (eg hg19)'
        print '-p followed by species (eg hsa)'
        print 'Optional arguments:'
        print '-g followed by file type (bed or gff). Use this if your file has the right format, but does not end with .bed, .gff or .gtf'
        print '-k followed by the read length you want to query. Only for single read length queries.'
        print '-m followed by the minimum read length you want to query. Only for queries in a range of read lengths.'
        print '-x followed by the maximum read length you want to query. Only for queries in a range of read lengths.'
        print '-W to get the uniqueness of the whole region, not requiring reads to fit within the region.'
        print 'Output modes:'
        print '--counter to get a count of unique and non-unique position and the fraction between them as output'
        print '--positions to get an array showing uniqueness in each position as output'
        print '--all to activate both --counter and --positions'
        print 'For single region queries, "all" is default output mode. For file queries, "counter" is default'
        


