"""
Mapping transcriptome to genome mappings reads in sam format.

"""

import os, sys, json, string
import optparse


__version__ = 0.1

class tx2genome:
    
    #chromosomeDict = { 'mm9':['chrM', 'chrX','chrY','chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19'],
    #               'hg19':['chrM', 'chrX','chrY','chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19', 'chr20', 'chr21','chr22']}

    def __init__(self, txfile, assembly):
        self.txfile = txfile
        self.assembly = assembly
        self._read_txfile()
        self.chromosomes = self.tx.get('allChromosomes').keys()

        
    def _read_txfile(self):
        with open(self.txfile, 'r') as f:
            self.tx = json.load(f)

    def convert_read(self, read, readmatchlength):
        genomicregions=[]
        txname = read[2]
        tx = self.tx.get(txname, None)
        if tx is None: 
            if txname in self.chromosomes:
                genomicread = [read[0],read[1],txname,str(int(read[3])),str(int(read[3])+readmatchlength-1)]
                return genomicread
            else:
                sys.stderr.write("Transcript name not in database: %s\n"%txname)
                return None
        txstrand = tx[1]
        readtxstart = int(read[3])#will be zero-based
        exoncount = len(tx[2])
        chrom = tx[0]

        # get genomic coordinate for read start
        if txstrand == '+':
            clength = 0
            for exonid in xrange(exoncount):
                elength = tx[2][exonid][1] - tx[2][exonid][0] + 1
                if readtxstart < clength + elength:
                    #Will create 0-based readstart coordinate.
                    genomicreadstart = tx[2][exonid][0] + (readtxstart-clength-1)
                    break
                clength += elength
            try:
                v = genomicreadstart
            except UnboundLocalError:
                sys.stderr.write("no genomicreadstart found for: %s\n" % read)
                return None

                
        else:
            clength=0
            for exonid in xrange(exoncount):
                elength = tx[2][exoncount-exonid-1][1] - tx[2][exoncount-exonid-1][0] + 1#start with the last exon since this is the first exon in a reverse strand transcript
                if readtxstart < clength + elength:
                    # pos within this exon segment
                    genomicreadend = tx[2][exoncount-exonid-1][1] - (readtxstart-clength)-1
                    break
                clength += elength
            try:
                v = genomicreadend
            except UnboundLocalError:
                sys.stderr.write("no genomicreadend found for: %s\n" % read)
                return None

        if txstrand == '+':
            #exonEnd -genomicreadstart+1 longer than the read => read does not span junction
            if tx[2][exonid][1]-1 - genomicreadstart +1 >= readmatchlength or exonid == len(tx[2])-1:
                # read contained within exon
                genomicreadend = genomicreadstart+readmatchlength-1   
                genomicregions.append( (genomicreadstart, genomicreadend) )
            else:                                                                                   
                # read spans one or more junctions
                inc_exon = 1
                exonmatchinglength = tx[2][exonid][1] - genomicreadstart
                remaininglength = readmatchlength - exonmatchinglength                                                                    
                genomicregions.append( (genomicreadstart, tx[2][exonid][2]) )
                
                newexonlength = tx[2][exonid+inc_exon][1] - tx[2][exonid+inc_exon][0] + 1                                                 
                while remaininglength > newexonlength and exonid+inc_exon+1 != len(tx[2]):
                    genomicregions.append((tx[2][exonid+inc_exon][0], tx[2][exonid+inc_exon][1]))
                    inc_exon += 1
                    remaininglength -= newexonlength
                    newexonlength = tx[2][exonid+inc_exon][1] -tx[2][exonid+inc_exon][0] + 1
                # remaining part
                #First minus:make exonStart zerobased, second minus: because exonStart is included
                genomicreadend = tx[2][exonid+inc_exon][0]-1 + remaininglength-1
                genomicregions.append((tx[2][exonid+inc_exon][0], genomicreadend))
            
        else:
            #Plus one: because genomicreadstart is included in length. >= because startPos can be exonStart. Minus one: make exonStart zerobased
            if genomicreadend - readmatchlength+1 >= tx[2][exoncount-exonid-1][0]-1:
                # read contained within exon
                genomicreadstart = genomicreadend - readmatchlength+1
                genomicregions.append( (genomicreadstart, genomicreadend) ) 

            else:
                inc_exon = 1
                exonmatchinglength = genomicreadend - (tx[2][exoncount-exonid-1][0]-1) + 1
                remaininglength = readmatchlength - exonmatchinglength
                genomicregions.append( (tx[2][exoncount-exonid-1][0], genomicreadend))
                newexonlength = tx[2][exoncount-exonid-1-inc_exon][1] - tx[2][exoncount-exonid-1-inc_exon][0] + 1
                while remaininglength > newexonlength:
                    genomicregions.append((tx[2][exoncount-exonid-1-inc_exon][0], tx[2][exoncount-exonid-1-inc_exon][1]))
                    inc_exon += 1
                    remaininglength -= newexonlength
                    newexonlength = tx[2][exoncount-exonid-1-inc_exon][1] - tx[2][exoncount-exonid-1-inc_exon][0] + 1
                # remaining part                                                                              
                exon_idx = exoncount-exonid-1-inc_exon
                if exon_idx >= 0:
                    genomicreadstart = tx[2][exon_idx][1]-1 - remaininglength + 1 #Beginning of next exon (1-based)-1 - remaininglength+1
                else:
                    genomicreadstart = genomicreadend - readmatchlength+1
                    print 'WARNING! Genomic read start outside transcript!'
                    print 'Transcript: %s, Transcript start: %i, Genomicreadstart: %i'%(txname, tx[2][0][0]-1, genomicreadstart)

                
        #genomicreadstart,genomicreadend = genomicregions
        genomicread = [read[0],txstrand,chrom,str(genomicreadstart),str(genomicreadend)]
        return genomicread
    
    def convert_transcriptomefile(self, txmapped, genomicmapped, readmatchlength):#genomicmapped_sam is the output-file
        txf = open(txmapped)
        gf = open(genomicmapped, 'w')
        for line in txf:
            read = line.strip('\n').split('\t')#name, strand, chrom, startPos
            try:
                genomicread = self.convert_read(read, readmatchlength)
                if not genomicread is None: 
                    gf.write('\t'.join(genomicread)+'\n')
            except OverflowError:
                sys.stderr.write('%s\n'%read)
        gf.close()

# CREATE JSON FILES

def refseq2json(refseq_f, json_f):
    tx = {}
    tx.setdefault('allChromosomes', {})
    fh = open(json_f, 'w+')
    for line in open(refseq_f):
        p = line.strip().split("\t")
        # save genomic coordinates for each transcript: genecoords = (chr, strand, [(start, end, length), (start, end), (start, end)])
        tx[p[1]] = (p[2], p[3], [(int(st)+1, int(ed), int(ed)-int(st)+1+1) for st, ed in zip(p[9].split(",")[:-1], p[10].split(",")[:-1])] )
        tx['allChromosomes'].setdefault(p[2], 0)
    json.dump(tx, fh)
    fh.close()

if __name__=='__main__':
    opts = optparse.OptionParser()
    opts.formatter=optparse.TitledHelpFormatter()
    opts.add_option('-i','--in-file',dest='inFile')
    opts.add_option('-a','--assembly', dest='ass')
    opts.add_option('-j','--json-file', dest='jsonFile')
    opts.add_option('-o','--out-file', dest='outFile')
    opts.add_option('-J', '--make-json', dest = 'makeJson', default = False, action = 'store_true')
    
    (o, args) = opts.parse_args()
    if not o.inFile is None and not o.outFile is None:
        if o.makeJson:
            refseq2json(o.inFile, o.outFile)#Here inFile is the annotation file.
        else:
            t = tx2genome(o.jsonFile, assembly)
            t.convert_transcriptomefile(o.inFile, o.outFile)
    else:
        o.help = True
    
    if len(sys.argv) < 2 or o.help:
        print 'Non-optional arguments'
        print '-a followed by assembly (eg hg19)'
        print '-i followed by input file. A bowtie output file if making conversions, an annotation file if making a new json-file'
        print '-j followed by a json file created by mappingConversionSimple.py. Contains exon-position info for transcripts'
        print '-o followed by output file'
        print '-J to create a json-file'
        