'''
Created on Jul 25, 2011

@author: helena_storvall
'''
import optparse
import os, sys, subprocess, time
import threading, Queue


def revComp(seq, useRNA=False):
    return 

def reverse_chromosome(chr_f, reverseFile):
    f = open(chr_f)
    print '\tCreating the reverse strand file %s %s'%(reverseFile, time.asctime())
    reverseDict = {'A':'T', 'C': 'G', 'G': 'T', 'T': 'A', 'N': 'N'} #OBS! G translates to T instead of C.
    header = f.readline()
    print header
    wholeChrom = f.read().replace('\n', '').upper()
    print wholeChrom[:10]
    reverseCompBisulf = [reverseDict.get(b, 'N') for b in wholeChrom]
    rf = open(reverseFile, 'w')
    rf.write(header.replace('\n', '')+'_bisulfite_reverse\n')
    while reverseCompBisulf:
        i = 0
        lineList = []
        while i < 50 and reverseCompBisulf:
            lineList.append(reverseCompBisulf.pop())
            i+=1
        rf.write(''.join(lineList)+'\n')
    rf.close()
    print '\tDone creating the reverse strand file %s %s'%(reverseFile, time.asctime())
    
def transform_file(chr_f, outDir, override, ass):
    myChr = chr_f.split('/')[-1].split('.')[0]
    print 'Creating forward bisulfite fasta file for %s %s'%(myChr, time.asctime())
    fwDir = os.path.join(outDir, 'forward')
    if not os.path.isdir(fwDir):
        os.mkdir(fwDir)
    rvDir = os.path.join(outDir, 'reverse')
    if not os.path.isdir(rvDir):
        os.mkdir(rvDir)
    out_f = os.path.join(fwDir, '%s.bisulfite.forward.fa'%(myChr))
    reverseFile = os.path.join(rvDir, '%s.bisulfite.reverse.fa'%(myChr))
    if os.path.isfile(out_f):
        if override:
            print '\tRegenerating the file %s'%out_f
            subprocess.call('rm %s'%out_f, shell = True)
        else:
            print '\tThe file %s has already been generated. Chose override option if you wish to re-generate it.'%out_f
            sys.exit()
    f = open(chr_f)
    of = open(out_f, 'w')
    for line in f:
        if line.startswith('>'): 
            of.write(line.replace('\n', '')+'_bisulfite_forward\n')
            continue
        newLine = ''
        for c in line:
            c = c.upper()
            if c == 'C':
                newLine += 'T'
            else:
                newLine += c
        of.write(newLine)
    of.close()
    f.close()
    print 'Creating reverse bisulfite fasta file for %s %s'%(myChr, time.asctime())
    reverse_chromosome(chr_f, reverseFile)
    print 'Done creating both bisulfite fasta files for %s %s'%(myChr, time.asctime())
    
def merge_files(outDir, fileList, ass, strand, mergedFile):
    print 'Merging files from %s strand %s'%(strand, time.asctime())
    mf = open(mergedFile, 'w')
    for fileName in fileList:
        f = open(fileName)
        for line in f:
            mf.write(line)
        f.close()
    mf.close()
    
def transform_dir(genome_d, outDir, override, ass, multi, mergeNames):           
    # Setup multi-threading support
    cmd_queue = Queue.Queue()
    num_threads = int(multi)
    def process_command(i, q):
        while True:
            cmd = q.get()
            print "[%i] %s" % (i, cmd)
            subprocess.call(cmd, shell=True)
            q.task_done()
    # set up some threads
    for i in range(num_threads):
        print 'starting thread %i '%i+time.asctime()
        worker = threading.Thread(target=process_command, args=(i, cmd_queue,))
        worker.setDaemon(True)
        worker.start()
        
    if not os.path.isdir(outDir):
        os.mkdir(outDir)
    infiles = [f for f in os.listdir(genome_d) if f.endswith(".fa")]
    print infiles
    if override:
        ovr = '-O'
    else:
        ovr = ''
    fwDir = os.path.join(outDir, 'forward')
    if not os.path.isdir(fwDir):
        os.mkdir(fwDir)
    rvDir = os.path.join(outDir, 'reverse')
    if not os.path.isdir(rvDir):
        os.mkdir(rvDir)
    allForward = []
    allReverse = []
    for in_f in infiles:
        if in_f.endswith('merged.fa'):continue
        in_fp = os.path.join(genome_d, in_f)
        myChr = in_f.split('.')[0]
        cmd_queue.put("python src/makeBisulfiteFasta.py -c %s -o %s -a %s %s" % (in_fp, outDir, ass, ovr))    
        allForward.append(os.path.join(fwDir, '%s.bisulfite.forward.fa'%(myChr)))#, '%s_bisulfite.fa'%(myChr)))
        allReverse.append(os.path.join(rvDir, '%s.bisulfite.reverse.fa'%(myChr)))#, '%s_r_bisulfite.fa'%(myChr)))
    cmd_queue.join()
    if mergeNames is None or mergeNames.count(',') == 0:
        fwName = os.path.join(outDir, '%s_%s_bisulfiteMerge.fa'%(ass, 'forward'))
        reName = os.path.join(outDir, '%s_%s_bisulfiteMerge.fa'%(ass, 'reverse'))
    else:
        fwName,reName = mergeNames.split(',')
    merge_files(outDir, allForward, ass, 'forward', fwName)
    merge_files(outDir, allReverse, ass, 'reverse', reName)
    #return forwardMerged, reverseMerged

if __name__=='__main__':
    opts = optparse.OptionParser()
    opts.formatter=optparse.TitledHelpFormatter()
    opts.add_option('-g','--genome-dir',dest='genome_d')
    opts.add_option('-c', '--chromosome', dest= 'chr_f')
    opts.add_option('-o','--out-dir', dest='outDir', default = 'bisulfiteFastaFiles')
    opts.add_option('-p','--processors', dest='multi', default = 1)
    opts.add_option('-m','--merge-names', dest='mergeNames')
    opts.add_option('-a','--assembly', dest='ass')
    opts.add_option('-O', '--override', dest = 'override', default = False, action = 'store_true')
    opts.add_option('--print-help', dest = 'help', default = False, action = 'store_true')
    (o, args) = opts.parse_args()
    
    if not o.outDir is None:
        if not o.chr_f is None:
            print '\tMaking bisulfite genome for %s '%o.chr_f+time.asctime()
            transform_file(o.chr_f, o.outDir, o.override, o.ass)
            print '\tFinished with file: %s '%o.chr_f+time.asctime()
        elif not o.genome_d is None and not o.multi is None:
            print 'Making bisulfite genome for multiple chromosomes'+time.asctime()
            transform_dir(o.genome_d, o.outDir, o.override, o.ass, o.multi, o.mergeNames)
            print 'Finished all chr in %s '%o.genome_d+time.asctime()
        else:
            o.help = True
    
    if len(sys.argv) < 2 or o.help:
        print 'Non-optional arguments'
        print '-g followed by genome directory (input fasta files)'
        print '-c followed by chromosome fasta file'
        print '-a followed by assembly (eg hg19)'
        print '-s followed by species (eg hsa)'
        print 'Optional arguments'
        print '-o followed by output directory (default = transcriptomeFastaFiles)'
        print '-p followed by the number of processors to use'
        print '-O to regenerate files that already exist'
        