import os.path
from Bio import SeqIO
import os, sys
import re, commands
import time, getopt
from time import strftime, gmtime
from hashlib import md5
__author__="maslen"
__date__ ="$Jun 27, 2010 12:06:00 PM$"

analysisDirectory = None
sampleIdent = None

def readSequencefile(inputSeqFile, outputFASTAfile, fileFormat):
    dictFormat = {'sff':'sff-trim', 'fasta':'fasta', 'fastq':'fastq'}
    parseFormat = dictFormat.get(fileFormat)
    if not parseFormat:
        print 'Unable to determine file format for parsing'
        sys.exit(2)
    records = (rec for rec in SeqIO.parse(inputSeqFile, parseFormat) if len(rec.seq) > 99)
    SeqIO.write(records, outputFASTAfile, "fasta")

#def getUniqueSequenceList(uclustExecutable, outputUCLUSTfile, outputUNIQfile):
#    os.system(uclustExecutable)
#    inputFile = open(outputUCLUSTfile, 'r')
#    outputFile = open(outputUNIQfile, 'w')
#    #dump file for dev purposes only -> will this be audited??
#    #dumpFile = open("/nfs/nobackup/interpro/development/metagenomics/FKO5Z8U02.uniq.dump", 'w')
#    seedCharacter = 'S'
#    line = inputFile.readline()
#    while line:
#        if line.startswith(seedCharacter):
#            try:
#                outputFile.write(line.split()[8]+'\n')
#            except:
#                print 'Sequence ID not found:', line
#        #output has lines starting with 'H' and 'C' -> 'C' denotes clusters
#        #what needs auditing? e.g:
#        #if line.startswith('H'):
#        #    dumpFile.write(line.split()[8]+'\n')
#        line = inputFile.readline()
#    inputFile.close()
#    outputFile.flush()
#    outputFile.close()

def getUniqueSequenceList(uclustExecutable, outputUCLUSTfile):
    os.system(uclustExecutable)
    inputFile = open(outputUCLUSTfile, 'r')
    returnList = []
    #dump file for dev purposes only -> will this be audited??
    #dumpFile = open("/nfs/nobackup/interpro/development/metagenomics/FKO5Z8U02.uniq.dump", 'w')
    seedCharacter = 'S'
    line = inputFile.readline()
    while line:
        if line.startswith(seedCharacter):
            try:
                returnList.append(line.split()[8])
            except:
                print 'Sequence ID not found:', line
        #output has lines starting with 'H' and 'C' -> 'C' denotes clusters
        #what needs auditing? e.g:
        #if line.startswith('H'):
        #    dumpFile.write(line.split()[8]+'\n')
        line = inputFile.readline()
    inputFile.close()
    return returnList

#def produceUniqueFasta(outputFASTAfile, outputUNIQ_FASTAfile, outputUNIQfile):
#    uniqueIDdata = open(outputUNIQfile).read()
#    uniqueFastahandle = open(outputUNIQ_FASTAfile, "w")
#    #dump file for dev purposes only -> audited??
#    #uniqueFastaDumphandle = open("/nfs/nobackup/interpro/development/metagenomics/FKO5Z8U02_uniq_dump.fasta", "w")
#    for record in SeqIO.parse(outputFASTAfile, "fasta"):
#        if record.id in uniqueIDdata:
#            SeqIO.write(record, uniqueFastahandle, "fasta")
#        #audit those that are filtered out:
#        #else:
#        #    SeqIO.write(record, uniqueFastaDumphandle, "fasta")

def selectFastaSeqsFromList(inputFASTAfile, outputFASTAfile, inputIDlist):
    #dump file for dev purposes only -> audited??
    #uniqueFastaDumphandle = open("/nfs/nobackup/interpro/development/metagenomics/FKO5Z8U02_uniq_dump.fasta", "w")
    record_index = SeqIO.index(inputFASTAfile, "fasta")
    records = (record_index[id] for id in inputIDlist)
    SeqIO.write(records, outputFASTAfile, "fasta")

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def generateChunkedFastaFiles(fileHandle, outputDirectory, format, chunkSize):
    try:
        os.system("rm " + outDirectory + "/*")
    except:
        print("Chunked directory empty.")
    record_iter = SeqIO.parse(fileHandle, format)
    chunkedFiles = []
    for i, batch in enumerate(batch_iterator(record_iter, chunkSize)) :
        outFile = "%s/%s_%s.%s" % (outDirectory, batch[0].id, batch[-1:][0].id, format)
        handle = open(outFile, "w")
        SeqIO.write(batch, handle, format)
        #abspath may be overkill as currently using absolute paths
        chunkedFiles.append(os.path.abspath(outFile))
        handle.close()
    return chunkedFiles

def generateChunkedUniqueFastaFiles(fileHandle, outputDirectory, format, chunkSize):
    try:
        os.system("rm " + outDirectory + "/*")
    except:
        print("Chunked directory empty.")
    record_iter = SeqIO.parse(fileHandle, format)
    chunkedFiles = []
    for i, batch in enumerate(batch_iterator(record_iter, chunkSize)) :
        outFile = "%s/%s_%s.%s" % (outDirectory, batch[0].id, batch[-1:][0].id, format)
        handle = open(outFile, "w")
        SeqIO.write(batch, handle, format)
        #abspath may be overkill as currently using absolute paths
        chunkedFiles.append(os.path.abspath(outFile))
        handle.close()
        checkFastaDuplicates(outFile)
    return chunkedFiles

def removeFastaDuplicates(filein, fileout):
    dictSeqs = {}
    listIds = []
    for record in SeqIO.parse(filein, "fasta"):
        dictSeqs[str(record.seq)] = record.id
    for seq in dictSeqs.keys():
        listIds.append(dictSeqs.get(seq))
    selectFastaSeqsFromList(filein, fileout, listIds)

def checkFastaDuplicates(file):
    print 'Checking duplicates......'
    dictseq = {}
    dictmd5 = {}
    listmd5 = []
    for record in SeqIO.parse(file, "fasta"):
        listseqq = []
        if dictseq.has_key(record.id):
            listseqq = dictseq.get(record.id)
        listseqq.append(record.seq)
        dictseq[record.id] = listseqq
        seq_MD5 = md5(str(record.seq)).hexdigest().upper()
        listseqqq = []
        if dictmd5.has_key(seq_MD5):
            listseqqq = dictmd5.get(seq_MD5)
        listseqqq.append(record.id)
        dictmd5[seq_MD5] = listseqqq
    listmd5 = dictmd5.keys()
    for entry in listmd5:
        listy = dictmd5.get(entry)
        if len(listy) > 1:
            'MD5:', entry
            for fred in listy:
                print '\t', fred, dictseq.get(fred)

def generateLSFcommands(chunkedFiles, commandExec):
    lsfJobRandomName = randomWord(6)
    uniqueNameCount = 1
    listLsfJobNames = []
    for fileChunked in chunkedFiles:
        memUsage = '-M 4000 -R "rusage[mem=4000]"'
        uniqueJobName = "".join([lsfJobRandomName, "_", str(uniqueNameCount)])
        commandExecutable = commandExec.format(fileChunked)
        commandString = "bsub -q production -o {0}.out -e {0}.err {1} -J {2} {3}\n".format(fileChunked, memUsage, uniqueJobName, commandExecutable)
        os.system(commandString)
        listLsfJobNames.append(uniqueJobName)
        uniqueNameCount += 1
    time.sleep(30)
    return listLsfJobNames

def generateLSFcommands2(chunkedFiles, commandExec):
    lsfJobRandomName = randomWord(6)
    uniqueNameCount = 1
    dictLSFJobs = {}
    for fileChunked in chunkedFiles:
        memUsage = '-M 8000 -R "rusage[mem=8000]"'
        uniqueJobName = "".join([lsfJobRandomName, "_", str(uniqueNameCount)])
        commandExecutable = commandExec.format(fileChunked)
        commandString = "bsub -q production -o {0}.out -e {0}.err {1} -J {2} {3}\n".format(fileChunked, memUsage, uniqueJobName, commandExecutable)
        #os.system(commandString)
        dictLSFJobs[uniqueJobName] = [commandString, fileChunked]
        uniqueNameCount += 1
    #time.sleep(30)
    return dictLSFJobs

def randomWord(Length):
    from random import randint
    assert(1 <= Length <= 26) # Verify 'Length' is within range
    charlist = [c for c in "abcdefghijklmnopqrstuvwxyz"]
    for i in xrange(0, Length):
        other = randint(0, 25)
        charlist[i], charlist[other] = charlist[other], charlist[i] # Scramble list by swapping values
    word = ""
    for c in charlist[0:Length]: word += c
    return word.upper()

def checkJobsRunning(showOutput, listJobNames):
    #LSF_STATUS = ['RUN', 'PEND', 'DONE']
    totalJobNumber = len(listJobNames)
    totalCompleted = 0
    for jobName in listJobNames:
        status, output = commands.getstatusoutput('bjobs -J ' + jobName)
        if status != 0:
            print "An error occurred determining the jobs running status on LSF"
            sys.exit(2)
        if output == "".join(["Job <", jobName, "> is not found"]):
            totalCompleted += 1
    if totalCompleted == totalJobNumber:
        return False
    if showOutput:
        currentTime = time.strftime("%H:%M:%S", time.localtime())
        jobsRunning = totalJobNumber - totalCompleted
        if jobsRunning == 1:
            print currentTime, '-> there is 1 LSF job still running'
        else:
            print currentTime, '-> there are', str(jobsRunning), 'LSF jobs still pending/running'
    return True
#delete job if not running , resubmit up to 3 times -> store details in a dictionary

def monitorLSFJobsRunning(listJobNames, dictJobs, dictRestarts):
    LSF_STATUS = ['RUN', 'PEND', 'DONE']
    totalJobNumber = len(listJobNames)
    totalCompleted = 0
    for jobName in listJobNames:
        status, output = commands.getstatusoutput('bjobs -J ' + jobName)
        errorFile = dictJobs.get(jobName)[1] + '.err'
        jobOutFile = dictJobs.get(jobName)[1]+'.out'
        finalFile = dictJobs.get(jobName)[1]+'_out.tsv'
        lsfExecutable = dictJobs.get(jobName)[0]
        jobStatus = None
        if status != 0:
            print "An error occurred determining the jobs running status on LSF"
            sys.exit(2)
        if os.path.isfile(finalFile):
            totalCompleted += 1
        else:
            if output == "".join(["Job <", jobName, "> is not found"]):
                #tidyUpIprscanLSFJobs(errorFile, jobOutFile)
                numberRetries = dictRestarts.get(jobName)
                if numberRetries < 5:
                    os.system(lsfExecutable)
                    dictRestarts[jobName] = numberRetries + 1
                    print 'Restart number', numberRetries + 1, 'for Job', jobName, 'as it failed to complete.'
                    time.sleep(10)
                else:
                    print 'Maximum number of restarts reached - see', jobOutFile, 'for details.'
                    totalCompleted += 1
            else :
                try:
                    jobStatus = output.split('\n')[1:][0].split()[2]
                except:
                    print 'Unable to get LSF job status.'
                if jobStatus not in LSF_STATUS:
                    os.system('bkill -J ' + jobName)
                    numberRetries = dictRestarts.get(jobName)
                    if numberRetries < 5:
                        os.system(lsfExecutable)
                        dictRestarts[jobName] = numberRetries + 1
                        print 'Restart number', numberRetries + 1, 'for Job', jobName, 'as it had', jobStatus, 'status.'
                        time.sleep(10)
                    else:
                        print 'Maximum number of restarts reached - see', errorFile, 'for details.'
                        totalCompleted += 1
    if totalCompleted == totalJobNumber:
        return False
    currentTime = time.strftime("%H:%M:%S", time.localtime())
    jobsRunning = totalJobNumber - totalCompleted
    if jobsRunning == 1:
        print currentTime, '-> there is 1 LSF job still running'
    else:
        print currentTime, '-> there are', str(jobsRunning), 'LSF jobs still pending/running'
    return True

#def tidyUpIprscanLSFJobs(iprscanErrorFile, iprscanJobFile):
#    iprscanPattern = re.compile('iprscan-[0-9]{8}-[0-9]{8}')
#    outputText = None
#    iprscanJobName = []
#    try:
#        outputText = open(iprscanErrorFile).read().split()
#    except:
#        print 'No error file found.'
#    if outputText:
#        for i in outputText:
#            if re.match(iprscanPattern, i):
#                iprscanJobName.append(i)
#        if len(iprscanJobName) == 0:
#            print 'No jobs found to kill.'
#        else:
#            status, output = commands.getstatusoutput('bjobs -l')
#            if status == 0:
#                for line in output.split('\n')[1:]:
#                    for jobName in iprscanJobName:
#                        if line.startswith('Job <') and jobName in line:
#                            jobID = line.split()[1].replace('<', '').replace('>', '').replace(',', '')
#                        #jobName = line.split()[4].replace('<', '').replace('>', '').replace(',', '')
#                            os.system('bkill ' + jobID)
#        try:
#            os.system('rm ' + iprscanErrorFile)
#        except:
#            print 'Error file not found to delete.'
#        try:
#            os.system('rm ' + iprscanJobFile)
#        except:
#            print 'Job output file not found to delete.'
#        time.sleep(10)

def orfProcessingStep(getorfExecutable, orfExtractionExecutable):
    #need to discuss auditing.....
    os.system(getorfExecutable)
    os.system(orfExtractionExecutable)

def getUniquePfamFasta(pfamFastaFile):
    handle = open(pfamFastaFile, 'r')
    pfamData = list(set([line.split()[0] for line in (l for l in open(pfamFastaFile, 'r') if not l.startswith('#'))]))
    pfamData.sort()
    return pfamData

def getNameOnlyFasta(inputFastaFile, outputFastaFile):
    inputHandle = open(inputFastaFile, 'r')
    outputHandle = open(outputFastaFile, 'w')
    for line in inputHandle:
        if line.startswith('>'): outputHandle.write(line.split()[0]+'\n')
        else: outputHandle.write(line)
    inputHandle.close()
    outputHandle.close()

def initialiseFile(footer):
    return "{0}/{1}{2}".format(analysisDirectory, sampleIdent, footer)

def createREADME(readmeFile, inputFile):
    handle = open(readmeFile, 'w')
    today = strftime("%d-%m-%Y")
    current = strftime("%H:%M:%S")
    handle.write("The analyses were started on " + today + " at " + current + ".\n")
    handle.write("The file that was analysed was " + inputFile +  ".\n")
    handle.close()

def usage():
  print "usage: metagenomicsPipeline -s script_path -o output_path -p abs_file_path"

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:o:p:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    start = time.time()
    sampleIdent = None
    inputFilePath = None
    inputFileName = None
    outputPath = None
    scriptPath = None
    fileFormat = None

    for o, a in opts:
        if o == "-s":
            scriptPath = a
        elif o == "-o":
            outputPath = a
        elif o == "-p":
            inputFilePath = a

    if not inputFilePath or not outputPath or not scriptPath:
        usage()
        sys.exit(2)
    
    inputFileName = os.path.basename(inputFilePath)
    if not inputFileName:
        print 'Unable to retrieve file name from path provided, please check the inputs.'
        print usage()
        sys.exit(2)
    
    inputFormats = ['sff', 'srf', 'fastq', 'fasta']
    for formatTag in inputFormats:
        if inputFileName.endswith(formatTag):
            fileFormat = formatTag
    if not fileFormat:
        print "Unable to determine format of input or file, or format not recognised."
        sys.exit(2)
    
    #sampleIdent = ''.join(ch for ch in inputFileName.split('.'+fileFormat)[0] if ch.isalnum()).upper()
    sampleIdent = inputFileName.split('.'+fileFormat)[0].upper()
    if not sampleIdent:
        print "Unable to retrive sample id from file name - check that file name contains letters and/or digits."
        sys.exit(2)

    javaPath = "/ebi/research/software/Linux_x86_64/opt/java/jdk1.6/bin/java"
    analysisDirectory = ''.join([outputPath, sampleIdent])
    if not os.path.exists(analysisDirectory):
        os.system(" ".join(["mkdir", analysisDirectory]))

    dictFiles = {"readme":"{0}/{1}".format(analysisDirectory, "README"),
                 "stepSRF":initialiseFile("_srf.fasta"), "stepInit":initialiseFile(".fasta"), "stepUclust":initialiseFile(".uclust"), "stepInitUniq":initialiseFile(".uniq"),
                 "stepUniqFasta":initialiseFile("_uniq.fasta"), "stepMask":initialiseFile("_masked.fasta"), "stepOrf100":initialiseFile("_orf100.fasta"),
                 "stepOrf100200":initialiseFile("_orf100_200.fasta"), "stepOrf100200Name":initialiseFile("_orf100_200_nameonly.fasta"),
                 "stepPfamTbl":initialiseFile("_pfam.tbl"), "stepPfamFasta":initialiseFile("_pfam.fasta"), "stepPfamUniq":initialiseFile("_pfam_uniq.fasta"),
                 "stepI5tsv":initialiseFile("_I5.tsv")}

    scriptNameLocator = "{0}"
    i5Location = "/nfs/seqdb/production/interpro/development/metagenomics/interproscan/core/jms-implementation/target"
    dictExec = {"srfExec":"{0}/staden/bin/srf2fasta {1} > {2}".format(scriptPath, inputFilePath, dictFiles.get("stepSRF")),
                "uclustExec":"{0}/uclust3.0.617_i86linux32 --id 0.99 --usersort --nucleo --input {1} --uc {2}".format(scriptPath, dictFiles.get("stepInit"), dictFiles.get("stepUclust")),
                "getOrfExec":"{0}/getorf -minsize 100 -find 1 -table 11 {1} -outseq {2}".format(scriptPath, dictFiles.get("stepMask"), dictFiles.get("stepOrf100")),
                "orfExtractExec":"{0}/extract_orfs.pl {1} {2} 200 > {3}".format(scriptPath, dictFiles.get("stepMask"), dictFiles.get("stepOrf100"), dictFiles.get("stepOrf100200")),
                "maskerExec":"'RepeatMasker {0}'",
                "pfamExec":"'/ebi/production/interpro/binaries/64_bit_Linux/HMMER3.0b2/hmmsearch -Z 9421015 -E 100.0 --domtblout {0}_out.txt /ebi/production/interpro/data/members/pfam/24.0/Pfam-A.hmm {0} > /dev/null'",
                "i5Exec":"'{0} -jar -XX:+UseParallelGC -XX:+AggressiveOpts -XX:+UseFastAccessorMethods -Xms512M -Xmx2048M {1}/interproscan-5-dist.dir/interproscan-5.jar -m amqstandalone --appl PfamA --appl PRINTS -o {2}_out.tsv -i {2} -d'".format(javaPath, i5Location, scriptNameLocator)
                }

    createREADME(dictFiles.get("readme"), inputFilePath)

    outDirectory = "".join([analysisDirectory,'/chunked_fasta'])
    if not os.path.exists(outDirectory):
        os.system(" ".join(["mkdir", outDirectory]))
    os.system(" ".join(["chmod", "700", outDirectory]))

#    if fileFormat == 'srf':
#        os.system(dictExec.get("srfExec"))
#        inputFilePath = dictFiles.get("stepSRF")
#        fileFormat = "fasta"

#    readSequencefile(inputFilePath, dictFiles.get("stepInit"), fileFormat)
#    listUniqIDs = getUniqueSequenceList(dictExec.get("uclustExec"), dictFiles.get("stepUclust"))#, dictFiles.get("stepInitUniq"))
#    selectFastaSeqsFromList(dictFiles.get("stepInit"), dictFiles.get("stepUniqFasta"), listUniqIDs)
#    #produceUniqueFasta(dictFiles.get("stepInit"), dictFiles.get("stepUniqFasta"), dictFiles.get("stepInitUniq"))
#    fileHandle = open(dictFiles.get("stepUniqFasta"))
#
#    chunkedFiles = generateChunkedFastaFiles(fileHandle, outDirectory, "fasta", 10000)
#
#    listJobNames = generateLSFcommands(chunkedFiles, dictExec.get("maskerExec"))
#
#    areJobsRunning = checkJobsRunning(False, listJobNames)
#    while areJobsRunning:
#        time.sleep(60)
#        areJobsRunning = checkJobsRunning(True, listJobNames)
#    print 'LSF Submission --> All LSF jobs have completed.'
#
#    os.system("cat " + outDirectory + "/*.masked > " + dictFiles.get("stepMask"))
#
#    orfProcessingStep(dictExec.get("getOrfExec"), dictExec.get("orfExtractExec"))
#
#    fileHandle = open(dictFiles.get("stepOrf100200"))
#    chunkedFiles = generateChunkedFastaFiles(fileHandle, outDirectory, "fasta", 10000)
#    listJobNames = generateLSFcommands(chunkedFiles, dictExec.get("pfamExec"))
#
#    areJobsRunning = checkJobsRunning(False, listJobNames)
#    while areJobsRunning:
#        time.sleep(60)
#        areJobsRunning = checkJobsRunning(True, listJobNames)
#    print 'LSF Submission --> All LSF jobs have completed.'
#
#    os.system("cat " + outDirectory + "/*.txt > " + dictFiles.get("stepPfamTbl"))
#    uniqPfamFastaList = getUniquePfamFasta(dictFiles.get("stepPfamTbl"))
#    getNameOnlyFasta(dictFiles.get("stepOrf100200"), dictFiles.get("stepOrf100200Name"))
#    selectFastaSeqsFromList(dictFiles.get("stepOrf100200Name"), dictFiles.get("stepPfamFasta"), uniqPfamFastaList)

    removeFastaDuplicates(dictFiles.get("stepPfamFasta"), dictFiles.get("stepPfamUniq"))
    fileHandle = open(dictFiles.get("stepPfamUniq"))
    chunkedFiles = generateChunkedUniqueFastaFiles(fileHandle, outDirectory, "fasta", 7000)

    dictLSFJobs = generateLSFcommands2(chunkedFiles, dictExec.get("i5Exec"))
    listJobNames = dictLSFJobs.keys()
    dictRestarts = {}
    listJobNames.sort()
    for job in listJobNames:
        print dictLSFJobs.get(job)[0]
        dictRestarts[job] = 0
        os.system(dictLSFJobs.get(job)[0])
        time.sleep(5)

#    areJobsRunning = checkJobsRunning(False, listJobNames)
#    while areJobsRunning:
#        time.sleep(60)
#        areJobsRunning = checkJobsRunning(True, listJobNames)
#    print 'LSF Submission --> All LSF jobs have completed.'
#    os.system("cat " + outDirectory + "/*.tsv > " + dictFiles.get("stepI5tsv"))

    areJobsRunning = monitorLSFJobsRunning(listJobNames, dictLSFJobs, dictRestarts)
    while areJobsRunning:
        time.sleep(60)
        #print dictRestarts
        areJobsRunning = monitorLSFJobsRunning(listJobNames, dictLSFJobs, dictRestarts)
    print 'LSF Submission --> All LSF jobs have completed.'
    os.system("cat " + outDirectory + "/*_out.tsv > " + dictFiles.get("stepI5tsv"))

    end = time.time()
    timetaken = int(end-start)
    timeminutes =  float(timetaken)/60.0
    handle = open(dictFiles.get("readme"), 'a')
    handle.write("The complete set of analyses took " + str(timetaken) + " seconds, or %.1f minutes!\n" % timeminutes)