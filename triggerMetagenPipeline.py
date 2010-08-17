import os.path
import commands, os, getopt, sys
__author__="maslen"
__date__ ="$Jul 23, 2010 10:29:08 AM$"

def generateFileList(inputPath):
    fileList = []
    if os.path.isfile(inputPath):
        fileList.append(inputPath)
        return fileList
    elif os.path.isdir(inputPath):
        for file in os.listdir(inputPath):
            absfile = os.path.join(inputPath, file)
            if os.path.isfile(absfile):
                fileList.append(absfile)
    return fileList

def buildRunBsubCommand(pythonPath, outputPath, scriptPath, fileList):
  for fileName in fileList:
    bsubCommand = 'bsub -q production -M 8000 -R "rusage[mem=8000]" ' + pythonPath + ' ' + scriptPath + "/metagenomicsPipeline.py -s " + scriptPath + " -o " + outputPath + " -p " + fileName
    print bsubCommand
    os.system(bsubCommand)

def checkPathSlash(directory):
    if not directory.endswith('/'):
        directory += '/'
    return directory


def usage():
  print "usage: triggerMetagenPipeline -v python_path -o output_path -s script_path -p abs_file_path -> can be absolute path to a file, or to a directory of files to be analysed."

if __name__ == "__main__":

    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:o:s:p:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    inputPath = None
    outputPath = None
    scriptPath = None
    pythonPath = None

    for o, a in opts:
        if o == "-p":
            inputPath = a
        elif o == "-v":
            pythonPath = a
        elif o == "-s":
            scriptPath = a
        elif o == "-o":
            outputPath = a

    if not inputPath or not pythonPath or not outputPath or not scriptPath:
        usage()
        sys.exit(2)

    outputPath = checkPathSlash(outputPath)
    #scriptPath = checkPathSlash(scriptPath)

    fileList = generateFileList(inputPath)

    buildRunBsubCommand(pythonPath, outputPath, scriptPath, fileList)