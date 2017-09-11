import sys
import subprocess
import argparse
from os import listdir
from os.path import isfile, join
from parser_MAF_lessMem import catchExceptions

parser = argparse.ArgumentParser()
#parser.add_argument("iteration", type=int)
#parser.add_argument("mafDir", type=str)
#parser.add_argument("mafFile",type=str)
parser.add_argument("blockNumber", type=int)
parser.add_argument("speciesDir", type=str)

args = parser.parse_args()

#mafFiles = [join(args.mafDir,f) for f in listdir(args.mafDir) if isfile(join(args.mafDir, f))]
#mafFile = args.mafFile #mafFiles[args.iteration]

listSpeciesFiles = [join(args.speciesDir,f) for f in listdir(args.speciesDir) if isfile(join(args.speciesDir, f))]

blockNum = args.blockNumber

#print("converting maf files to temp files...")
lastLineA = False

writeFiles = list()
for _file in listSpeciesFiles:
    #print('file: {}'.format(_file))
    appendFile = open(_file, 'a')
    #print(type(appendFile))
    writeFiles.append(appendFile)

for line in sys.stdin:

        if line[0] not in {'#','a','s','i','e','q',' ','\n'}:
            #print("Ignoring {}. File not of maf format".format(readFile.name))
            break


        if line[0] == 'a':
            blockNum += 1
            lastLineA = True
            buf = line.split()[1]
            #read score, with and without a decimal place

            str_score = buf.split('=')[1]
            if '.' in str_score:
                score = int(float(str_score))
            else:
                score = int(str_score)

        if line[0] == 's':
            buf = line.split()
            ID = buf[1].split('.')

            try:
                ID[1]
            except IndexError:
                continue #skip lines that do not have chromosome information
            else:
                i = 0
                for writeFile in writeFiles:
                    i+=1
                    numSearches = len(listSpeciesFiles)

                    writeFileName = writeFile.name.split('/')[-1]
                    if writeFileName.startswith(ID[0]):#ex: does dm6_temp.bed start with dm6
                        writeFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
                                        ID[1],ID[0]+'_'+str(blockNum), buf[2], buf[3], buf[4], lastLineA, score))
                        #               chr   species_blocknum         start   length  strand  isRefSpecies allignmentScore
                    elif i > numSearches:
                        raise Exception("the species of a maf allignment was not in the list of species given\n"+\
                                        "allignment species: {}".format(ID[0]))
            lastLineA = False

for openFile in writeFiles:
    openFile.close()

print(blockNum)
'''
print("i: {}".format(args.iteration))
if args.iteration+1 < len(mafFiles):
    print("calling new iteration")
    print("blockNum: {}".format(blockNum))
    catchExceptions(mafFile)
    subprocess.call("zcat "+mafFile+" | python3 testRecur.py "+str((args.iteration+1))+" "+args.mafDir+" "+str(blockNum)+" "+args.speciesDir, shell=True)
else:
    print("exiting out of recur: {}+1 not < {}".format(args.iteration,len(mafFiles)))
    print("blockNum: {}".format(blockNum))
    pass
'''
