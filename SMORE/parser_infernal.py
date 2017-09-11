
#This program takes an infernal output file and turns it into a list of genes
from os import listdir
from os.path import isfile, join
from dataStructures_lessMem import Gene
import subprocess

def parseInfernal(filePath,  geneList=list()):
    '''
    Takes a directory filled with infernal files. parses the infernal files.
    Stores all of the genes of interest marked significant by infernal.
    Returns list of the genes of interest.

    filePathIn = directory filled with infernal files
    '''

    if len(listdir(filePath)) <=0:
        raise Exception("no Infernal outputs to parse")
    else:
        for _file in listdir(filePath):
            filePlusPath = join(filePath,_file)
            geneList.append(list())
            if isfile(filePlusPath):
                try:
                    f = open(filePlusPath, 'r')
                    f.readline()
            
                except IOError:
                    raise Exception("cannot open {}".format(_file))
                except UnicodeDecodeError:
                    raise Exception("try using an unzipped file")
                else:
                    f.seek(0)
                    for line in f:
                        if line.startswith('# INFERNAL'):
                            infernalVersion = line[2:].strip("#\n") #cut off the leading '# ' and strip the end line

                        elif line[0] == '>':
                            chromosome = line.split()[1]
                            lineList = readGene(f)
                            if lineList[7] == "!":

                                #                                               panTro                      chr3         100          10           
                                geneList[len(geneList)-1].append(Gene(f.name.split('/')[-1].split('.')[0], chromosome, lineList[1], lineList[2],\
                                                                      lineList[3], lineList[4], lineList[5], lineList[6], score=lineList[8]))
                                #                                       +             23       2nd structure   Sequence    score
                            else:
                                break
                    f.close()
                    
                    subprocess.call("rm "+filePlusPath, shell=True)
    #print("done")
    return geneList, infernalVersion

def getGeneLength(start, end):
    '''
    returns the length of the gene
    
    if we are on the + strand we read left to right so end is bigger so subtract start from end.
    If we are on the - strand we read right to left so start is bigger so subtract end from start.
    '''
    if end > start:
        return end - start
    return start - end

def getGeneStart(_from, to, strand):
    '''
    takes the start and end coords of a gene, assuming measurement strand was +
    converts to strand specific measuring. 
    (- strand genes measured from - strand start, + strand genes measured from + strand start)

    if strand given is +:
       return start
    if strand given is -:
       return end 
    '''
    if strand == '+':
        return int(_from)
    return int(to)

def readGene(_file):
    #skip two lines
    _file.readline()
    _file.readline()
    #read len, !or?, score, strand, rank, start(+ strand) and end(+strand)
    line1 = _file.readline().split()
    length = getGeneLength(int(line1[9]),int(line1[10]))
    #skip two lines
    _file.readline()
    _file.readline()
    #read secondary structure + get number of spaces before secondary structure appears starts
    line2 = _file.readline()
    spaceCount = 0
    for char in line2:
        if char == ' ':
            spaceCount += 1
        else:
            break
    line2 = line2.split()
    #skip two lines
    _file.readline()
    _file.readline()
    #skip junk before sequence
    _file.read(spaceCount)
    #read sequence
    line3 =_file.readline().split()
    if len(line3) > 2:
        sequence = ''.join(line3[:-1]) #join everything in the line together except the last space separated element
    else:
        sequence = line3[0]
    #       chromosome             start                          length    strand     id num
    return (line3[0], getGeneStart(line1[9],line1[10],line1[11]), length, line1[11], int(line1[0].strip('()')),\
            line2[0], sequence, line1[1], line1[3])
    #secondary structure, sequence  !or?  score
    
if __name__ == '__main__':
    parseInfernal("../test/infernalIn")
