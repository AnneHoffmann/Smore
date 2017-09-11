'''
This is the shell script that will run the gene evolution pipeline.

It will take the following arguments:

  1)File pathway to the Genomes that are to be compared
  2)File pathway to the Multiple sequence allignments of those genomes
      These allignments must be labeled identically to the genomes (ex: genomename: panTro, multiple sequence alignment name: panTro)
      The allignment format should be:
        -MAF (file extension '.maf')
          -uncompressed OR
          -compressed with gzip (file extension '.maf.gz', '.maf.Z', '.maf.bz2')

  3)File pathway to a directory where you would like the output to be placed
  4)Name of the reference species used in the multiple sequence alignment as it appears in the multiple sequence alignment / genomes (these should be the same)

Optional Arguments:



  -sg:Genes used in analysis will be found by cmsearch. User must provide a model of the genes they are looking for
  -og:Genes used in analysis are those provided by the user. There should be one file for every species. The files should be named:
          'panTro.bed' or 'dm6.bed', etc.
        Files should be in bed format, with one gene per line of the file.
  -q :Remove a percentage of the genes 

'''
import argparse #python module dealing with command line arguments
import sys
import subprocess #python module allowing other non-python files to be called
from os import listdir #listdir returns a list of everything in a directory
from os.path import isfile, join #join joins a file extension with a file name
from parser_MAF_lessMem import maf2bed #maf2bed parses maf files and writes bed files
from parser_infernal import parseInfernal #parses infernal files and returns a list of Genes
from parser_genes import parseGenes #parses user provied gene files and returns a list of Genes
from helperFunctions import *

parser = argparse.ArgumentParser()
parser.add_argument("multiSeq", type=str, help="directory where the multiple sequence allignments are")
parser.add_argument("outPutDir", type=str, help="directory where you want the output files to be placed")
parser.add_argument("pathToRepo", type=str, help="absolute path of the git repo")
parser.add_argument("referenceSpecies", type=str, help="the name of the reference species as it appears in the "\
                    + "multiple sequence alignment.")
parser.add_argument("python", type=str, help="path to python 3.0 or greater")
parser.add_argument("-sg","--search_genes", action="store", nargs=3,help="search for genes using infernal cmsearch. "\
                    + "First argument must be a file pathway to a RNA Model. "\
                    + "Second argument must be a file pathway to the Genomes."\
                    + "Third argument must be a file path to the infernal.")
parser.add_argument("-og","--own_genes", action="store", help="use prespecified genes in analysis. File path "\
                    + "provided should be a directory. Files in this directory should be named for the species "\
                    + "that it concerns. The names of the species should match those in the multiple sequence "\
                    + "alignments exactly. The files should contain the genes of interest in BED format.")
parser.add_argument("-q","--quality", action="store", default=0, type=int, help="a number "\
                    + "specifying the percentage of blocks that will be thrown away. These blocks will be the lowest "\
                    + "quality blocks. Warning: throwing away too many blocks may lead to an inconclusive analysis.")
parser.add_argument("-incE", action="store", type=float, help="consider sequences <= this E-value threshold as "\
                    + "significant")
parser.add_argument("-incT", action="store", type=float, help="consider sequences >= this score threshold as significant")

args = parser.parse_args()

'''
handel irregularites in args
'''
argList = [args.multiSeq, args.outPutDir, args.pathToRepo]
args.multiSeq, args.outPutDir, args.pathToRepo = endSlash(argList)

#Call helper functions
#full path is now already in argument
#args.pathToRepo = makeFullPath(args.pathToRepo)
makeOutPutDir(args.outPutDir)

#Catch exceptions 
args.quality = float(args.quality)
if args.quality < 0 or args.quality > 99:
    raise Exception("quality value must be a number between 0 and 99 corresponding to the percentage of "\
                    "blocks that will be thrown away.")

if args.own_genes != None and args.search_genes != None:
    raise Exception("both own_genes and search_genes given. own_genes and search_genes are mutually exclusize.")

if len(listdir(args.multiSeq)) == 0:
    raise Exception("No multiple sequence alignments in {}".format(args.multiSeq))

if args.search_genes != None:
    #print("search Genes")

    args.genomes = args.search_genes[1]
    args.model = args.search_genes[0]
    args.infernalPath = args.search_genes[2]
    
    args.genomes, args.infernalPath = endSlash([args.genomes, args.infernalPath])

    genomeFiles = [join(args.genomes,f) for f in listdir(args.genomes) if isfile(join(args.genomes,f))]
    listOfSpecies = makeSpeciesList([args.referenceSpecies], args.genomes)
    
    geneObjects, versionInfo = infernal(args.outPutDir, genomeFiles, args.model, args.incE, args.incT, args.infernalPath)
    #print(len(geneObjects))
    #print("chromo: {}".format(geneObjects[0].chromosome))
    #print("species: {}".format(geneObjects[0].species))
    
elif args.own_genes != None:
    #print("own Genes")

    geneObjects, listOfSpecies = parseGenes(args.own_genes, [args.referenceSpecies])
    #print("search Genes")
    numgenes = len(geneObjects)
    numspecs = len(listOfSpecies)
    #print("num genes: {}".format(numgenes))
    #print("num species: {}".format(numspecs))
    versionInfo = "User provided genes"
    
else:
    raise Exception("either own_genes or search_genes must be given")


#print("num species: {}".format(len(geneObjects)))
'''
g = 0
for _list in geneObjects:
    g += len(_list)
    print(_list[0])
print("num genes: {}".format(g))
'''
maf2bed(args.multiSeq, args.outPutDir, args.pathToRepo, geneObjects, listOfSpecies, args.quality, versionInfo)

