This is the readme to the SMORE pipeline, a Synteny Modulator Of Repetitive Elements.
usage: ./smore <subcommand> [options]

For more details, please see the SMORE manual.

The smore subcommands are:

    bake    The subcommand bake combines the subcommands prep and toast in
            order to easily start and run the pipeline completely. Hence,
	    parameters for bake are the same as the combined parameters
	    for prep and toast.

    prep    This program will sort genetic elements in between genomic anchors
            based on MultiZ alignments. The genetic elements are
            taken from a list given as input or retrieved based on a covariance
	    model as input for infernal.

    toast   This program will take the prep-output and calculate the numbers
            for evolutionary events at the given phylogenetic tree.
            For a more detailed output, use --verbose or any combination
            of --clus, --graph, --aln. Please note that the verbose version
	    might take significantly longer.

    mix     This subcommand can be used after running SMORE prep. SMORE mix
            only produces a list of genetic clusters. This can be used to
	    test different joining methods. Additionally, SMORE mix can split
	    the number of clusters in several disjoint lists such that the
	    succeeding subcommands can be run in parallel on the disjoint
	    subsets to speed up running time of the program for large
	    data sets. In case SMORE mix is splitting the list of clusters,
	    it will output a command list that can be called to continue
	    running the pipeline. The next subcommand to be called is
	    SMORE roast.

    roast   SMORE roast starts from a list of genetic clusters and will output
            lists of genetic events that can be further proceeded as an input
	    for SMORE eat. SMORE roast can be used for large data sets in order
	    to split and parallelize the process. The next subcommand in the
	    pipeline will be SMORE eat.

    eat     This subcommand is used after SMORE roast and will take list(s) of
            evolutionary events and output a phylogenetic tree with event and
	    element counts. SMORE eat is able to summarize the outputs of
	    parallelized runs such that all data is recombined at the
	    phylogenetic tree.


Subcommands and their options:

    #general options:
      --help|-h         print help page
      --version|-v      print version information
      --citation	print citation information
      --contact         print contact information
      --helpout         print explanation for structure and visualization of
	                output files

    #options needed in all subcommands
      ##obligatory parameter
      --out|-o PATH     folder where to write the output files. If the folder
	                does not exist, it will be created
      ##optional parameter
      --tool|-t PATH    path to where the smore tool is located. This is only
	                needed, if the tools is in a special location.
      --python PATH     path to python, only needed when the python path is not
	                located in the environment variable
      --perl PATH       path to perl, only needed when the perl path is not
	                located in the environment variable

    #options for prep
      ##General options
      --ref|-r REF      name of the reference species as given in MultiZ files
      --maf|-m PATH     folder containing all MultiZ alignments. Files can be
	                compressed.
      --filter NUM      optional: remove a percentage of the lowest scoring
	                blocks base on the MAF scores; number between 0 and 100
      ##CM mode
      --cm|-c FILE      covariance model file, input for infernal
      --genomes|-g PATH folder with genomes of the species used as an input to
	                infernal to scan the genomes for genetic elements
			specified by the CM. Filenames should match species
			names in MultiZ files.
      --incE NUM        optional parameter for infernal, e-value threshold
      --incT NUM        optional parameter for infernal, bitscore threshold
      --infernal PATH   optionel: path to infernal, only needed when the
			cmsearch path is not located in the environment variable
      --pseudo NUM      optional parameter, bitscore threshold that defines
	                pseudogenes
      ##Genelist mode
      --genes FILE      a list of genetic elements as input.
			Format (tab separated):
			
			chromosome start_pos end_pos species strand type
			pseudogene sequence structure comment
			
			-fields to be omitted are filled with NA
			-species name should be the same as in MultiZ files
			-type is optional (e.g. tRNA type Met)
			-pseudogene (optional): true or false
			-structure is optional, dot-bracket-notation
			-for optional sequence, take loci mode
			-comment is optional
      ##Loci mode
      --genomes|-g PATH folder with genomes of the species used to retrieve
			the sequence of the genetic elements in the loci list
			based on their coordinates. Filenames should match
			species names in loci list files.
      --loci FILE       list with genetic elements without sequence
			Format (tab separated):
			
			chromosome start_pos end_pos species strand type
			pseudogene
			
			-fields to be omitted are filled with NA
			-species name should be the same as in MultiZ files
			-type is optional (e.g. tRNA type Met)
			-pseudogene (optional): true or false

      ##tRNAscan-SE mode
      --genomes|-g PATH folder with genomes of the species used to retrieve
			the sequence of the genetic elements in the loci list
			based on their coordinates. Filenames should match
			species names in loci list files.
      --trna            option to turn on the trna mode, hence the program will
			search for tRNAs in the provided genomes using tRNAscan-SE
      --trnascan PATH   optional, path to tRNAscan-SE installation. Only needed
			if the path to tRNAscan-SE is not stored in the environment


    #options for toast
      --prep PATH       folder with output from smore prep run
      --seqsim|-s NUM   percentage of sequence similarity to be considered
	                homolog sequences, default 0.8
      --strucsim|-p NUM percentage of structure similarity to be considered
	                homolog sequences, default 0.8
      --newick FILE     tree in newick format containing the species included
	                in the output of smore prep.
      --id FILE		in case the species' identifier do not fit, this file
	 		can be used to automatically translate the names
			Format (tab separated):
			current_name_in_tree name_to_translate_to
      --join OPT	the way of how original clusters should be joined,
			either none, strict or relaxed; default: relaxed
      --nomiss          do not check for missing anchors. in case of large
                        clusters, using --nomiss will speed up the runtime
      --verbose         print all intermediary files, thus clusters, graphs
                        and duplication alignments
      --clus            print only intermediary cluster files
      --graph           print only intermediary graph files
      --aln             print only intermediary duplication alignment files
			

    #options for mix
      --prep PATH       folder with output from smore prep run
      --species FILE    file listing all species one per line. Species
                        identifier must be the same as in newick or id
			and prep output data.
      --seqsim|-s NUM   percentage of sequence similarity to be considered
	                homolog sequences, default 0.8
      --strucsim|-p NUM percentage of structure similarity to be considered
	                homolog sequences, default 0.8
      --newick FILE     tree in newick format containing the species included
	                in the output of smore prep.
      --id FILE		in case the species' identifier do not fit, this file
	 		can be used to automatically translate the names
			Format (tab separated):
			current_name_in_tree name_to_translate_to
      --max NUM         maximal number of clusters to be included in the
                        following steps of the analysis. In case there exist
			more clusters, the program automatically splits the
			data set and creates a command list for the following
			steps of the pipeline. Default: 50000.
      --join LEVEL	the way of how original clusters should be joined,
			either none, strict or relaxed; default: relaxed
      --nomiss          do not check for missing anchors. in case of large
                        clusters, using --nomiss will speed up the runtime
      --verbose         print all intermediary files, thus clusters
      --clus            print only intermediary cluster files
      
    #options for roast
      --in|-i PATH      output file of smore mix containing a list of
                        clusters
      --species FILE    file listing all species one per line. Species
                        identifier must be the same as in newick or id
			and prep output data.
      --seqsim|-s NUM   percentage of sequence similarity to be considered
	                homolog sequences, default 0.8
      --strucsim|-p NUM percentage of structure similarity to be considered
	                homolog sequences, default 0.8
      --newick FILE     tree in newick format containing the species included
	                in the output of smore prep.
      --id FILE		in case the species' identifier do not fit, this file
	 		can be used to automatically translate the names
			Format (tab separated):
			current_name_in_tree name_to_translate_to
      --nomiss          do not check for missing anchors. in case of large
                        clusters, using --nomiss will speed up the runtime
      --verbose         print all intermediary files,thus graphs
                        and duplication alignments
      --graph           print only intermediary graph files
      --aln             print only intermediary duplication alignment files

    #options for eat
      --prep PATH       folder with output from smore prep run
      --mix PATH        folder with output from smore mix run
      --roast PATH      folder with output from smore roast run
      --newick FILE     tree in newick format containing the species included
	                in the output of smore prep.
      --id FILE		in case the species' identifier do not fit, this file
	 		can be used to automatically translate the names
			Format (tab separated):
			current_name_in_tree name_to_translate_to
      --nomiss          do not check for missing anchors. in case of large
                        clusters, using --nomiss will speed up the runtime



#General Output:
    data_iTOL:   folder giving files that can be uploaded to itol.embl.de
                 The file called F0tree.txt is uploaded at
                 http://itol.embl.de/upload.cgi.
                 The remaining files can be added using drag and drop
                 into the browser window. This will result in a
                 interactive visualization of the resulting tree.
		 A legend is added automatically. All nodes in the tree will
                 have unique names, thus some nodes might have names such as
                 'innerNode0' because it was added automatically.

    OutTree.txt: the resulting tree in newick format with numbers at the
                 nodes given in brackets. This format can be used to
                 visualize the tree with newick compatible programs.

geneticEvents.txt: File listing all genetic events counted during the 
                   analysis. The numbers are sorted by event and node of
                   the tree. The file includes a event called 'Other'. This
                   will give the difference of genetic elements between the
                   total amount and the elements used in the analysis.
                   For a successful run of the pipeline, the numbers should
                   be 0.

allClusters_original.txt 
and allClusters_joined.txt: These two files contain lists of clusters
                            showing which elements are contained together
                            in one cluster before and after joining.

*_errors: for each part of the smore pipeline, there is a file giving errors
          that happened during the run. If no errors occured, the file is
          empty.

list_cographs.txt  and
list_noncographs.txt: these files contain statistics about graphs that were
                      cographs from the beginning or had to be edited in
                      order to become a cograph. The tables list number of
                      nodes, number of edges and density of the graphs.

remoldings.txt and
inremoldings.txt: These files contain genetic elements that (a) have
                  highly similar sequences but different types or (b)
                  have the same types but clearly distinct sequences.
						

allTypes.txt and
allPseudoTypes.txt: These two files list the different types of genetic
                    elements for all species and functional or pseudogenized
                    genes. This can be used to analyse the distribution of
                    different types of genetic elements.

Additional files:
- for each species, there are files listing their genetic elements as
they were used as intermediate files.
- there are files for each of the evolutionary events where counts are
listed.
- for singletons, there is also a listing about types and pseudogenes
- for elements which could not be sorted in between genomic anchors (nones)
there is listings about types and pseudogenes.



Verbose output: the verbose version of toast will output three additional
files (if --verbose) or any combination of them (if --clus, --graph, --aln)
for each cluster, named with left and right anchor numbers to match all three
files. They will be in three different folders: cluster, graph and
duplication_alignment. The files contain the specific structures of the
cluster in each step of the analysis and can be used to gain a deeper
insight into the data.