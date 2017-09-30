![SMORE](http://www.bioinf.uni-leipzig.de/~bsarah/smore_logo3.png "")

# SMORE - Synteny Modulator of Repetitive Elements.

This is the readme to the SMORE pipeline, Synteny Modulator Of Repetitive Elements.

Usage: ./smore [general_options]

The general options are:
        --help|-h       print help page
	--version|-v    print version information
	--citation	print citation information
	--contact       print contact information
        --helpout       print explanation for structure and visualization of
	                output files

For subcommands, use: ./smore <subcommand> [options]

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



For more details, please see the SMORE manual or the SMORE help pages.

For subcommands' help pages:
./smore <subcommand> --help

