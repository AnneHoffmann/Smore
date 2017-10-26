#!/usr/bin/perl -w

#smoretoast takes smoreprep output, generates cluster, joins them, created graphs, checks the graphs for cograph structure, creates alignments and counts events
#intermediary files are not printed but some statistics will be included in the summary
#


use IO::Uncompress::Gunzip qw($GunzipError);
use Data::Dumper;
use Getopt::Long qw(GetOptions);
use strict;
use warnings;


##options for all modes
my $toolpath;
my $outpath;
my $pythonpath;
my $perlpath;
##options for toast and bake
my $pathtocam; #path to output of smore prep (genes folder)
my $seqsim=0.8;
my $strucsim=0.8;
my $newicktree;
my $joinmode="relaxed";
my $specieslist;
my $nocheck;
my $errors;
##options for verbose
my $printall;
my $printg;
my $printa;
my $printc;


GetOptions(
    #all modes
    'tool|t=s' => \$toolpath,
    'out|o=s' => \$outpath,
    'python=s' => \$pythonpath,
    'perl=s' => \$perlpath,
##options for toast and bake
    'prep=s' => \$pathtocam, #path to output of smore prep (genes folder)
    'seqsim|s=f' => \$seqsim,
    'strucsim|p=f' => \$strucsim,
    'newick=s' => \$newicktree,
    'join=s' => \$joinmode,
    'species=s' => \$specieslist,
    'nomiss' => \$nocheck,
    'err=s' => \$errors,
##options for verbose
    'verbose' => \$printall,
    'graph' => \$printg,
    'aln' => \$printa,
    'clus' => \$printc 
    ) or die "error in smoretoast";


my $mode = 1;

my $tmpfilelist = "$outpath\/tmpfiles";
my $tmplscmd = "ls $pathtocam\/temp \> $tmpfilelist";
my @outlscmd = readpipe("$tmplscmd");

my $del2check = "$outpath\/tmpfile\_deletionsCheck\.txt";
my $pseudel2check = "$outpath\/tmpfile\_pseudeletionsCheck\.txt";


my %blocks = ();
my %leftblocks = ();
my %species = ();
my %pseuspecies = ();
my %nones = ();
my %pseudonones = ();
my %nonesNT = ();
my %pseudononesNT = ();

my $remlist = "$outpath\/remoldings\.txt";
my $inremlist = "$outpath\/inremoldings\.txt";
my $remcmd = "touch $remlist";
my $inremcmd = "touch $inremlist";
readpipe("$remcmd");
readpipe("$inremcmd");

open(my $outsr,">>", $remlist);
open(my $outsi,">>", $inremlist);


my %allTypes = (); #hash with spec_type -> num
my %pallTypes = (); #hash with spec_type -> num

my %types = (); #this is only the types, not the species

my %block2spec = ();
my @allspec = ();
my %neighbors = (); #species_anchor -> leftneighbor_rightneighbor or -1 

my %spec2mindist = ();  #in order to find the max/min distance between blocks in species
my %spec2maxdist = ();
my %spec2sumdist = ();
my %spec2numdist = ();

open FA,"<$specieslist" or die "can't open $specieslist\n";

while(<FA>){
    chomp;
    my $precurfile = $_;
    my $curfile = "$pathtocam\/genes\/$precurfile";
    #get species name
    my @S1 = split '\/', $curfile;
    my @S2 = split '\.', $S1[-1];
    my $spec = $S2[0];
    push @allspec, $spec;
    #insert in %specis
    $species{$spec} = 0;
    $spec2maxdist{$spec} = 0;
    $spec2mindist{$spec} = 1000000;
    $spec2sumdist{$spec} = 0;
    $spec2numdist{$spec} = 0;
    $pseuspecies{$spec} = 0;
    #sort file and remove blanks
#    my $curfile_nb = "$outpath\/$spec\_nb\.bed";
#    my $curfile_sorted = "$outpath\/$spec\_sorted\.bed";
#    my $curfile_sorteduniq = "$outpath\/$spec\_sorteduniq\.bed";
#    my $cmdnoblanks = "awk \'NF\' $curfile \> $curfile_nb";
    #sort by chr, elemstart, elemend
#    my $col1 = 1;
#    my $col3 = 3;
#    my $col4 = 4;
#    my $coln = "n";
#    my $cmdsort = "sort -k$col1,$col1 -k$col3,$col3$coln -k$col4,$col4$coln $curfile_nb \> $curfile_sorted";
#    my $cmduniq = "uniq $curfile_sorted > $curfile_sorteduniq";
#    readpipe("$cmdnoblanks");
#    readpipe("$cmdsort");
#    readpipe("$cmduniq");
    #this describes the neighbors
    my $chr2 = "-1";
    my $bl2 = "-1";
    my $chr1 = "-1";
    my $bl1 = "-1";
    #read file to blocks/nones
    open CF, "<$curfile" or die "can't open $curfile\n";
    while(<CF>){
	chomp;
	my $line = $_;
	if($line=~/^#/){next;}
	if($line eq ""){next;}
	my @F = split '\t', $line;

	if($mode==1)
	{
	    my $tkey = "$spec\_$F[-3]";
	    if(exists($types{$F[-3]})){$types{$F[-3]}++;}
	    else{$types{$F[-3]}=1;}
	    if($F[-2] eq "t" || $F[-2] eq "T" ||$F[-2] eq "TRUE" || $F[-2] eq "True" || $F[-2] eq "true" || $F[-2] eq "1"){
		$pseuspecies{$spec}++;
		##add to allTypes hash
		if(exists($pallTypes{$tkey})){$pallTypes{$tkey}++;}else{$pallTypes{$tkey}=1;}
	    }
	    else{
		$species{$spec}++;
		##add to allTypes hash
		if(exists($allTypes{$tkey})){$allTypes{$tkey}++;}else{$allTypes{$tkey}=1;}
	    }
	}
	else{
	    $species{$spec}++;
	}
	my $none = "None";
	if($F[5] eq $none || $F[6] eq $none){
	    my $nonekey;
	    if($mode==1){
		my $type = $F[-3];
		$nonekey = "$spec\_$type";
		my $pseu = $F[-2];
		if($pseu eq "TRUE" || $pseu eq "True" || $pseu eq "true" || $pseu eq "1"){
		    if(exists($pseudonones{$nonekey})){
			$pseudonones{$nonekey}++;
		    }
		    else{
			$pseudonones{$nonekey}=1;
		    }
		    if(exists($pseudononesNT{$spec})){
			$pseudononesNT{$spec}++;
		    }
		    else{
			$pseudononesNT{$spec}=1;
		    }
		}
		else{
		    if(exists($nones{$nonekey})){
			$nones{$nonekey}++;
		    }
		    else{
			$nones{$nonekey}=1;
		    }
		    if(exists($nonesNT{$spec})){
			$nonesNT{$spec}++;
		    }
		    else{
			$nonesNT{$spec}=1;
		    }
		}
	    }
	    else{
		$nonekey = "$spec";
		if(exists($nones{$nonekey})){
		    $nones{$nonekey}++;
		}
		else{
		    $nones{$nonekey}=1;
		}
		if(exists($nonesNT{$spec})){
		    $nonesNT{$spec}++;
		}
		else{
		    $nonesNT{$spec}=1;
		}
	    }
	    next;
	}
	else{
	    ##store neighbor relations of blocks in each species
	    my $key;
	    my $leftkey;
	    my $rightkey;
	    if($F[5] < $F[6]){
		$key = "$F[5]\_$F[6]";
		$leftkey = "$F[5]";
		$rightkey = "$F[6]";
	    }
	    else{
		$key = "$F[6]\_$F[5]";
		$leftkey = "$F[6]";
		$rightkey = "$F[5]";
	    }
	    if(exists($blocks{$key})){
		my $tmp = $blocks{$key};
		$blocks{$key} = "$tmp\;$line";
	    }
	    else{
		$blocks{$key} = "$line";
	    }
	    if(exists($leftblocks{$leftkey})){
		my $ltmp = $leftblocks{$leftkey};
		$leftblocks{$leftkey} = "$ltmp\;$line";
	    }
	    else{
		$leftblocks{$leftkey} = "$line";
	    }
	    if($chr1 eq "-1"){
		$chr1 = $F[1];
		$bl1 = $leftkey;
		my $nkey = "$spec\_$bl1";
		my $nval = "$bl2\_$rightkey";
		$neighbors{$nkey} = $nval;
		$chr2 = $chr1;
		$bl2 = $bl1;
		$chr1 = $F[1];
		$bl1 = $rightkey;

	    }
	    elsif($chr1 eq $F[1]){
		if($bl1 ne $leftkey){
		    my $nkey = "$spec\_$bl1";
		    my $nval = "$bl2\_$leftkey";
		    $neighbors{$nkey} = $nval;
		    my $mkey = "$spec\_$leftkey";
		    my $mval = "$bl1\_$rightkey";
		    $neighbors{$mkey} = $mval;
		    $chr2 = $chr1;
		    $bl2 = $bl1;
		    $chr1 = $F[1];
		    $bl1 = $rightkey;
		}
		else{
		    my $mkey = "$spec\_$leftkey";
		    my $mval = "$bl2\_$rightkey";
		    $neighbors{$mkey} = $mval;
		    $chr2 = $chr1;
		    $bl2 = $bl1;
		    $chr1 = $F[1];
		    $bl1 = $rightkey;		    
		}
	    }
	    elsif($chr1 ne $F[1]){
		my $nkey = "$spec\_$bl1";
		my $nval = "$bl2\_-1";
		$neighbors{$nkey} = $nval;
		my $mkey = "$spec\_$leftkey";
		my $mval = "-1\_$rightkey";
		$neighbors{$mkey} = $mval;
		$bl2 = $leftkey;
		$chr2 = $F[1];
		$bl1 = $rightkey;
		$chr1 = $F[1];
	    }
	    else{
	    }
	}
    }
}
my $numspec = scalar @allspec;
print "Number of species: $numspec \n";

my $typenum = scalar (keys %types);
print "Number of different element types: $typenum \n";

my $orignum = scalar (keys %blocks);
print "Number of original clusters: $orignum \n";

##analyse original blocks for distances between elements
foreach my $block (keys %blocks){
    #output a hash with distance for each species that exist in the cluster and have a max and min pos
    my %dists = getMaxDist($blocks{$block});
    foreach my $s (keys %dists){
	if($dists{$s} == 0){next;}
	if(exists($spec2mindist{$s})){if($dists{$s} < $spec2mindist{$s}){$spec2mindist{$s}=$dists{$s};}}
	else{$spec2mindist{$s}=$dists{$s};}
	if(exists($spec2maxdist{$s})){if($dists{$s} > $spec2maxdist{$s}){$spec2maxdist{$s}=$dists{$s};}}
	else{$spec2maxdist{$s}=$dists{$s};}
	if(exists($spec2sumdist{$s})){$spec2sumdist{$s}+=$dists{$s};}
	else{$spec2sumdist{$s}=$dists{$s};}
	if(exists($spec2numdist{$s})){$spec2numdist{$s}++;}
	else{$spec2numdist{$s}=1;}
    }
}

print "The following numbers show the maximal, minimal and average distances between
elements in original clusters for each species (counted in nt).\n";
print "species\tmaxdist\tmindist\taverage\tnum_of_counted_clusters\n";
#calculate average distances
foreach my $t (keys %spec2numdist){
    my $avdist = 0;
    if($spec2numdist{$t}>0){
	$avdist = sprintf "%.2f",$spec2sumdist{$t}/$spec2numdist{$t};
    }
    if($spec2mindist{$t}==1000000){$spec2mindist{$t} = 0;}
    print "$t\t$spec2maxdist{$t}\t$spec2mindist{$t}\t$avdist\t$spec2numdist{$t}\n";
    #reset the hashes in order to use it for joined clusters
    $spec2maxdist{$t} = 0;
    $spec2mindist{$t} = 1000000;
    $spec2sumdist{$t} = 0;
    $spec2numdist{$t} = 0;
}





#join cluster by hash num
my %joinedblocks = ();

my $curblock = "";
my $curstart = -5;#a
my $curend = -1;#a

my @spezi = keys %species;
my $joined = 0;

if($joinmode eq "relaxed" || $joinmode eq "strict"){
    $joined = 1;
    foreach my $lb (sort { $a <=> $b} keys %leftblocks){
	my $curb = $leftblocks{$lb};
	my $rb = getRightBlock($curb); #this is the right block, -1 if there are several rightblocks, then it doesnt work 
	if($curstart == -5){
	    $curblock = $curb;
	    $curstart = $lb;
	    $curend = $rb;
	    next;
	}
	#as keys are sorted, if the blocks are adjacent and joinable, lb should be the same as curend
	#adjacent doesn't have to be directly adjacent but all elements in both clusters should be adjacent in all species
	my $isgood = 1;
	#change checkBlocks based on the created blocktable
	#grep for blocknum in each sorted file (for each species) and check if the blocknums are adjacent in all of them
	if($rb>0 && $curend > 0){
	    if($joinmode eq "relaxed"){
		my $futureblocks = "$curblock\;$curb";
		@spezi = getSpecList($futureblocks);
	    }
	    for(my $sp =0;$sp < scalar @spezi;$sp++){#in this line, we could only take species that are really involved in the current clusters
		my $curlb = "$spezi[$sp]\_$lb";
		if(exists($neighbors{$curlb})){
		    my $spn = $neighbors{$curlb};
		    my @Q = split '_', $spn;
		    if($Q[1] eq $rb || $Q[1] == $rb){
			$isgood=1;
		    }
		    else{$isgood=0;last;}
		}
		else{$isgood=0;last;}
	    }	
	}
	if($rb>0 && $isgood == 1 && $curend > 0){
	    #should be joinable
	    $curend = $rb;
	    $curblock = "$curblock\;$curb";
	}
	else{
	    #add the old block to joined clusters
	    if($curend == -1){
		#check if the elements are still neighbors
		#if not separately add each element in curblock
		my %rb2clus = ();
		my @F = split ';', $curblock;
		for(my $f=0;$f<scalar @F;$f++){
		    my @G = split '\t', $F[$f];
		    my $rkey;
		    if($G[5] < $G[6]){$rkey = $G[6]}
		    else{$rkey = $G[5]}
		    if(exists($rb2clus{$rkey})){$rb2clus{$rkey}="$rb2clus{$rkey}\;$F[$f]";}		  
		    else{$rb2clus{$rkey} = $F[$f];}
		}
		#curstart = $lb
		my $curclus = "";
		my $curend = -1;
		my $neighborsum = 0; #if this is equals scalar @rks, then we don't have neighbors
		my @rks = sort { $a <=> $b} keys %rb2clus;
		for(my $r=0;$r<scalar @rks;$r++){
		    my $noneighbor = 0;
		    for(my $k=$r+1;$k<scalar @rks;$k++){
			for(my $sp =0;$sp < scalar @spezi;$sp++){
			    my $curlb = "$spezi[$sp]\_$rks[$r]";
			    if(exists($neighbors{$curlb})){
				my $spn = $neighbors{$curlb};
				my @Q = split '_', $spn;
				if($Q[1] eq $rb || $Q[1] == $rb){
				}
				else{
				    $neighborsum += 1;
				    last;
				}
			    }
			    else{
				$neighborsum += 1;
				last;
			    }
			}	
		    }
		}
		
		if($neighborsum < scalar @rks){ ##check again, but we first check the prints
		    my $nkey = "$lb\_$rks[-1]";
		    $joinedblocks{$nkey} = join(';',values %rb2clus);
		}
		else{
		    for(my $e=0;$e<scalar @rks;$e++){
			my $newkey = "$lb\_$rks[$e]";
			if(exists($joinedblocks{$newkey})){
			    $joinedblocks{$newkey}="$joinedblocks{$newkey}\;$rb2clus{$rks[$e]}";
			}
			else{$joinedblocks{$newkey}="$rb2clus{$rks[$e]}";}
		    }
		}
	    }
	    else{
		#all elements in curblock are in one cluster
		my $newkey2 = "$curstart\_$curend";
		
		if(exists($joinedblocks{$newkey2})){
		    $joinedblocks{$newkey2}="$joinedblocks{$newkey2}\;$curblock";
		}
		else{$joinedblocks{$newkey2}=$curblock;}
	    }
	    #set the new block
	    $curblock = $curb;
	    $curstart = $lb;
	    $curend = $rb;	
	}
	
    }
    #add the last cluster to hash
    if($curstart > 0){
	if($curend == -1){
	    my @F = split ';', $curblock;
	    for(my $f=0;$f<scalar @F;$f++){
		my @G = split '\t', $F[$f];
		my $rkey;
		if($G[5] < $G[6]){$rkey = $G[6]}
		else{$rkey = $G[5]}
		my $nowkey = "$curstart\_$rkey";
		if(exists($joinedblocks{$nowkey})){$joinedblocks{$nowkey}="$joinedblocks{$nowkey}\;$F[$f]";}		  
		else{$joinedblocks{$nowkey} = $F[$f];}
	    }
	}
	else{
	    my $newkey2 = "$curstart\_$curend";
	    if(exists($joinedblocks{$newkey2})){
		$joinedblocks{$newkey2}="$joinedblocks{$newkey2}\;$curblock";
	    }
	    else{$joinedblocks{$newkey2}=$curblock;}
	}
    }
}

##after joining or no joining

my $sumelemorig = 0;
my $jclusfile = "$outpath\/allClusters_original\.txt";
open(my $outjc,">>",$jclusfile);
my $origcount = 0;
foreach my $bk (keys %blocks){
    my @B = split ';', $blocks{$bk};
    my $clussize = scalar @B;
    #write to allClustersList
    print $outjc ">clus$origcount $clussize\n";
    $origcount++;
    $sumelemorig += $clussize;
    for(my $cb = 0;$cb < $clussize;$cb++){
	print $outjc "$B[$cb]\n";
    }
}
close $outjc;
my $origavsize = 0;
if($orignum > 0){
    $origavsize = sprintf "%.2f",$sumelemorig/$orignum;
}
print "Average number of elements of original clusters: $origavsize \n";

if($joined == 1){    
    my $newclusnum = scalar (keys %joinedblocks);
    print "Number of joined clusters: $newclusnum \n";
    %blocks = %joinedblocks;

    ##analyse original blocks for distances between elements
    foreach my $block (keys %blocks){
	#output a hash with distance for each species that exist in the cluster and have a max and min pos
	my %dists = getMaxDist($blocks{$block});
	foreach my $s (keys %dists){
	    if($dists{$s} == 0){next;}
	    if(exists($spec2mindist{$s})){if($dists{$s} < $spec2mindist{$s}){$spec2mindist{$s}=$dists{$s};}}
	    else{$spec2mindist{$s}=$dists{$s};}
	    if(exists($spec2maxdist{$s})){if($dists{$s} > $spec2maxdist{$s}){$spec2maxdist{$s}=$dists{$s};}}
	    else{$spec2maxdist{$s}=$dists{$s};}
	    if(exists($spec2sumdist{$s})){$spec2sumdist{$s}+=$dists{$s};}
	    else{$spec2sumdist{$s}=$dists{$s};}
	    if(exists($spec2numdist{$s})){$spec2numdist{$s}++;}
	    else{$spec2numdist{$s}=1;}
	}
    }
    
    print "The following numbers show the maximal, minimal and average distances between
elements in joined clusters for each species (counted in nt).\n";
    print "species\tmaxdist\tmindist\taverage\tnum_of_counted_clusters\n";
    #calculate average distances
    foreach my $t (keys %spec2numdist){
	my $avdist = 0;
	if($spec2numdist{$t}>0){
	    $avdist = sprintf "%.2f",$spec2sumdist{$t}/$spec2numdist{$t};
	}
	if($spec2mindist{$t}==1000000){$spec2mindist{$t} = 0;}
	print "$t\t$spec2maxdist{$t}\t$spec2mindist{$t}\t$avdist\t$spec2numdist{$t}\n";
	#reset the hashes 
	$spec2maxdist{$t} = 0;
	$spec2mindist{$t} = 1000000;
	$spec2sumdist{$t} = 0;
	$spec2numdist{$t} = 0;
    }
}

##print None hash to file
my $noneout = "$outpath\/nones.txt";
my $pseunoneout = "$outpath\/pseunones.txt";
my $nonestr = "";
my $pseudononestr = "";
foreach my $k (keys %nones){
    $nonestr = "$nonestr$k\t$nones{$k}\n";
}
foreach my $k (keys %pseudonones){
    $pseudononestr = "$pseudononestr$k\t$pseudonones{$k}\n";
}
open(my $outnone,">>",$noneout);
print $outnone $nonestr;
close $outnone;
open(my $outpnone, ">>", $pseunoneout);
print $outpnone $pseudononestr;
close $outpnone;

#print allTypes to file
my $typesout = "$outpath\/allTypes\.txt";
my $ptypesout = "$outpath\/allPseudoTypes\.txt";
open(my $tyout, ">>", $typesout);
foreach my $ty (keys %allTypes){
    print $tyout "$ty\t$allTypes{$ty}\n";
}
close $tyout;
open(my $ptyout, ">>", $ptypesout);
foreach my $pty (keys %pallTypes){
    print $ptyout "$pty\t$pallTypes{$pty}\n";
}
close $ptyout;


my $allnonestr = "";
foreach my $k (keys %nonesNT){
    $allnonestr = "$allnonestr$k\-$nonesNT{$k}\=";
}
if(scalar (keys %nonesNT) == 0){$allnonestr = "=";}
$allnonestr = "$allnonestr\!";
foreach my $k (keys %pseudononesNT){
    $allnonestr = "$allnonestr$k\-$pseudononesNT{$k}\=";
}

my $allspecstr = "";
foreach my $t (keys %species){
    $allspecstr = "$allspecstr$t\-$species{$t}\=";
}
$allspecstr = "$allspecstr\!";
foreach my $t (keys %pseuspecies){
    $allspecstr = "$allspecstr$t\-$pseuspecies{$t}\=";
}


##verbose options
my $clusfolder = "$outpath\/cluster";
if($printc || $printall){
    if(!(-e $clusfolder)){
	my $mkcmd = "mkdir $clusfolder";
	readpipe("$mkcmd");
    }
}
my $graphfolder = "$outpath\/graph";
if($printg || $printall){
    if(!(-e $graphfolder)){
	my $mkcmd2 = "mkdir $graphfolder";
	readpipe("$mkcmd2");
    }
}
my $alnfolder = "$outpath\/duplication_alignments";
if($printa || $printall){
    if(!(-e $alnfolder)){
	my $mkcmd3 = "mkdir $alnfolder";
	readpipe("$mkcmd3");
    }
}




my %dupevents =();
my %matevents =();
my %insevents =();
my %misevents =(); #missing data=no anchors
my %delevents = ();#deletions
my %psevents = (); #insertions of pseudogenes
my %psedels = ();
my %psemis = ();
my %pseins = ();
my %singletons = (); #hash species -> number
my %pseusingles = (); #include types
my %singletonsNT = (); #hash species -> number
my %pseusinglesNT = (); #include types

##go on with cluster hash
##sort keys of blocks and join adjacent ones?

my $cglist = "$outpath\/list\_cographs\.txt";
my $ncglist = "$outpath\/list\_noncographs\.txt";
my $cgcmd = "touch $cglist";
my $ncgcmd = "touch $ncglist";
readpipe("$cgcmd");
readpipe("$ncgcmd");


open(my $outcg,">>",$cglist);
open(my $outn,">>",$ncglist);
my $graphheader = "graph \t nodenum \t edgenum \t corr. edgenum \t density\n";
print $outcg $graphheader;
print $outn $graphheader;
close $outcg;
close $outn;


my $eventlist = "$outpath\/allClusters_joined\.txt";
my $listcmd = "touch $eventlist";
readpipe("$listcmd");
open(my $outg, ">>",$eventlist);
my $cluscount = 0;
my $sumelems = 0;
foreach my $k (keys %blocks){
    my @A = split '_', $k;
    my $leftanchor = $A[0];
    my $rightanchor = $A[1];
    my @B = split ';', $blocks{$k};
    my $clussize = scalar @B;
    #verbose option for cluster
    my $outcf;
    if($printc || $printall){
	my $clusfile = "$clusfolder\/cluster$leftanchor\_$rightanchor\.clus";
	open($outcf, ">>", $clusfile);
    }
    #write to allClustersList
    if($joinmode eq "relaxed" || $joinmode eq "strict"){
	print $outg ">clus$cluscount $clussize\n";
	$cluscount++;
	$sumelems+=$clussize;
	for(my $cb = 0;$cb < $clussize;$cb++){
	    print $outg "$B[$cb]\n";
	    if($printc || $printall){print $outcf "$B[$cb]\n";}
	}
    }
    if(scalar @B == 1){
	##sort single element cluster into hash
	my @F = split '\t', $B[0];
	my @C = split '_', $F[1];
	my $curspec = $C[0];
	my $singlekey;
	if($mode==1){
	    my $type = $F[-3];
	    $singlekey = "$curspec\_$type";
	    my $pseu = $F[-2];
	    if($pseu eq "T" || $pseu eq "TRUE" || $pseu eq "True" || $pseu eq "true" || $pseu eq "1"){
		if(exists($pseusingles{$singlekey})){
		    $pseusingles{$singlekey}++;
		}
		else{
		    $pseusingles{$singlekey}=1;
		}
		if(exists($pseusinglesNT{$curspec})){
		    $pseusinglesNT{$curspec}++;
		}
		else{
		    $pseusinglesNT{$curspec}=1;
		}
	    }
	    else{
		if(exists($singletons{$singlekey})){
		    $singletons{$singlekey}++;
		}
		else{
		    $singletons{$singlekey} = 1;
		}
		if(exists($singletonsNT{$curspec})){
		    $singletonsNT{$curspec}++;
		}
		else{
		    $singletonsNT{$curspec} = 1;
		}
		
	    }
	}
	else{
	    $singlekey = "$curspec";
	
	    if(exists($singletons{$singlekey})){
		$singletons{$singlekey}++;
	    }
	    else{
		$singletons{$singlekey} = 1;
	    }
	    if(exists($singletonsNT{$curspec})){
		$singletonsNT{$curspec}++;
	    }
	    else{
		$singletonsNT{$curspec} = 1;
	    }
	}
	next;
    }
    #have a tmpfile for blocks string
    my $tmpfile0 = "$outpath\/tmpdata0";
    open(my $outtmp0, ">>", $tmpfile0);
    print $outtmp0 $blocks{$k};
    close $outtmp0;
    my $tmpfile1;
    if($printg || $printall){$tmpfile1 = "$graphfolder\/graph$k\.gr";}
    else{$tmpfile1 = "$outpath\/tmpdata1";}

    #build graph, tmpfile1 will contain the output
    my $graphcmd = "perl $toolpath\/buildEdgeList_fast.pl $tmpfile0 $mode $toolpath $strucsim $seqsim $tmpfile1 2>> $errors";
    my @graphout = readpipe("$graphcmd");
    my $nodenum = $graphout[0];
    
    #check graph for cograph or not and edit slightly

    my $checkcmd = "perl $toolpath\/checkGraph_fast.pl $tmpfile1  $k $seqsim $strucsim $mode $cglist $ncglist 2>> $errors";
    
    ##argument list is getting too long...create a temporary file

    my @newoutgraph = readpipe("$checkcmd");
    my $newgraphstr = $newoutgraph[0];

    my $toprint = 0;
    if($printa || $printall){$toprint = 1;}
    my $filetoprint = "$alnfolder\/alignment$k\.aln";

    #create alignment    
    my $alncmd = "perl $toolpath\/createAlignments_fast2.pl $tmpfile1 $outpath $toolpath $seqsim $strucsim $mode $typenum $newicktree $leftanchor $rightanchor $del2check $pseudel2check $toprint $filetoprint $errors 2>> $errors";

    my @outaln = readpipe("$alncmd"); #this array contains: dup mat insertion pseins pseudomatch deletion missinganchor missinanchorpseu deletionpseu rems inrems
    my $rmcmd;
    if($printg || $printall){$rmcmd = "rm $tmpfile0";}
    else{$rmcmd = "rm $tmpfile1 $tmpfile0";}
    readpipe("$rmcmd");
    #the elements in each pos are separated by = and key and value are separated by -
    #starting with =!
    #if there are no counts, the string only consists of =, thats why, we start with 1 in the loops

    my @dups = split '=', $outaln[0];
    for(my $d=1;$d<scalar @dups;$d++){
	my @dt = split '-', $dups[$d];
	if(exists($dupevents{$dt[0]})){$dupevents{$dt[0]}+=$dt[1];}
	else{$dupevents{$dt[0]} = $dt[1];}
    }
    my @mats = split '=', $outaln[1];
    for(my $a=1;$a<scalar @mats;$a++){
	my @at = split '-', $mats[$a];
	if(exists($matevents{$at[0]})){$matevents{$at[0]}+=$at[1];}
	else{$matevents{$at[0]} = $at[1];}	
    }
    my @ins = split '=', $outaln[2];
    for(my $n=1;$n<scalar @ins;$n++){
	my @nt = split '-', $ins[$n];
	if(exists($insevents{$nt[0]})){$insevents{$nt[0]}+=$nt[1];}
	else{$insevents{$nt[0]} = $nt[1];}	
    }
    my @pins = split '=', $outaln[3];
    for(my $pn=1;$pn<scalar @pins;$pn++){
	my @pnt = split '-', $pins[$pn];
	if(exists($pseins{$pnt[0]})){$pseins{$pnt[0]}+=$pnt[1];}
	else{$pseins{$pnt[0]} = $pnt[1];}	
    }
    my @pmats = split '=', $outaln[4];
    for(my $pa=1;$pa<scalar @pmats;$pa++){
	my @pat = split '-', $pmats[$pa];
	if(exists($psevents{$pat[0]})){$psevents{$pat[0]}+=$pat[1];}
	else{$psevents{$pat[0]} = $pat[1];}	
    }
    my @dels = split '=', $outaln[5];
    for(my $dd=1;$dd<scalar @dels;$dd++){
	my @ddt = split '-', $dels[$dd];
	if(exists($delevents{$ddt[0]})){$delevents{$ddt[0]}+=$ddt[1];}
	else{$delevents{$ddt[0]} = $ddt[1];}
    }
    my @mis = split '=', $outaln[6];
    for(my $m=1;$m<scalar @mis;$m++){
	my @ms = split '-', $mis[$m];
	if(exists($misevents{$ms[0]})){$misevents{$ms[0]}+=$ms[1];}
	else{$misevents{$ms[0]} = $ms[1];}
    }
    my @pmis = split '=', $outaln[7];
    for(my $pm=1;$pm<scalar @pmis;$pm++){
	my @pms = split '-', $pmis[$pm];
	if(exists($psemis{$pms[0]})){$psemis{$pms[0]}+=$pms[1];}
	else{$psemis{$pms[0]} = $pms[1];}
    }
    my @pdels = split '=', $outaln[8];
    for(my $pdd=1;$pdd<scalar @pdels;$pdd++){
	my @pddt = split '-', $pdels[$pdd];
	if(exists($psedels{$pddt[0]})){$psedels{$pddt[0]}+=$pddt[1];}
	else{$psedels{$pddt[0]} = $pddt[1];}
    }
    my @currems = split '=', $outaln[9];
    for(my $sr =0;$sr<scalar @currems;$sr++){
	if($currems[$sr] eq "" || $currems[$sr] eq " " || $currems[$sr] eq "\n"){next;}
	print $outsr "$currems[$sr]\n";
    }
    my @curinrems = split '=', $outaln[10];
    for(my $ri =0;$ri<scalar @curinrems;$ri++){
	if($curinrems[$ri] eq "" || $curinrems[$ri] eq " " || $curinrems[$ri] eq "\n"){next;}
	print $outsi "$curinrems[$ri]\n";
    }
    
    
}
close $outsi;
close $outsr;
close $outg;
##print singletons to file and create allsinglestr
my $singleout = "$outpath\/singletons.txt";
my $pseusingleout = "$outpath\/pseusingletons.txt";
my $singlestr = "";
my $psinglestr = "";
foreach my $s (keys %singletons){
    $singlestr = "$singlestr$s\t$singletons{$s}\n";
}
foreach my $s (keys %pseusingles){
    $psinglestr = "$psinglestr$s\t$pseusingles{$s}\n";
}
open(my $outs, ">>", $singleout);
open(my $outt, ">>", $pseusingleout);
print $outs $singlestr;
print $outt $psinglestr;
close $outs;
close $outt;
my $allsinglestr = "";
foreach my $s (keys %singletonsNT){
    $allsinglestr = "$allsinglestr$s\-$singletonsNT{$s}\=";
}
$allsinglestr = "$allsinglestr\!";
foreach my $s (keys %pseusinglesNT){
    $allsinglestr = "$allsinglestr$s\-$pseusinglesNT{$s}\=";
}



###do the check for missing anchors here!
my %tmphash;
if(! $nocheck){    
    ##first, open and read the tmp files
    open TF,"<$tmpfilelist" or die "can't open $tmpfilelist\n";
    while(<TF>){
	chomp;
	my $curtmpfile = $_;
	my @Tmpname = split '_', $curtmpfile;
	my $tmpspec = $Tmpname[0];
	my $CTF = IO::Uncompress::Gunzip->new( "$pathtocam\/temp/$curtmpfile" )
	    or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
	while(<$CTF>){
	    my $tmpline = $_;
	    my @TL = split '\t', $tmpline;
	    my @PTL = split '_', $TL[1];
	    my $tan = $PTL[1];
	    $tmphash{$tmpspec}{$tan} = 1;
	}
    }
}
##now, we check if anchors exist with
open DC,"<$del2check" or die "can't open $del2check \n";
while(<DC>){
    my $l2c = $_;
    my @ll2c = split '\t', $l2c;
    my @anchors = split '_', $ll2c[0];
    my @lspec = split ',',$ll2c[1];
    my $mini = $ll2c[2];
    my @missings = ();
    my @dells = ();
    for(my $l=0;$l< scalar @lspec;$l++){
	if(! $nocheck && exists($tmphash{$lspec[$l]})){
	    if(exists($tmphash{$lspec[$l]}{$anchors[0]})
	       && exists($tmphash{$lspec[$l]}{$anchors[1]})){
		push @dells, $lspec[$l];
	    }
	    else{
		push @missings, $lspec[$l];
	    }
	}
	else{print STDERR "species should exist: $lspec[$l] \n";}
    }
    if(scalar @missings > 0){
	my $misstr = join(',', @missings);
	if(exists $misevents{$misstr}){$misevents{$misstr} += $mini;}
	else{$misevents{$misstr} = $mini;}
    }
    if(scalar @dells > 0){
	my $delstr = join(',',@dells);
	if(exists $delevents{$delstr}){$delevents{$delstr} += $mini;}
	else{$delevents{$delstr} = $mini;}
    }
   
}


open PC,"<$pseudel2check" or die "can't open $pseudel2check \n";
while(<PC>){
    my $l2c = $_;
    my @ll2c = split '\t', $l2c;
    my @anchors = split '_', $ll2c[0];
    my @lspec = split ',',$ll2c[1];
    my $mini = $ll2c[2];
    my @missings = ();
    my @dells = ();
    for(my $l=0;$l< scalar @lspec;$l++){
	if(! $nocheck && exists($tmphash{$lspec[$l]})){
	    if(exists($tmphash{$lspec[$l]}{$anchors[0]})
	       && exists($tmphash{$lspec[$l]}{$anchors[1]})){
		push @dells, $lspec[$l];
	    }
	    else{
		push @missings, $lspec[$l];
	    }
	}
	else{print STDERR "species should exist: $lspec[$l] \n";}
    }
    if(scalar @missings > 0){
	my $misstr = join(',', @missings);
	if(exists $psemis{$misstr}){$psemis{$misstr} += $mini;}
	else{$psemis{$misstr} = $mini;}
    }
    if(scalar @dells > 0){
	my $delstr = join(',',@dells);
	if(exists $psedels{$delstr}){$psedels{$delstr} += $mini;}
	else{$psedels{$delstr} = $mini;}
    }
   
}

my $rmtmpcmd = "rm $pseudel2check $del2check $tmpfilelist";
readpipe("$rmtmpcmd");


#get output files which are the input to countEvents
my $matchout = "$outpath\/matches\.txt";
my $duplout = "$outpath\/duplications\.txt";
my $insout = "$outpath\/insertions\.txt";
my $pseout = "$outpath\/pseudomatches\.txt";
my $psemisout = "$outpath\/pseudomissing\.txt";
my $psedelout = "$outpath\/pseudodeletions\.txt";
my $pseinsout = "$outpath\/pseudoinsertions\.txt";
my $delout = "$outpath\/deletions\.txt";
my $misout = "$outpath\/missinganchors\.txt";


open(my $outd,">>",$duplout);

foreach my $du (sort keys %dupevents) {
    my @devs = split ',', $du;
    @devs = sort @devs;
    my $dustr = join(',',@devs);
    print $outd "$dustr\t$dupevents{$du}\n";
}
close $outd;

open(my $outm,">>",$matchout);
foreach my $ma (sort keys %matevents) {
    my @mevs = split ',', $ma;
    @mevs = sort @mevs;
    my $mastr = join(',',@mevs);
    print $outm "$mastr\t$matevents{$ma}\n";
}
close $outm;

open(my $outi,">>",$insout);
foreach my $in (sort keys %insevents) {
    my @ievs = split ',', $in;
    @ievs = sort @ievs;
    my $instr = join(',',@ievs);
    print $outi "$instr\t$insevents{$in}\n";
}
close $outi;


open(my $outpi,">>",$pseinsout);
foreach my $pin (sort keys %pseins) {
    my @pievs = split ',', $pin;
    @pievs = sort @pievs;
    my $pinstr = join(',',@pievs);
    print $outpi "$pinstr\t$pseins{$pin}\n";
}
close $outpi;


open(my $outp,">>",$pseout);
foreach my $mi (sort keys %psevents) {
    my @mevs = split ',', $mi;
    @mevs = sort @mevs;
    my $mistr = join(',',@mevs);
    print $outp "$mistr\t$psevents{$mi}\n";
}
close $outp;

open(my $oute,">>",$delout);
foreach my $dl (sort keys %delevents) {
    my @dls = split ',', $dl;
    @dls = sort @dls;
    my $dlstr = join(',',@dls);
    print $oute "$dlstr\t$delevents{$dl}\n";
}
close $oute;

open(my $outn2,">>",$misout);
foreach my $ms (sort keys %misevents) {
    my @miss = split ',', $ms;
    @miss = sort @miss;
    my $msstr = join(',',@miss);
    print $outn2 "$msstr\t$misevents{$ms}\n";
}
close $outn2;

open(my $outpm,">>",$psemisout);
foreach my $pm (sort keys %psemis) {
    my @pmiss = split ',', $pm;
    @pmiss = sort @pmiss;
    my $pmsstr = join(',',@pmiss);
    print $outpm "$pmsstr\t$psemis{$pm}\n";
}
close $outpm;

open(my $outpd,">>",$psedelout);
foreach my $pd (sort keys %psedels) {
    my @pdls = split ',', $pd;
    @pdls = sort @pdls;
    my $pdlstr = join(',',@pdls);
    print $outpd "$pdlstr\t$psedels{$pd}\n";
}
close $outpd;


#count events
##singletoncount = spec-num=spec-num...!spec-pnum=spec-pnum=...
##totelem and nonestr as singleton

##create summarypath and itolout
my $cmditol = "mkdir $outpath\/data_iTOL";
readpipe("$cmditol");
my $cmdsum = "touch $outpath\/geneticEvents\.txt";
readpipe("$cmdsum");
my $cmdtree = "touch $outpath\/OutTree\.txt";
readpipe("$cmdtree");
my $countcmd = "perl $toolpath\/countEvents.pl $newicktree $allsinglestr $matchout $duplout $insout $pseout $psemisout $psedelout $pseinsout $delout $misout $outpath\/OutTree\.txt $outpath\/geneticEvents\.txt $allspecstr $allnonestr $outpath\/data_iTOL 2>> $errors";
readpipe("$countcmd");

my $avnum=0;
if($cluscount>0){$avnum = sprintf "%.2f",$sumelems/$cluscount;}
print "Average number of elements per cluster in joined cluster: $avnum\n";

#cograph and noncograph information
my $cgnodenumsum = 0;
my $cgedgenumsum = 0;
my $cgdensitysum = 0;
my $cggraphnum = 0;
open CGF,"<$cglist" or die "can't open $cglist \n";
while(<CGF>){
    my $cgline = $_;
    if($cgline =~ /^graph/){next;}
    my @CG = split '\t', $cgline;
    $cggraphnum++;
    $cgnodenumsum+=$CG[1];
    $cgedgenumsum+=$CG[2];
    $cgdensitysum+=$CG[4];
}



my $ncgnodenumsum = 0;
my $ncgedgenumsum = 0;
my $ncgdensitysum = 0;
my $ncggraphnum = 0;
open NCGF,"<$ncglist" or die "can't open $ncglist \n";
while(<NCGF>){
    my $ncgline = $_;
    if($ncgline =~ /^graph/){next;}
    my @NCG = split '\t', $ncgline;
    $ncggraphnum++;
    $ncgnodenumsum+=$NCG[1];
    $ncgedgenumsum+=$NCG[2];
    $ncgdensitysum+=$NCG[4];
}

my $avcgnodenum = 0;
my $avcgedgenum = 0;
my $avcgdensity = 0;

if($cggraphnum>0){
    $avcgnodenum = sprintf "%.2f",$cgnodenumsum/$cggraphnum;
    $avcgedgenum = sprintf "%.2f",$cgedgenumsum/$cggraphnum;
    $avcgdensity = sprintf "%.2f",$cgdensitysum/$cggraphnum;
}


my $avncgnodenum = 0;
my $avncgedgenum = 0;
my $avncgdensity = 0;

if($ncggraphnum>0){
    $avncgnodenum = sprintf "%.2f",$ncgnodenumsum/$ncggraphnum;
    $avncgedgenum = sprintf "%.2f",$ncgedgenumsum/$ncggraphnum;
    $avncgdensity = sprintf "%.2f",$ncgdensitysum/$ncggraphnum;
}


print "\n";
print "Information on graph structures and corrections:\n";
print "Number of cographs $cggraphnum, with on average\n";
print "$avcgnodenum nodes, $avcgedgenum edges and a density of $avcgdensity.\n";
print "Number of non-cographs $ncggraphnum, with on average\n";
print "$avncgnodenum nodes, $avncgedgenum edges and a density of $avncgdensity.\n";
print "All non-cographs were corrected to obtain a cograph structure.\n";
print "\n";


sub getRightBlock{
    #get THE right block. if it is several different ones, return -1 as the cluster seems to be too diverse
    my @inp = @_;
    my $cluster = $inp[0];
    my @c = split ';', $cluster;
    my $curmin = -1;
    for(my $i=0;$i<scalar @c;$i++){
	my @F = split '\t', $c[$i];
	my $rightanchor;
	if($F[5] < $F[6]){$rightanchor = $F[6];}
	else{$rightanchor = $F[5];}
	if($curmin==-1){$curmin = $rightanchor;}
	elsif($curmin != $rightanchor){return -1;}
	else{}
    }
    return $curmin;
}

sub getSpecList{
    my @inp = @_;
    my $cluslist = $inp[0];
    my @C = split ';', $cluslist;
    my %spec = ();
    for(my $c=0;$c<scalar @C;$c++){
	my @F = split '\t', $C[$c];
	my @G = split '_', $F[1];
	my $sp = $G[0];
	if(exists($spec{$sp})){}
	else{$spec{$sp}=1;}
    }
    return (keys %spec);
}

sub checkBlocks{
    my @inp = @_;
    ##check for intersecting species
    my %spec1 = ();
    my @block1 = split ';', $inp[0];
    for(my $i=0;$i<scalar @block1;$i++){
	my @F = split '\ţ', $block1[$i];
	my @S = split '_', $F[1];
	my $spec = $S[0];
	$spec1{$spec}=1;
    }
    my @block2 = split ';', $inp[1];
    for(my $j=0;$j<scalar @block2;$j++){
	my @G = split '\t', $block2[$j];
	my @T = split '_', $G[1];
	my $spec2 = $T[0];
	if(exists($spec1{$spec2})){return 1;}
    }
    return 0;
}

sub getMaxDist{
    #input is separated by ; and then by tab
    my @inp = @_;
    my %spec2max = ();
    my %spec2min = ();
    my @F = split ';', $inp[0];
    my $minpos = -1;
    my $maxpos = 0;
    for(my $i=0;$i<scalar @F;$i++){
	my @G = split '\t', $F[$i];
	my $start;
	my $end;
	my @H = split '_', $G[1];
	my $spec = $H[0];
	if($G[2] < $G[3]){
	    $start = $G[2];
	    $end = $G[3];
	}
	else{
	    $start = $G[3];
	    $end = $G[2];
	}
	if(exists($spec2min{$spec})){if($start < $spec2min{$spec}){$spec2min{$spec}=$start;}}
	else{$spec2min{$spec}=$start;}

	if(exists($spec2max{$spec})){if($end > $spec2max{$spec}){$spec2max{$spec}=$end;}}
	else{$spec2max{$spec}=$end;}
    }
    my %spec2dist = ();
    ##check all species that have an entry in both and return the distances
    foreach my $k (keys %spec2max){
	if(exists($spec2min{$k})){
	    $spec2dist{$k} = $spec2max{$k} - $spec2min{$k};
	}
    }
    return %spec2dist;
}
