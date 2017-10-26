#!/usr/bin/perl -w

# mix runs between prep and toast
# hence, it takes prep output and produces a gene list based on the join parameter
# additionally, if the number of clusters is very high, mix will partition the clusters into sevearl files with up to a maximum number of clusters given by the user. The output will include a list of commands to call the next step of the pipeline -smoreroast- based on the cluster partition.
# in this way, jobs can run in parallel started independently by the user

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
my $joinmode="relaxed";
my $specieslist;
my $maxclusnum = 50000; #maximal number of cluster for a roast run
#options for creating the roast command(s)
my $seqsim=0.8;
my $strucsim=0.8;
my $newicktree;
my $errors;
my $nocheck;
my $printc;
my $printall;


GetOptions(
    #all modes
    'tool|t=s' => \$toolpath,
    'out|o=s' => \$outpath,
    'python=s' => \$pythonpath,
    'perl=s' => \$perlpath,
##options for toast and bake
    'prep|i=s' => \$pathtocam, #path to output of smore prep (genes folder)
    'join=s' => \$joinmode,
    'species=s' => \$specieslist,
    'max=f' => \$maxclusnum,
    #options for creating roast commands
    'seqsim|s=f' => \$seqsim,
    'strucsim|p=f' => \$strucsim,
    'newick=s' => \$newicktree,
    'err=s' => \$errors,
    'nomiss' => \$nocheck,
    ##options for verbose
    'verbose' => \$printall,
    'clus' => \$printc 
    ) or die "error in smoremix";


my $mode = 1;

my %blocks = ();
my %leftblocks = ();
my %species = ();
my %pseuspecies = ();
my %nones = ();
my %pseudonones = ();
my %nonesNT = ();
my %pseudononesNT = ();

my %allTypes = (); #hash with spec_type -> num
my %pallTypes = (); #hash with spec_type -> num

my %block2spec = ();
my @allspec = ();
my %neighbors = (); #species_anchor -> leftneighbor_rightneighbor or -1 

my %types = (); #this is only the types, not the species

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
#    print STDERR "SORT: $cmdsort \n";
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
	    #	print STDERR "join a,b,c: $curstart, $curend, $rb\n";
	    $curend = $rb;
	    $curblock = "$curblock\;$curb";
	    #	print STDERR "$curblock\n";
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
		#	    print STDERR "curend=-1, lb: $curstart, rb: $rb\n";
		#	    print STDERR Dumper(\%rb2clus);
		#	    print STDERR "\n";
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
		    #		print STDERR "join overlapping cluster\n";
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
    my @K = split '_', $bk;
    my @B = split ';', $blocks{$bk};
    my $clussize = scalar @B;
    #write to allClustersList
    print $outjc ">clus$origcount $clussize $K[0] $K[1]\n";
    $origcount++;
    $sumelemorig += $clussize;
    for(my $cb = 0;$cb < $clussize;$cb++){
	print $outjc "$B[$cb]\n";
    }
}
close $outjc;
my $origavsize=0;
if($orignum>0){$origavsize = sprintf "%.2f",$sumelemorig/$orignum;}
print "Average number of elements of original clusters: $origavsize \n";

my $newclusnum = 0;
if($joined == 1){    
    $newclusnum = scalar (keys %joinedblocks);
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

my $clusfolder = "$outpath\/cluster";
if($printc || $printall){
    if(!(-e $clusfolder)){
	my $mkcmd = "mkdir $clusfolder";
	readpipe("$mkcmd");
    }
}

my $eventlist = "$outpath\/allClusters_joined\.txt";
if($joined==1){
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
	print $outg ">clus$cluscount $clussize $leftanchor $rightanchor\n";
	$cluscount++;
	$sumelems+=$clussize;
	for(my $cb = 0;$cb < $clussize;$cb++){
	    print $outg "$B[$cb]\n";
	    if($printc || $printall){print $outcf "$B[$cb]\n";}
	}
    }
    my $avnum=0;
    if($cluscount>0){$avnum = sprintf "%.2f",$sumelems/$cluscount;}
    print "Average number of elements per cluster: $avnum\n";
}

##check number of clusters and do partitioning
my $dopart = 0;
my $cluslist;
if($joined == 1){
    if($newclusnum > $maxclusnum){$dopart =1;}
    $cluslist = $eventlist;
}
else{
    if($orignum > $maxclusnum){$dopart=1;}
    $cluslist = $jclusfile;
}


my $cmdlist = "$outpath\/commandlist\.txt";
open(my $outcmd, ">>", $cmdlist);

my $optstr = "";
if($nocheck){$optstr = "$optstr --nomiss";}


if($dopart == 0){
    ##we are done, print toast command and exit
    my $rstcmd = "smore roast --tool $toolpath --out $outpath --python $pythonpath --perl $perlpath --in $cluslist -s $seqsim -p $strucsim --newick $newicktree --species $specieslist $optstr\n";
    print $outcmd $rstcmd;
    
}
else{
    ##create a subfolder
    ##partition the data
    #go through blocks and print in extra files
    #create toast commands and write into a list
    my $partfolder = "$outpath\/clusterlists";
    if(-e $partfolder){}
    else{
	my $cmd42= "mkdir $partfolder";
	my @out42 = readpipe("$cmd42");
    }
    my $cluscount = 0;
    my $partcount = 0;
    my $partfile = "$partfolder\/clusList\_part$partcount\.txt";
    my $outg2;
    open($outg2,">>",$partfile);
    foreach my $k (keys %blocks){
	if($cluscount >= $maxclusnum){
	    close $outg2;
	    my $rstcmd = "smore roast --tool $toolpath --out $outpath --python $pythonpath --perl $perlpath --in $partfile -s $seqsim -p $strucsim --newick $newicktree --species $specieslist $optstr\n\n";
	    print $outcmd $rstcmd;
	    $partcount++;
	    $partfile = "$partfolder\/clusList\_part$partcount\.txt";
	    open($outg2,">>",$partfile);
	    $cluscount = 0;
	}
	my @A = split '_', $k;
	my $leftanchor = $A[0];
	my $rightanchor = $A[1];
	my @B = split ';', $blocks{$k};
	my $clussize = scalar @B;
	#write to allClustersList
	print $outg2 ">clus$cluscount $clussize $leftanchor $rightanchor\n";
	$cluscount++;
	for(my $cb = 0;$cb < $clussize;$cb++){
	    print $outg2 "$B[$cb]\n";
	}
    }
    #print last command
    close $outg2;
    my $rstcmd = "smore roast --tool $toolpath --out $outpath --python $pythonpath --perl $perlpath --in $partfile -s $seqsim -p $strucsim --newick $newicktree --species $specieslist $optstr\n\n";
    print $outcmd $rstcmd;
}


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


__END__

##old version of getMaxDist
sub getMaxDist{
    #input is separated by ; and then by tab
    my @inp = @_;
    my @F = split ';', $inp[0];
    my $minpos = -1;
    my $maxpos = 0;
    for(my $i=0;$i<scalar @F;$i++){
	my @G = split '\t', $F[$i];
	my $start;
	my $end;
	if($G[2] < $G[3]){
	    $start = $G[2];
	    $end = $G[3];
	}
	else{
	    $start = $G[3];
	    $end = $G[2];
	}
	if($minpos == -1 || $start < $minpos){$minpos = $start;}
	if($end > $maxpos){$maxpos = $end;}
    }
    my $dist = $maxpos - $minpos;
    return $dist;
}
