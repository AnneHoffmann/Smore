#!/usr/bin/perl -w

# smorebake does the same as smoretoast
# but additionally prints all intermediary output files,
# e.g. cluster, graphs, duplication alignments


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
my $skipg;
my $skipa;
my $skipc;
my $nocheck;
my $errors;

GetOptions(
    #all modes
    'tool|t=s' => \$toolpath,
    'out|o=s' => \$outpath,
    'python=s' => \$pythonpath,
    'perl=s' => \$perlpath,
##options for toast and bake
    'prep|i=s' => \$pathtocam, #path to output of smore prep (genes folder)
    'seqsim|s=s' => \$seqsim,
    'strucsim|p=s' => \$strucsim,
    'newick=s' => \$newicktree,
    'join=s' => \$joinmode,
    'nograph' => \$skipg,
    'noaln' => \$skipa,
    'noclus' => \$skipc,
    'nomiss' => \$nocheck,
    'species=s' => \$specieslist,
    'err=s' => \$errors
    ) or die "error in smoretoast";



my $optstr = "";
if($skipa){$optstr = "$optstr --noaln";}
if($skipc){$optstr = "$optstr --noclus";}
if($skipg){$optstr = "$optstr --nograph";}
if($nocheck){$optstr = "$optstr --nomiss";}

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

my $eventlist = "$outpath\/allClusters_joined\.txt";
my $listcmd = "touch $eventlist";
readpipe("$listcmd");
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

my %block2spec = ();
my @allspec = ();
my %neighbors = (); #species_anchor -> leftneighbor_rightneighbor or -1 

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
    $pseuspecies{$spec} = 0;
    #sort file and remove blanks
    my $curfile_nb = "$outpath\/$spec\_nb\.bed";
    my $curfile_sorted = "$outpath\/$spec\_sorted\.bed";
    my $cmdnoblanks = "awk \'NF\' $curfile \> $curfile_nb";
    #sort by chr, elemstart, elemend
    my $col1 = 1;
    my $col3 = 3;
    my $col4 = 4;
    my $coln = "n";
    my $cmdsort = "sort -k$col1,$col1 -k$col3,$col3$coln -k$col4,$col4$coln $curfile_nb \> $curfile_sorted";
#    print STDERR "SORT: $cmdsort \n";
    readpipe("$cmdnoblanks");
    readpipe("$cmdsort");

    #this describes the neighbors
    my $chr2 = "-1";
    my $bl2 = "-1";
    my $chr1 = "-1";
    my $bl1 = "-1";
    #read file to blocks/nones
    open CF, "<$curfile_sorted" or die "can't open $curfile_sorted\n";
    while(<CF>){
	chomp;
	my $line = $_;
	if($line=~/^#/){next;}
	if($line eq ""){next;}
	my @F = split '\t', $line;

	if($mode==1)
	{
	    my $tkey = "$spec\_$F[-3]";
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


my $orignum = scalar (keys %blocks);
print "ORIG CLUSNUM: $orignum \n";


my $numspec = scalar @allspec;
print "Number of species: $numspec \n";

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

my $origavsize = $sumelemorig/$orignum;
print "Average size of original clusters: $origavsize \n";

    if($joined == 1){    
	my $newclusnum = scalar (keys %joinedblocks);
	print "JOINED CLUSNUM: $newclusnum \n";
	%blocks = %joinedblocks;
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

#print "totelemnum: $allspecstr \n";

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
my $clusfolder = "$outpath\/cluster";
if(! $skipc && !(-e $clusfolder)){
    my $mkcmd = "mkdir $clusfolder";
    readpipe("$mkcmd");
}
my $graphfolder = "$outpath\/graph";
if(! $skipg && !(-e $graphfolder)){
    my $mkcmd2 = "mkdir $graphfolder";
    readpipe("$mkcmd2");
}
my $alnfolder = "$outpath\/duplication_alignments";
if(! $skipa && !(-e $alnfolder)){
    my $mkcmd3 = "mkdir $alnfolder";
    readpipe("$mkcmd3");
}


##sort keys of blocks and join adjacent ones?
#my $clusnum = scalar (keys %blocks);
#print "NUM CLUSTER: $clusnum \n";
open(my $outg, ">>",$eventlist);
my $cluscount = 0;
my $sumelems = 0;
foreach my $k (keys %blocks){
    my @A = split '_', $k;
    my $leftanchor = $A[0];
    my $rightanchor = $A[1];
    my @B = split ';', $blocks{$k};
    my $clussize = scalar @B;
    #write to allClustersList
    print $outg ">clus$cluscount $clussize\n";
    my $outcf;
    if(! $skipc){
	my $clusfile = "$clusfolder\/cluster$leftanchor\_$rightanchor\.clus";
	open($outcf, ">>", $clusfile);
    }
    $cluscount++;
    $sumelems+=$clussize;
    for(my $cb = 0;$cb < $clussize;$cb++){
	print $outg "$B[$cb]\n";
	if(! $skipc){print $outcf "$B[$cb]\n";}
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
    if(! $skipg){$tmpfile1 = "$graphfolder\/graph$k\.gr";}
    else{$tmpfile1 = "$outpath\/tmpdata1";}
    
    #build graph, tmpfile1 will contain the output
    my $graphcmd = "perl $toolpath\/buildEdgeList_fast.pl $tmpfile0 $mode $toolpath $strucsim $seqsim $tmpfile1";
#    print "GRAPHCMD: $graphcmd \n";
    readpipe("$graphcmd");
    

    #check graph for cograph or not and edit slightly
    my $cglist = "$outpath\/list\_cographs\.txt";
    my $ncglist = "$outpath\/list\_noncographs\.txt";
    my $cgcmd = "touch $cglist";
    my $ncgcmd = "touch $ncglist";
    readpipe("$cgcmd");
    readpipe("$ncgcmd");

    my $checkcmd = "perl $toolpath\/checkGraph_fast.pl $tmpfile1  $k $seqsim $strucsim $mode $cglist $ncglist";

    ##argument list is getting too long...create a temporary file

    my @newoutgraph = readpipe("$checkcmd");
    my $newgraphstr = $newoutgraph[0];
    
    #sort the non edges graphs
    my $toprint = 0;
    if(! $skipa){$toprint = 1;}
    my $filetoprint = "$alnfolder\/alignment$k\.aln";
    my $alncmd = "perl $toolpath\/createAlignments_fast2.pl $tmpfile1 $outpath $toolpath $seqsim $strucsim $mode 0 $newicktree $leftanchor $rightanchor $del2check $pseudel2check $toprint $filetoprint";
    #before countEvents, check for deletions with the files created here
    #create alignment
    #my $alncmd = "perl $scriptpath\/createAlignments_fast.pl $tmpfile1 $outpath $pathtonw $seqsim $strucsim $mode 0 $nwtree $inpath\/temp $leftanchor $rightanchor";

#    print STDERR "ALNCMD: $alncmd \n";
    my @outaln = readpipe("$alncmd"); #this array contains: dup mat insertion pseins pseudomatch deletion missinganchor missinanchorpseu deletionpseu
    my $rmcmd;
    if(! $skipg){$rmcmd = "rm $tmpfile0";}
    else{$rmcmd = "rm $tmpfile1 $tmpfile0";}
    readpipe("$rmcmd");
    my $alen = scalar @outaln;
#    print "OUTALN: $alen\n";
    my $alnput = join('_',@outaln);
#    print "ALNOUT: $alnput \n";
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
	print $outsr "$currems[$sr]\n";
    }
    my @curinrems = split '=', $outaln[10];
    for(my $ri =0;$ri<scalar @curinrems;$ri++){
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
	#   open CTF,"<$curtmpfile" or die "can't open $curtmpfile\n";
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

#if(exists($tmphash{$species})){
#if(exists($tmphash{$species}{$anchor})){
##add to deletions
#}
#else{
##add to missings
#}

#		if(scalar @missingtmp > 0){
#		    my $misstr = join(',', @missingtmp);
#		    if(exists $misevents{$misstr}){$misevents{$misstr} += $mini;}
#		    else{$misevents{$misstr} = $mini;}
#		}
#		if(scalar @missingdel > 0){
#		    my $delstr = join(',',@missingdel);
#		    if(exists $delevents{$delstr}){$delevents{$delstr} += $mini;}
#		    else{$delevents{$delstr} = $mini;}
#		}





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

open(my $outn,">>",$misout);
foreach my $ms (sort keys %misevents) {
    my @miss = split ',', $ms;
    @miss = sort @miss;
    my $msstr = join(',',@miss);
    print $outn "$msstr\t$misevents{$ms}\n";
}
close $outn;

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

my $err = "$outpath\/errors";
my $cmderr = "touch $err";
readpipe("$cmderr");

##create summarypath and itolout
my $cmditol = "mkdir $outpath\/data_iTOL";
readpipe("$cmditol");
my $cmdsum = "touch $outpath\/geneticEvents\.txt";
readpipe("$cmdsum");
my $cmdtree = "touch $outpath\/OutTree\.txt";
readpipe("$cmdtree");
my $countcmd = "perl $toolpath\/countEvents.pl $newicktree $allsinglestr $matchout $duplout $insout $pseout $psemisout $psedelout $pseinsout $delout $misout $outpath\/OutTree\.txt $outpath\/geneticEvents\.txt $allspecstr $allnonestr $outpath\/data_iTOL 2>> $err";
#print "COUNTING: $countcmd\n";
readpipe("$countcmd");



print "Program run ended.\n";
print "Number of clusters: $cluscount\n";
    my $avnum = $sumelems/$cluscount;
print "Average number of elements per cluster: $avnum\n";

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
