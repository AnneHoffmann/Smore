#!/usr/bin/perl -w

## perl createAlignments.pl graphfile outpath pathtonw secsim strsim mode numdifftypes outfolder match dupl ins pseudo pseumis pseudels pseuins dels missing newicktree pathoTemp summary remodlingsout inremoldingsout

use Data::Dumper;
use strict;
use warnings;
use List::Util qw(any none first min max);
use POSIX qw(strftime);

my $file = shift;
my $outpath=shift;
my $pathtonw = shift;
my $seqlim = shift;
my $strlim = shift;
my $mode = shift;
my $numdifftypes = shift;
my $nwtree = shift;
my $leftanchor = shift;
my $rightanchor = shift;
my $del2check = shift;
my $pseudel2check = shift;
my $toprint = shift;
my $alnfile = shift;

if(! $toprint){$toprint = 0;}

my $aout;
if($toprint == 1){
    open($aout,">>",$alnfile);
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


my %remoldings = ();
# The pairs of elements (separated with ':') are defined 
#as orthologs based on the similarity score but have distinct types according to the input";
#this is only reported IFF there are more than just one type!

my %inremoldings = ();
#The pairs of elements are defined as orthologs as they have the same types but the similarity is below the orthology threshold.
#this is only reported IFF there are more than just one type!

my %elemcount = (); #count how many elements of which species occur during the creation of alignments

open(my $outd2c,">>", $del2check);
open(my $outp2c,">>", $pseudel2check);


##Counting: build graphs based on duplication alignments
##get connected components and use them as counting "entity"

##one letter codes needed
##are all the chars ok? ##no - or ~ as they appear in the alignment
my @letters = ('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','0','1','2','3','4','5','6','7','8','9');#,'!','#','@','$','%','^','&','*','(',')','=','+','>','<','?','[',']','|'); 
my $letcount = 0;
my $letnum = scalar @letters;

my %node2letter = ();
my %species = ();
my %start2node = ();
my %pgenes = (); ##pseudogenes
my %spec2pseudo = ();   
my @vertices = ();
my @arcs = ();


open FA,"<$file" or die "can't open $file\n";

while(<FA>){
    chomp;
    
    my $line = $_;
    if($line eq ""){next;}
    my @F = split ' ', $line;
    my $n1 = $F[0]; #node 1
    my $n2 = $F[1]; #node 2
    my $seqsim = $F[2];
    my $strsim = $F[3];
    
    my @G = split '_', $n1;
    my @G2 = split '_', $n2;
    my $isp1 = 0; #check if genes are pseudo or not
    my $isp2 = 0;
    
    if($mode == 0){
	#nodes look like: chr_spec_id_start_end_strand_pseudo
	my $spec = $G[(scalar @G) - 6];
	my $spec2 = $G2[(scalar @G2) - 6];
	$species{$spec} = "";
	$species{$spec2} = "";
	$spec2pseudo{$spec} = "";
	$spec2pseudo{$spec2} = "";
	##if we find a pseudogenes, add it to the pseudogene vector and remove it from the analysis
	my $p1 = $G[(scalar @G) - 1];	    
	if($p1 eq "P"){if(exists($pgenes{$n1})){}else{$pgenes{$n1}=1;} $isp1 = 1;}
	my $startvec = $G[(scalar @G) - 4] + (0.0001 * $G[(scalar @G) - 5]);
	if(exists $start2node{$startvec}){}
	else{$start2node{$startvec} = $n1;}
	my $p2 = $G2[(scalar @G2) - 1];
	if($p2 eq "P"){if(exists($pgenes{$n2})){}else{$pgenes{$n2}=1;}$isp2 = 1;}
	my $startvec2 = $G2[(scalar @G2) - 4] + (0.0001 * $G2[(scalar @G2) - 5]);
	if(exists $start2node{$startvec2}){}
	else{$start2node{$startvec2} = $n2;}
    }
    else{
	#nodes look like: chr_spec_id_start_end_strand_type_pseudo
	#	    chr2A_gorGor3_276_9333523_9333524_-_Undet_TRUE
	#	    chr6_GL000252v2_alt_hg38_1276_10814083_10814084_-_Ala_FALSA
	my $spec = $G[(scalar @G) - 7];
	$species{$spec} = "";
	my $spec2 = $G2[(scalar @G2) - 7];
	$species{$spec2} = "";
	$spec2pseudo{$spec} = "";
	$spec2pseudo{$spec2} = "";
	##if we find a pseudogenes, add it to the pseudogene vector and remove it from the analysis
	my $p1 = $G[(scalar @G) - 1];
	my $p2 = $G2[(scalar @G2) - 1];
	if($p1 eq "T" || $p1 eq "t" || $p1 eq "True" || $p1 eq "true"
	   || $p1 eq "1" || $p1 eq "TRUE"){if(exists($pgenes{$n1})){}else{$pgenes{$n1}=1;}$isp1 = 1;}
	my $startvec = $G[(scalar @G) - 5] + (0.0001 * $G[(scalar @G) - 6]);
	if(exists $start2node{$startvec}){}
	else{$start2node{$startvec} = $n1;}
	if($p2 eq "T" || $p2 eq "t" || $p2 eq "True" || $p2 eq "true"
	   || $p2 eq "1" || $p2 eq "TRUE"){if(exists($pgenes{$n2})){}else{$pgenes{$n2}=1;} $isp2 = 1;}
	my $startvec2 = $G2[(scalar @G2) - 5] + (0.0001 * $G2[(scalar @G2) - 6]);
	if(exists $start2node{$startvec2}){}
	else{$start2node{$startvec2} = $n2;}
	#remoldings
	if($numdifftypes > 1){
	    ##check types:
	    my $t1 = $G[(scalar @G) - 2];
	    my $t2 = $G2[(scalar @G2) - 2];
	    
	    if($seqsim >= $seqlim && $strsim >= $strlim && $t1 ne $t2){ 
		##sequence similarity is high but types seem to be different
		my $remstr = "";
		if($n1 le $n2){$remstr = "$n1\:$n2";}
		else{$remstr = "$n2\:$n1";}
		if(exists($remoldings{$remstr})){$remoldings{$remstr}++;}
		else{$remoldings{$remstr}=1;}
	    }
	    if($t1 eq $t2 && $seqsim < $seqlim){
		##equal types but low sequence similarity
		my $inremstr = "";
		if($n1 le $n2){$inremstr = "$n1\:$n2";}
		else{$inremstr = "$n2\:$n1";}
		if(exists($inremoldings{$inremstr})){$inremoldings{$inremstr}++;}
		else{$inremoldings{$inremstr}=1;}
	    }
	}
    }
    
    if($letcount>=scalar @letters){$letcount = 0;}
    
    if($seqsim >= $seqlim && $strsim >= $strlim){
	##similarity fits, nodes get the same letter
	if (exists $node2letter{$n1}) {
	    if(exists $node2letter{$n2}){
		if($node2letter{$n1} eq $node2letter{$n2}){}
		else{$node2letter{$n2} = $node2letter{$n1};} ##should fit as we iterate over all nodes.
	    }
	    else{
		$node2letter{$n2} = $node2letter{$n1};
	    }
	} 
	else{
	    if(exists $node2letter{$n2}){
		$node2letter{$n1} = $node2letter{$n2};
	    }
	    else{
		$node2letter{$n1} = $letters[$letcount];
		$node2letter{$n2} = $letters[$letcount];
		if($letcount >= $letnum){$letcount=0;}
		else{$letcount++;}
	    }
	}
    }
    else{
	##similarity does not fit, nodes get different letters.
	if (exists $node2letter{$n1}) {}
	else{$node2letter{$n1} = $letters[$letcount];
	     if($letcount >= $letnum){$letcount=0;}
	     else{$letcount++;}
	}
	if (exists $node2letter{$n2}) {}
	else{$node2letter{$n2} = $letters[$letcount];
	     if($letcount >= $letnum){$letcount=0;}
	     else{$letcount++;}
	}
    }
    
}#done reading input file


my $specstr = "";
$specstr = join(',',(keys %species));


my $null = "0";
my $eins = "1";


##do a special case if there is only one species!
if(scalar (keys %species) == 1){
    #start2node contains all current nodes
    #pgenes contains all pseudogenes
    #now add the pseudogenes to pseudinsertions and the remaining genes to insertions and return
    my $numpgenes = scalar (keys %pgenes);
    my $numgenes = scalar (keys %start2node);
    my $numrealgenes = $numgenes - $numpgenes;
    my @SP = keys %species;
    if(exists($insevents{$SP[0]})){$insevents{$SP[0]}+=$numrealgenes;}
    else{$insevents{$SP[0]}=$numrealgenes;}
    if(exists($pseins{$SP[0]})){$pseins{$SP[0]}+=$numpgenes;}
    else{$pseins{$SP[0]}=$numpgenes;}
}
else{
    ##for every species, sort nodes in the correct order and create the sequences to be aligned
    foreach my $k1 (sort { $a <=> $b } keys(%start2node) ) {
	my $curnode = $start2node{$k1};
	my @H = split "_", $curnode;
	my $spec = "";
	my $state = "";
	if($mode == 0){
	    $spec = $H[(scalar @H) - 6];
	    $state = $H[(scalar @H) - 1];
	}
	else{
	    $spec = $H[(scalar @H) - 7];
	    $state = $H[(scalar @H) - 1];
	}
	my $pstr = $spec2pseudo{$spec};
	##write numbers in the same order as chars in sequence, 0 for normal, 1 for pseudo
	if($mode == 0){ ##cm mode
	    if($state eq "N"){
		$spec2pseudo{$spec} = "$pstr$null";
	    }
	    else{
		$spec2pseudo{$spec} = "$pstr$eins";
	    }
	}
	else{
	    if($state eq "T" || $state eq "t" || $state eq "True" || $state eq "true"
	       || $state eq "1" || $state eq "TRUE"){
		$spec2pseudo{$spec} = "$pstr$eins";
	    }
	    else{$spec2pseudo{$spec} = "$pstr$null";}
	}
	my $sstr = $species{$spec};
	$species{$spec} = "$sstr$node2letter{$curnode}";
    }

    ##create graphs out of the duplciation alignments, edges are where the columns fit
    ##nodes are: letter_spec_id!!!
    ##edges are alphabetically sorted e1 < e2 (unique edges), written with e1 e2
    
    my $lseq1sum = 0;
    my $jump = 0;
    ##get alignments for each species against each other
    foreach my $k (keys %species){
	my $lseq1 = $species{$k};
	my $lseq1len = length($lseq1);
	$lseq1sum += $lseq1len;
	if(exists($elemcount{$k})){$elemcount{$k}+=$lseq1len;}
	else{$elemcount{$k}=$lseq1len;}
	##skip if it is the same species as there cannot be an edge
	##problem, if the cluster consists of only nodes of the same species!
	##thus, add those elements directly to the insertions
	if(scalar (keys %species) == 1){
	    my $psseq = $spec2pseudo{$k};
	    my $zcount = $psseq =~ tr/0/0/; #normal genes
	    my $ecount = $psseq =~ tr/1/1/; #pseudogenes
	    if(exists $insevents{$k}){$insevents{$k} += $zcount;}
	    else{$insevents{$k} = $zcount;}
	    if(exists($psevents{$k})){$psevents{$k}+=$ecount;}
	    else{$psevents{$k}=$ecount;}
	    $jump=1;
	    last;
	}
	if($toprint == 1){print $aout "\@$k\t$species{$k}\n";}
	my $match=0; #m
	my $pseu=0; #s
	my $dupl=0; #d
	my $ins=0; #i
	my $del=0; #l
	foreach my $k2 (keys %species){
	    if($k eq $k2) {next;}
	    my $lseq2 = $species{$k2};
	    my $lseq2len = length($lseq2);
            my $newls1;
            my $newls2;
            my $remls1;
            my $remls2;
            my $altnwtwice = 0;
            if($lseq1len >= 900 && $lseq2len < 900){
                $newls1 = substr($lseq1,0,900);
                $newls2 = $lseq2;
                $remls1 = substr($lseq1,900);
                $remls2 = "-" x (length($remls1));
            }
            elsif($lseq2len >= 900 && $lseq1len < 900){
                $newls2 = substr($lseq2,0,900);
                $newls1 = $lseq1;
                $remls2 = substr($lseq2,900);
                $remls1 = "-" x (length($remls2));
            }
            elsif($lseq2len < 900 && $lseq1len < 900){
                $newls2 = $lseq2;
                $newls1 = $lseq1;
                $remls2 = "";
                $remls1 = "";
            }
            else{ ##both longer than 900...
                print STDERR "both altNW seqs are longer than 900: lseq1: $lseq1len, lseq2: $lseq2len \n";
                $altnwtwice = 1;
                $newls2 = substr($lseq2,0,900);
                $newls1 = substr($lseq1,0,900);
                $remls2 = substr($lseq2,900);
                $remls1 = substr($lseq1,900);
            }
            ###ALTNW DOES !NOT! WORK WITH SEQUENCES LONGER THAN 1000!!!
            ##split the long sequences as it is unprobable that the first of seq1 and the
            ##last of seq 2 align
	    my $cmd1 = "$pathtonw/altNW 1 1 \"$newls1\" \"$newls2\"";
	    my @out1 = readpipe("$cmd1");
	    if($toprint == 1){print $aout ">$k2\t";}
	    if($toprint == 1){print $aout "@out1";}
	    chomp(@out1);
	    ##count duplication etc for each pos with ~ or - and write it down
	    ##such that we take the max for event (all letters together)
	    my $m=0;
	    my $l=0;
	    my $s=0;
	    my $i=0;
            my $l1out = "";
            my $l2out = "";
            if($altnwtwice == 0){
                $l1out = "$out1[1]$remls1";
                $l2out = "$out1[2]$remls2";
            }
            else{
                my $cmd2 = "$pathtonw/altNW 1 1 \"$remls1\" \"$remls2\"";
                my @out2 = readpipe("$cmd2");
                $l1out = "$out1[1]$out2[1]";
                $l2out = "$out1[2]$out2[2]";
            }
	    my @ref = split '',$l1out;
	    my @rseq = split '', $lseq1;
	    my @oth = split '',$l2out;
	    my @oseq = split '', $lseq2;
	    if(scalar @ref != scalar @oth){
		print STDERR "alignment does not fit! lseq1, lseq2: $lseq1, $lseq2 \n";
		foreach my $n (keys %node2letter){
		    print STDERR "$n\t$node2letter{$n}\n";
		}
		print STDERR "END ERRORMESSAGE\n";
		next;
	    }
	    my $rcount =0;
	    my $ocount=0;
	    for(my $r=0;$r<scalar @ref;$r++){
		if($ref[$r] eq '-' || $ref[$r] eq '~'){}
		else{$rcount++;}
		if($oth[$r] eq '-' || $oth[$r] eq '~'){}
		else{$ocount++;}
	    }
	    #compare sequence lengths with seq len in alignment
	    if($lseq1len != $rcount){print STDERR "lengths don't fit! sequence: $lseq1, alnseq: $out1[1]\n";}
	    if($lseq2len != $ocount){print STDERR "lengths don't fit! sequence: $lseq2, alnseq: $out1[2]\n";}
	    if($lseq1len != $rcount || $lseq2len != $ocount){next;}
	    
	    my @sp12pseudo = split '', $spec2pseudo{$k};  ##this is the sequence for species k for pseudogene or not
	    my @sp22pseudo = split '', $spec2pseudo{$k2};	    
	    my $tild = "~";
	    my $mins = "-";
	    my $rz = 0;
	    my $oz = 0;
	    #should have the same length as it is an alignment
	    #build a graph between the chars of the aligned sequences
	    for(my $z=0;$z < scalar @ref;$z++){
		#create nodes and edges
		my $curc = $ref[$z];
		my $v1;
		my $v2;
		my $ed;
		#do not declare pe1 and p2 here, use the array index directly when needed?
		my $pe1 = $sp12pseudo[$rz];
		my $pe2 = $sp22pseudo[$oz];
		
		if($curc eq $oth[$z])
		{
		    ##add $z behind the nodes names to make them unique
		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		    if($v1 lt $v2){$ed = "$v1 $v2";}
		    else{$ed = "$v2 $v1";}
		    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
		    $rz++;
		    $oz++;

		}
		elsif($curc eq $tild){
		    ##add an edge from the last letter in the current to the current letter in the other
		    my $pe1p = $sp12pseudo[$rz-1];
		    my $smrz = $rz-1;
		    $v1 = "$rseq[$smrz]\_$k\_$pe1p\_$smrz"; ##v1 should already exist!
		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1; print STDERR "v1 should exist! $v1 \n";}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		    if($v1 lt $v2){$ed = "$v1 $v2";}
		    else{$ed = "$v2 $v1";}
		    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
		    $oz++;
		}
		elsif($curc eq $mins){
		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		    $oz++;
		}
		elsif($oth[$z] eq $mins){
		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    $rz++;
		}
		elsif($oth[$z] eq $tild){
		    my $pe2p = $sp22pseudo[$oz-1];
		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    my $smoz = $oz-1;
		    $v2 = "$oseq[$smoz]\_$k2\_$pe2p\_$smoz"; ##v2 should already exist!
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2; print STDERR "v2 should exist! $v2 \n";}
		    if($v1 lt $v2){$ed = "$v1 $v2"; }
		    else{$ed = "$v2 $v1";}
		    if(none {$_ eq $ed} @arcs){push @arcs, $ed;}
		    $rz++;
		}
		else{
		    $v1 = "$curc\_$k\_$pe1\_$rz";
		    $v2 = "$oth[$z]\_$k2\_$pe2\_$oz";
		    if(none {$_ eq $v1} @vertices){push @vertices, $v1;}
		    if(none {$_ eq $v2} @vertices){push @vertices, $v2;}
		    $rz++;
		    $oz++;
		}
	    }
	}

    }
    ##if jump==1, no graph could be built
    if($jump==1){
	print STDERR "No graph could be build..\n";
    }
    ##graph for this cluster is built.
    ##get connected components
    my $vertices = join(',',@vertices);
    my $curvertnum = scalar @vertices;
    if($lseq1sum != $curvertnum){print STDERR "Sequence sum does not fit vertex num: seq: $lseq1sum, vertex: $curvertnum, vertices: $vertices\n";}
    my $arcs = join(',',@arcs);
    my @CCs = connectedComponents($vertices,$arcs);
    
    ##count events for each CC
    my $perccnum = scalar @CCs; #add all nodes of each cc for this graph to see if it fits the curvertnum
    my $ccnodecount = 0;
    my $cccount = 0;
    for(my $ii = 0; $ii < scalar @CCs; $ii++){
	my %spe2count = ();
	my %pse2count = ();
	my @verts = split ' ', $CCs[$ii];
	$ccnodecount += scalar @verts;
	for(my $jj=0;$jj< scalar @verts; $jj++){
	    my @SP = split '_', $verts[$jj];
	    my $spe = $SP[1];
	    my $pseudi = $SP[2];
	    if($pseudi eq "1"){
		if(exists $pse2count{$spe}){
		    $pse2count{$spe}++;
		}
		else{
		    $pse2count{$spe}=1;
		}
	    }
	    else{
		if(exists $spe2count{$spe}){
		    $spe2count{$spe}++;
		}
		else{
		    $spe2count{$spe}=1;
		}
	    }
	    
	}
	##counting: only matches and dup (and pseudogenes?), losses later when adding to the tree
	##add elements from the singleton clusters that were sorted
	##and specify the number of elements per species for the None cluster
	my $spstr = join(',',sort (keys %spe2count));
	my $spnum = scalar (keys %spe2count);	
	my @vals = values %spe2count;
	my $vnum = scalar @vals;
	my $mini = min @vals;
	if(scalar @vals == 0){}
	elsif(scalar @vals == 1){##singleton/insertion
	    if(exists $insevents{$spstr}){$insevents{$spstr} += $vals[0];}
	    else{$insevents{$spstr} = $vals[0];}
	}
	elsif(none {$_ != $vals[0]} @vals)#check if all values are equal
	{
	    if(exists $matevents{$spstr}){$matevents{$spstr} += $vals[0];}
	    else{$matevents{$spstr} = $vals[0];}
	}
	else{##count duplications
	    if(exists $matevents{$spstr}){$matevents{$spstr} += $mini;}
	    else{$matevents{$spstr} = $mini;}
	    foreach my $kk (keys %spe2count){
		my $diff = $spe2count{$kk} - $mini; 
		if($diff > 0){
		    if(exists $dupevents{$kk}){$dupevents{$kk} += $diff;}
		    else{$dupevents{$kk} = $diff;}
		}
	    }
	}
	my @pseuvals = values %pse2count;
	my $psnum = scalar @pseuvals;
	my $psestr = "";
	my $pmin = min @pseuvals;
	if(scalar @pseuvals > 0){
	    $psestr = join(',',sort (keys %pse2count));
	    ##do the same distinguishing for the pseudogenes as for matching
	    if(scalar @pseuvals == 1){#singleton
		if(exists $pseins{$psestr}){$pseins{$psestr} += $pseuvals[0];}
		else{$pseins{$psestr} = $pseuvals[0];}
	    }
	    elsif(none {$_ != $pseuvals[0]} @pseuvals)#all the same
	    {
		if(exists $psevents{$psestr}){$psevents{$psestr} += $pseuvals[0];}
		else{$psevents{$psestr} = $pseuvals[0];}
	    }
	    else{#add the minimum common value for the complete psestr and the remainings as singletons, as the duplications have been counted already
		$pmin = min @pseuvals;
		if(exists $psevents{$psestr}){$psevents{$psestr} += $pmin;}
		else{$psevents{$psestr} = $pmin;}
		foreach my $ppp (keys %pse2count){
		    my $pdiff = $pse2count{$ppp} - $pmin;
		    if($pdiff > 0){
			if(exists $psevents{$ppp}){$psevents{$ppp} += $pdiff; }
			else{$psevents{$ppp} = $pdiff;}	    
		    }
		}
	    

	    }
	}
	#do check for missing data or deletions here
	
	##if cluster >= 2 nodes && num_species >= 2
	##check if all species below the LCA appear in this cluster
	##if a species does not appear, check if
	##the missing species have the anchoring blocks
	##if yes: write the combination of species as a deletion(like matches)
	##if no: ignore this element in the species (add to missing data list)
	
	#do this for every CC as we are looking at genetic events of homologs
	#steps:
	#collect set of species
	#missing species = find LCA and missing species below (sub)
	#check if anchors exist
	#thus here: counting of deletion events and missing data
	if($spnum > 1){
	    #species: comma-separated in $spstr
	    my @missingspecs = getMissingSpecs($spstr,$nwtree);
	    if(scalar @missingspecs > 0){
		my $missstr = join(',',@missingspecs);
		my $d2cstr = "$leftanchor\_$rightanchor\t$missstr\t$mini\n";
		print $outd2c $d2cstr;

	    }
	}

	#do the same for pseudogenes in order to detect missing data and deletions for pseudogenes
	#extra files for pseudo: -singletons, ins, del, ex, mis


	if($psnum > 0){
	    	    #species: comma-separated in $psestr
	    my @pmissingspecs = getMissingSpecs($psestr,$nwtree);
	    if(scalar @pmissingspecs > 0){
		my $pmissstr = join(',',@pmissingspecs);
		my $p2cstr = "$leftanchor\_$rightanchor\t$pmissstr\t$pmin\n";
		print $outp2c $p2cstr;
	    }
	}
	

	
    }#done check for every CC
    #check if ccnodecount == curvertnum
    if($ccnodecount != $curvertnum){
	print STDERR "numbers don't fit! previous nodenum: $curvertnum, now: $ccnodecount\n";
    }
}#end else if more than one species

##sort the entries in each entry for the hashes alphabetically

my $outstring = "";
#instead of writing to files, we put all information in a string and return it
#there is not too much information as its only the counts for one graph

foreach my $du (sort keys %dupevents) {
    my @devs = split ',', $du;
    @devs = sort @devs;
    my $dustr = join(',',@devs);
    $outstring = "$outstring\=$dustr\-$dupevents{$du}";
}
$outstring = "$outstring\n";

foreach my $ma (sort keys %matevents) {
    my @mevs = split ',', $ma;
    @mevs = sort @mevs;
    my $mastr = join(',',@mevs);
    $outstring = "$outstring\=$mastr\-$matevents{$ma}";
}
$outstring = "$outstring\n";


foreach my $in (sort keys %insevents) {
    my @ievs = split ',', $in;
    @ievs = sort @ievs;
    my $instr = join(',',@ievs);
    $outstring = "$outstring\=$instr\-$insevents{$in}";
}
$outstring = "$outstring\n";


foreach my $pin (sort keys %pseins) {
    my @pievs = split ',', $pin;
    @pievs = sort @pievs;
    my $pinstr = join(',',@pievs);
    $outstring = "$outstring\=$pinstr\-$pseins{$pin}";
}
$outstring = "$outstring\n";


foreach my $mi (sort keys %psevents) {
    my @mevs = split ',', $mi;
    @mevs = sort @mevs;
    my $mistr = join(',',@mevs);
    $outstring = "$outstring\=$mistr\-$psevents{$mi}";
}
$outstring = "$outstring\n";

foreach my $dl (sort keys %delevents) {
    my @dls = split ',', $dl;
    @dls = sort @dls;
    my $dlstr = join(',',@dls);
    $outstring = "$outstring\=$dlstr\-$delevents{$dl}";
}
$outstring = "$outstring\n";

foreach my $ms (sort keys %misevents) {
    my @miss = split ',', $ms;
    @miss = sort @miss;
    my $msstr = join(',',@miss);
    $outstring = "$outstring\=$msstr\-$misevents{$ms}";
}
$outstring = "$outstring\n";


foreach my $pm (sort keys %psemis) {
    my @pmiss = split ',', $pm;
    @pmiss = sort @pmiss;
    my $pmsstr = join(',',@pmiss);
    $outstring = "$outstring\=$pmsstr\-$psemis{$pm}";
}
$outstring = "$outstring\n";

foreach my $pd (sort keys %psedels) {
    my @pdls = split ',', $pd;
    @pdls = sort @pdls;
    my $pdlstr = join(',',@pdls);
    $outstring = "$outstring\=$pdlstr\-$psedels{$pd}";
}
$outstring = "$outstring\n";




my @krems = keys %remoldings;
if(scalar @krems > 0){
#    print $outsr "The following pairs of elements (separated with ':') are defined 
#as orthologs based on the similarity score but have distinct types according to the input:\n";
    for(my $i=0;$i< scalar @krems; $i++){
	$outstring = "$outstring$krems[$i]\=";
    }
}
$outstring="$outstring\n";

my @kinrems = keys %inremoldings;

if(scalar @kinrems > 0){
#    print $outsi "The following pairs of elements (separated with ':') are defined 
#as orthologs as they have the same types but they sequence similarity is below the given threshold:\n";
    for(my $i=0;$i< scalar @kinrems; $i++){
	$outstring = "$outstring$kinrems[$i]\=";
    }
}

$outstring="$outstring\n";

print $outstring;

sub connectedComponents{

    my @inp = @_;
    my @uniqedges = split ',', $inp[1];
    my @permnodes1 = split ',', $inp[0];
    my @permnodes2=();

    for(my $i=0;$i<scalar @uniqedges;$i++){
	my @E = split ' ', $uniqedges[$i];
	my $n1 = $E[0];
	my $n2 = $E[1];
	my $id1 =-1;
	my $id2 =-1;
	for(my $j=0;$j<scalar @permnodes1;$j++)
	{
	    if(index($permnodes1[$j],$n1)!= -1) #remember current position of node
	    {
		$id1 = $j;
	    }
	    if(index($permnodes1[$j],$n2)!= -1) #remember current position of node
	    {
		$id2 = $j;
	    }
	}
	for(my $k=0;$k<scalar @permnodes1;$k++)
	{
	    if($k != $id1 && $k != $id2)#node is not in this edge
	    {
		push @permnodes2, $permnodes1[$k];
	    }
	}
	if($id1 != $id2){#node occur in this edge and are in one cc
	    my $str = "$permnodes1[$id1] $permnodes1[$id2]";
	    push @permnodes2, $str;
	}
	else{#node do not belong to the same cc
	    push @permnodes2, $permnodes1[$id1];
	}

	@permnodes1 = ();
	@permnodes1 = @permnodes2;
	@permnodes2 = ();
    }
    
    return @permnodes1;
}



#find LCA
#get diffset of leafs_under_LCA - spstr
sub getMissingSpecs{
    my @inp = @_;
    my $spstr = $inp[0];
    my $treefile = $inp[1];

    open TF,"<$treefile" or die "can't open $treefile\n";
    
    my $tree="";
    
    while(<TF>){
	chomp;
	$tree = $_;
	last;
    }


    
    my @output = ();
    my @species = split ',', $spstr;
    if(scalar @species <= 1){return @output;}


    if($tree eq ""){print STDERR "createAlignments: tree format doesn't fit!\n"; exit 1;}
    ##split the tree into an array of its elements
    my @T = (); ##tree with its components
    my @N = (); ##at each position there is a number showing opening brackets - closing brackets before this position, except for ( ) , ; then -2
    my @L = (); #leaves
    my @Lids = ();
    my @tr = split '', $tree;
    my $brackets = 0;
    my $tmp = "";
    for(my $i = 0; $i < scalar @tr; $i++)
    {
	if($tr[$i] eq ')' || $tr[$i] eq '(' || $tr[$i] eq ',' || $tr[$i] eq ';'){
	    if($tmp ne ""){
		push @T, $tmp; 
		push @N, $brackets;
		#	    print "leaves? $T[(scalar @T) -2] \n";
		if($T[(scalar @T) -2] ne ")"){ #leaves
		    push @L, $tmp;
		    push @Lids, (scalar @T) -1;
		}
		$tmp="";
	    }
	    push @T, $tr[$i];
	    push @N, -2;
	    if($tr[$i] eq '('){$brackets++;}
	    if($tr[$i] eq ')'){$brackets--;}
	}
	else{
	    $tmp = "$tmp$tr[$i]";
	}
    }

    my @leafbutnospec = (); #insert ids of leaves in T: Lids
    for(my $l=0;$l< scalar @L;$l++){
	my $isnotspec = 1;
	for(my $s = 0;$s < scalar @species;$s++){
	    if($species[$s] eq $L[$l]){
		$isnotspec = 0;
		last;
	    }
	}
	if($isnotspec == 1){
	    push @leafbutnospec, $Lids[$l];
	}
    }
    
    my $treestr = join('=',@T);
    my $leafbuitnotspecstr = join('=',@leafbutnospec);
    my $lca = findLCA($treestr,$spstr);
    #lca is the id in T
    #check for all elements in leafbutnospec if lca is on their path to the root
    #if yes, put in output
    for(my $z=0;$z<scalar @leafbutnospec;$z++){
	my $num = $N[$leafbutnospec[$z]];
	for(my $t=$leafbutnospec[$z];$t<scalar @T;$t++){
	    if($N[$t] != $num){next;}
	    if($T[$lca] eq $T[$t]){
		push @output, $T[$leafbutnospec[$z]];
		last;
	    }
	    $num--;
	    if($num < 0){last;}
	}
    }

    my $outputstr = join('=',@output);
    return @output;
}





sub findLCA{
    my @inp= @_;
    #tree and leaf string is separated by =
    my @T = split '=', $inp[0];
    my @L = split ',', $inp[1]; #names of the species
    my @Ltmp = split ',', $inp[1];
    
    my @output = ();
    my $rootid = -1;
    
    my @N = ();
    my @allleaves = ();
    my $brackets = 0;
    my $maxbracket = 0;
    for(my $i = 0; $i < scalar @T; $i++)
    {
	if($T[$i] eq ')' || $T[$i] eq '(' || $T[$i] eq ',' || $T[$i] eq ';'){
	    
	    push @N, -2;
	    if($T[$i] eq '('){$brackets++;}
	    if($T[$i] eq ')'){$brackets--;}
	    if($brackets > $maxbracket){$maxbracket = $brackets;}
	}
	else{
	    push @N, $brackets;
	    if($brackets == 0){$rootid = $i;}
	    if($i>0){
		if($T[$i-1] eq ')'){
		    ##this is an inner node
		}
		else{
		    push @allleaves, $T[$i];
		}
	    }	
	}
    }

    if($rootid == -1){return @output;}
    my @Lids = (); #ids of the leaves in T and N
    for(my $l=0;$l<scalar @L;$l++){
	for(my $t=0;$t<scalar @T;$t++){
	    if($L[$l] eq $T[$t]){
		push @Lids, $t;
	    }
	}	
    }

    my @Ltmpids = @Lids;
    my @L2tmp = ();
    my @L2tmpids = ();
    
    #######################################
    my $pstr1="";
    my $pstr2="";

    ##try to find curlca by always looking at two species at the same time
    ##then, get set of nodes from the pathes from the current leaves to the root and intersect them
    ##if there are several nodes in the intersection, take the one with the highest N
    my $curlca = $Lids[0]; #start with species at $S[0]
    for(my $k=1; $k < scalar @L; $k++){
	my @path1 = ();
	my @path2 = ();
	push @path1, $curlca;
	push @path2, $Lids[$k];
	my $num1 = $N[$curlca];
	for(my $t4=$curlca;$t4<scalar @T;$t4++){
	    #got up in the tree and check the bracketnums
	    if($N[$t4] != $num1){next;}
	    $num1--; #bracketnum
	    if($num1 < 0){push @path1, $rootid;last;}
	    push @path1, $t4;
	}

	my $num2 = $N[$Lids[$k]];
	for(my $t5=$Lids[$k];$t5<scalar @T;$t5++){
	    if($N[$t5] != $num2){next;}
	    $num2--;
	    if($num2 < 0){push @path2, $rootid;last;}
	    push @path2, $t5;
	}

	if(scalar @path1 == 0){}
	else{$pstr1 = join("=",@path1);}
	if(scalar @path2 == 0){}
	else{$pstr2 = join("=",@path2);}
	if(scalar @path1 == 0 || scalar @path2 == 0){
	    next;
	}

	#do the intersection of both
	##as we go from leaves to root, we take the first element that we have in common
	my @pdiff = ();
	for(my $a1 = 0; $a1 < scalar @path1;$a1++){
	    for(my $b1 = 0; $b1 < scalar @path2; $b1++){
		if($path1[$a1] == $path2[$b1]){
		    push @pdiff, $path2[$b1];
#		    last; #this can be deleted later on in order to get all overlapping nodes
		}
	    }

	}

	my $pdiffstr = join('=',@pdiff);
	if(scalar @pdiff == 0){
	    print STDERR "error in createAlignments, findLCA; pathes: $pstr1; $pstr2\n";
	}
	elsif(scalar @pdiff == 1){
	    $curlca = $pdiff[0];
	}
	else{
	    #choose the curlca with the highest N
	    my $curplca = $N[$pdiff[0]];
	    my $curpid = $pdiff[0];
	    for(my $pd=1;$pd<scalar @pdiff;$pd++){
		if($N[$pdiff[$pd]] >= $curplca){
		    $curplca = $N[$pdiff[$pd]];
		    $curpid = $pdiff[$pd];
		}
	    }
	    $curlca = $curpid;
	}
    }

    if($N[$curlca] < 0){
	print STDERR "createAlignments, findLCA: strange current LCA\n";
    }

    return $curlca;
}
