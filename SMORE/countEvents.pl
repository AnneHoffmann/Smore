#!/usr/bin/perl -w

# countEvents will receive counts of events happening in different numbers and combinations of species
# the events will be added at the tree by adding an event at the lca of all species and add deletions in missing species
# deletions will be collected in a list in order to later check if they are real deletions or based on missing data

use Data::Dumper;
use strict;
use warnings;
use List::Util qw(any none first);
use POSIX qw(strftime);


my $treefile = shift;
my $singletoncount = shift;
my $matches = shift;
my $dupl = shift;
my $ins = shift;
my $pseudo = shift; #matches
my $pseumis = shift;
my $pseudels = shift;
my $pseuins = shift;
my $dels = shift;
my $misses = shift;
my $treeout = shift;
my $summary = shift;
my $totelemstr = shift; ##for later, printing iTOL files, includes the total number of pseudogenes
my $nonestr = shift; ##for later, printing iTOL filesincludes the total number of pseudogenes
my $iTOLout = shift; ##for later, printing iTOL files

open(my $outt,">>",$treeout);
open(my $outs,">>",$summary);

my %plusnodes = (); ##nodes of the tree with numbers (+)
my %minusnodes = (); ##nodes of the tree with numbers (-)

my %duplications = ();
my %insertions = ();
my %pseudos = (); #matches

my %singletons = ();
my %missing_data = ();

open TF,"<$treefile" or die "can't open $treefile\n";

my $tree="";

while(<TF>){
    chomp;
    $tree = $_;
    last;
}
    
if($tree eq ""){print STDERR "countEvents: tree format doesn't fit!\n"; exit 1;}

print $outs "tree: $tree \n";

##split the tree into an array of its elements
my @T = (); ##tree with its components
my @N = (); ##at each position there is a number showing opening brackets - closing brackets before this position, except for ( ) , ; then -2
my @L = (); #leaves
my @tr = split '', $tree;
my $brackets = 0;
my $tmp = "";
for(my $i = 0; $i < scalar @tr; $i++)
{
    if($tr[$i] eq ')' || $tr[$i] eq '(' || $tr[$i] eq ',' || $tr[$i] eq ';'){
	if($tmp ne ""){
	    push @T, $tmp; 
	    push @N, $brackets;
	    $plusnodes{$tmp}=0;
	    $minusnodes{$tmp}=0;
	    $duplications{$tmp} = 0;
	    $insertions{$tmp} = 0;
	    $pseudos{$tmp} = 0;
	    if($T[(scalar @T) -2] ne ")"){ #leaves
		push @L, $tmp;
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


#define root ids
my $rootid=-1;
for(my $nn = (scalar @N) -1;$nn>=0;$nn--){
    if($N[$nn]==0){$rootid =$nn;last;}
}

if($rootid == -1){
    print STDERR "countEvents: no root in the tree!\n"; exit;
}

##pseudogenes

open PS,"<$pseudo" or die "can't open $pseudo\n";

while(<PS>){
    chomp;
    my $psline = $_;
    if($psline eq ""){next;}
    my @PS = split '\t', $psline;
    if(scalar @PS < 2){print STDERR "misformatted input line (pseudo), will be skipped! line= $psline \n"; next;}
    my @SS = split ',', $PS[0];
    my $psnum = $PS[1];
    if(scalar @SS == 0){next;}
    elsif(scalar @SS == 1){
	$pseudos{$SS[0]} += $psnum;
    }
    else{
	my $Tstr = join('=', @T);
	my $Sstr = join('=',@SS);
	my @pslcas = findLCA($Tstr,$Sstr);
	for(my $i=0;$i<scalar @pslcas;$i++){
	    if(exists($pseudos{$T[$pslcas[$i]]})){$pseudos{$T[$pslcas[$i]]} += $psnum;}
	    else{$pseudos{$T[$pslcas[$i]]} = $psnum;}
	}
    }
}


my %totcounts = (); #only for debugging if numbers in input files (matches...) fit

##MATCHES
##analyse match nodes first as they include several species


open FA,"<$matches" or die "can't open $matches\n";

while(<FA>){
    chomp;
    my $line = $_;
    if($line eq ""){next;}
    my @F = split '\t', $line; ##should have at least 2! entries
    if(scalar @F < 2){print STDERR "misformatted input line (matches), will be skipped! line= $line \n"; next;}
    
    my @S = split ',', $F[0]; ##species
    my $num = $F[1];
    my @nums = (); ##get bracket numbers for each species
    my @ids = (); ##indices in the tree array
    for(my $s = 0;$s<scalar @S;$s++){
	if(exists($totcounts{$S[$s]})){$totcounts{$S[$s]}+=$num;}
	else{$totcounts{$S[$s]}=$num;}
    }
    
    my $maxbracket = 0;
    
    ##Find lowest common ancestor (lca) for the species    
    ##sort numbers, indices, species
    for(my $j = 0;$j < scalar @S; $j++){
	my $idx;
	for(my $jj = 0; $jj < scalar @T; $jj++){
	    if($T[$jj] eq $S[$j]){
		$idx = $jj;
		last;
	    }
	}
	if($N[$idx] > $maxbracket){$maxbracket = $N[$idx];}
	push @ids, $idx;
	push @nums, $N[$idx];
    }
    
    my $Treestr = join('=', @T);
    my $speciestr = join('=',@S);


#####################only adding nodes#######################
# unse findLCA instead of findParentofAll as findParent includes the deletions, but we want to have them explicitly    
#    my @leafsToAdd = findParentOfAll($Treestr,$speciestr);
    my @leafsToAdd = findLCA($Treestr,$speciestr);
    for(my $i=0;$i<scalar @leafsToAdd; $i++){
	if(exists($plusnodes{$T[$leafsToAdd[$i]]})){$plusnodes{$T[$leafsToAdd[$i]]} += $num;}
	else{$plusnodes{$T[$leafsToAdd[$i]]} = $num;}
    }
}

#deletions
############only deleting leaves############################
open FD,"<$dels" or die "can't open $dels\n";

while(<FD>){
    chomp;
    my $line = $_;
    if($line eq ""){next;}
    my @F = split '\t', $line; ##should have at least 2! entries
    if(scalar @F < 2){print STDERR "misformatted input line (dels), will be skipped! line = $line \n"; next;}
    my @S = split ',', $F[0]; ##species
    my $num = $F[1];

    if(scalar @S == 1){
	if(exists($minusnodes{$S[0]})){$minusnodes{$S[0]} += $num;}
	else{$minusnodes{$S[0]} = $num;}
	next;
    }
    
    my $Treestr = join('=', @T);
    my $speciestr = join('=',@S);
    my @leafsToDel = findParentOfAll($Treestr,$speciestr);
    for(my $i=0;$i<scalar @leafsToDel; $i++){
	if(exists($minusnodes{$leafsToDel[$i]})){$minusnodes{$leafsToDel[$i]} += $num;}
	else{$minusnodes{$leafsToDel[$i]} = $num;}
    }
}


##pseudogenes
my %psdels = ();

open FDP,"<$pseudels" or die "can't open $pseudels\n";

while(<FDP>){
    chomp;
    my $pline = $_;
    if($pline eq ""){next;}
    my @FP = split '\t', $pline; ##should have at least 2! entries
    if(scalar @FP < 2){print STDERR "misformatted input line (pseudel), will be skipped! line = $pline\n"; next;}
    my @PS = split ',', $FP[0]; ##species
    my $pnum = $FP[1];


    if(scalar @PS == 1){
	if(exists($psdels{$PS[0]})){$psdels{$PS[0]} += $pnum;}
	else{$psdels{$PS[0]} = $pnum;}
	next;
    }
    my $Treestr2 = join('=', @T);
    my $speciestr2 = join('=',@PS);


    my @pseuleafsToDel = findParentOfAll($Treestr2,$speciestr2);
    for(my $i=0;$i<scalar @pseuleafsToDel; $i++){
	if(exists($psdels{$pseuleafsToDel[$i]})){$psdels{$pseuleafsToDel[$i]} += $pnum;}
	else{$psdels{$pseuleafsToDel[$i]} = $pnum;}
    }
}

my %pseusingles = ();

#singletoncount
my @allS = split '!', $singletoncount;
my @SC = split "=", $allS[0];
for(my $s = 0; $s < scalar @SC; $s++){
    if($SC[$s] eq ""){next;}
    my @stmp = split "-", $SC[$s];
    $singletons{$stmp[0]} = $stmp[1];
}

if(scalar @allS > 1){
my @SP = split "=", $allS[1];
for(my $sp = 0; $sp < scalar @SP; $sp++){
    if($SP[$sp] eq ""){next;}
    my @sptmp = split "-", $SP[$sp];
    $pseusingles{$sptmp[0]} = $sptmp[1];
}
}

##DUPLICATIONS
open FD,"<$dupl" or die "can't open $dupl\n";

while(<FD>){
    chomp;
    my $dline = $_;
    if($dline eq ""){next;}
    my @D = split '\t', $dline; ##should have at least 2! entries
    if(scalar @D < 2){print STDERR "misformatted input line in duplications, will be skipped! line = $dline \n"; next;}
    if(exists($duplications{$D[0]})){$duplications{$D[0]} += $D[1];}
    else{$duplications{$D[0]} = $D[1];}
    if(exists($totcounts{$D[0]})){$totcounts{$D[0]}+=$D[1];}
    else{$totcounts{$D[0]}=$D[1];}
}

    
##INSERTIONS
open FI,"<$ins" or die "can't open $ins\n";

while(<FI>){
    chomp;
    my $iline = $_;
    if($iline eq ""){next;}
    my @I = split '\t', $iline; ##should have at least 2! entries
    if(scalar @I < 2){print STDERR "misformatted input line in insertions, will be skipped! line = $iline \n"; next;}
    if(exists($insertions{$I[0]})){ $insertions{$I[0]}+= $I[1];}
    else{ $insertions{$I[0]} = $I[1];}
    if(exists($totcounts{$I[0]})){$totcounts{$I[0]}+=$I[1];}
    else{$totcounts{$I[0]}=$I[1];}
}

##pseuins
my %psein = ();

open FPI,"<$pseuins" or die "can't open $pseuins\n";

while(<FPI>){
    chomp;
    my $piline = $_;
    if($piline eq ""){next;}
    my @PI = split '\t', $piline; ##should have at least 2! entries
    if(scalar @PI < 2){print STDERR "misformatted input line in pseudoinsertions, will be skipped! line = $piline \n"; next;}
    if(exists($psein{$PI[0]})){$psein{$PI[0]} += $PI[1];}
    else{$psein{$PI[0]} = $PI[1];}
    if(exists($totcounts{$PI[0]})){$totcounts{$PI[0]}+=$PI[1];}
    else{$totcounts{$PI[0]}=$PI[1];}
}

#missing data
open FM,"<$misses" or die "can't open $misses\n";

while(<FM>){
    chomp;
    my $mline = $_;
    if($mline eq ""){next;}
    my @M = split '\t', $mline;
    if(scalar @M < 2){print STDERR "misformatted input line in missing, will be skipped! line = $mline \n"; next;}
    my @specs = split ',', $M[0];
    my $num = $M[1];
    for(my $s = 0;$s<scalar @specs;$s++){
	if(exists($missing_data{$specs[$s]})){$missing_data{$specs[$s]}+=$num;}
	else{$missing_data{$specs[$s]}=$num;}

    }

}



##pseumis
my %psemis = ();

open FPM,"<$pseumis" or die "can't open $pseumis\n";

while(<FPM>){
    chomp;
    my $pmline = $_;
    if($pmline eq ""){next;}
    my @PM = split '\t', $pmline; ##should have at least 2! entries
    if(scalar @PM < 2){print "misformatted input line in pseudomissing, will be skipped! line = $pmline \n"; next;}
    my @pmspecs = split ',', $PM[0];
    for(my $pm = 0; $pm < scalar @pmspecs;$pm++){
	if(exists($psemis{$pmspecs[$pm]})){$psemis{$pmspecs[$pm]} += $PM[1];}
	else{$psemis{$pmspecs[$pm]} = $PM[1];}
    }
}



##enter all the numbers in the tree
##each node has a tuple [a+i,-b,dx,py]
##with a = matches, i = insertions, b = deletions, 
##x = duplications whereas d stays as letter for indication, 
##y as pseudogenes whereas p stays as letter for indication



my %totnumbers = ();
my %totpseudos = ();

my @allT = split '!', $totelemstr;
my @GT = split '=', $allT[0];
for(my $gt=0;$gt < scalar @GT; $gt++){
    if($GT[$gt] eq ""){next;}
    my @HT = split '-', $GT[$gt];
    $totnumbers{$HT[0]} = $HT[1];
}

my @GTP = split '=', $allT[1];
for(my $gtp=0;$gtp < scalar @GTP; $gtp++){
    if($GTP[$gtp] eq ""){next;}
    my @HTP = split '-', $GTP[$gtp];
    $totpseudos{$HTP[0]} = $HTP[1];
}




my %nonenums = ();
my %pseunons = ();
my @allG = split '!', $nonestr;
my @GG = split '=', $allG[0];
for(my $gg=0;$gg<scalar @GG;$gg++){
    if($GG[$gg] eq ""){next;}
    my @HG = split '-', $GG[$gg];
    $nonenums{$HG[0]} = $HG[1];
}

if(scalar @allG > 1){
my @GP = split '=', $allG[1];
for(my $gp=0;$gp<scalar @GP;$gp++){
    if($GP[$gp] eq ""){next;}
    my @HP = split '-', $GP[$gp];
    $pseunons{$HP[0]} += $HP[1];
}
}






my @newT = ();
for(my $t=0; $t < scalar @N; $t++){
    my $vert = $T[$t];
    if($N[$t] >= 0){ ##not a -2, thus a node in the tree
	##check that the values we wanna get from the hash really exist
	if(exists($plusnodes{$vert})){}else{$plusnodes{$vert}=0;}
	if(exists($insertions{$vert})){}else{$insertions{$vert}=0;}
	my $allplus = $plusnodes{$vert} + $insertions{$vert};
	if(exists($minusnodes{$vert})){}else{$minusnodes{$vert}=0;}
	if(exists($totnumbers{$vert})){}else{$totnumbers{$vert} = 0;}
	if(exists($duplications{$vert})){}else{$duplications{$vert}=0;}
	if(exists($singletons{$vert})){}else{$singletons{$vert}=0;}
	if(exists($pseudos{$vert})){}else{$pseudos{$vert}=0;}
	if(exists($nonenums{$vert})){}else{$nonenums{$vert}=0;}
	if(exists($missing_data{$vert})){}else{$missing_data{$vert}=0;}
	if(exists($pseusingles{$vert})){}else{$pseusingles{$vert}=0;}
	if(exists($psdels{$vert})){}else{$psdels{$vert}=0;}
	if(exists($psein{$vert})){}else{$psein{$vert}=0;}
	my $pseuplus = $pseudos{$vert} + $psein{$vert};
	if(exists($psemis{$vert})){}else{$psemis{$vert}=0;}
	if(exists($totpseudos{$vert})){}else{$totpseudos{$vert}=0;}
	if(exists($pseunons{$vert})){}else{$pseunons{$vert}=0;}
	my $newname = "$T[$t]\[t$totnumbers{$vert}\($totpseudos{$vert}\)|i$allplus\($pseuplus\)\|l$minusnodes{$vert}\($psdels{$vert}\)\|d$duplications{$vert}\|s$singletons{$vert}\($pseusingles{$vert}\)\|n$nonenums{$vert}\($pseunons{$vert}\)|m$missing_data{$vert}\($psemis{$vert}\)\]";
	push @newT, $newname;
    }
    else{
	push @newT, $T[$t];
    }
}

my %numdiffs = ();
my %pnumdiffs = ();

###############Test numbers, thus sum up all the numbers for each species and check the difference to the total number of elems
for(my $ll=0;$ll<scalar @L;$ll++){
    #go through all leaves
    #add up all the numbers from nodes that have a N-val < N[$ll] and are not leaves, each number appears once!!!
    #find ID of leave in T
    my $curid;
    for(my $t0=0;$t0<scalar @T;$t0++){
	if($L[$ll] eq $T[$t0]){
	    $curid = $t0;
	    last;
	}
    }
    my $allelems = 0;
    my $allpseudos = 0;
    my $elemsum = 0;
    my $pseudosum = 0;
    my $curval = $N[$curid];
    #go through Tnew and N and add up all the numbers
    for(my $t1 = $curid;$t1<scalar @newT;$t1++){
	if($N[$t1] != $curval){next;}
	$curval--;
	my $curel = $newT[$t1];
	my @split1 = split '\[', $curel;
	my @split1b = split '\]', $split1[1];
	my @split2 = split '\|', $split1b[0];
	for(my $s2 = 0;$s2 < scalar @split2; $s2++){
	    #letters are: t,i,-l,d,n,s
	    my @CUR = split '\(', $split2[$s2];
	    my $curnum = $CUR[0];
	    my $psnum = 0;
	    if(scalar @CUR == 2){$psnum = substr($CUR[1],0,length($CUR[1])-1);}
	    if($curnum =~ /^t/){
		$allelems += substr($curnum,1);
		$allpseudos += $psnum;
	    }#total
	    if($curnum =~ /^i/){
		$elemsum += substr($curnum,1);
		$pseudosum += $psnum;
	    }#ins
	    if($curnum =~ /^l/){
		$elemsum -= substr($curnum,1);
		$pseudosum -= $psnum;
	    }#loss
	    if($curnum =~ /^d/){$elemsum += substr($curnum,1);}#dupl
	    if($curnum =~ /^s/){
		$elemsum += substr($curnum,1);
		$pseudosum += $psnum;
	    }#single	    
	    if($curnum =~ /^n/){
		$elemsum += substr($curnum,1);
		$pseudosum += $psnum;
	    }#none
	    if($curnum =~ /^m/){
		$elemsum -= substr($curnum,1);
		$pseudosum -= $psnum;
	    }#missing data
	}
    }
    my $diff = $allelems-$elemsum;
    my $pdiff = $allpseudos-$pseudosum;
    $numdiffs{$L[$ll]} = $diff;
    $pnumdiffs{$L[$ll]} = $pdiff;
}

my $tstr = join ("", @newT);

print $outt "Newick tree with numbers specifying event counts at its nodes.
Each node name contains the following numbers in event specification:
[ta(p)|ib(p)|lc(p)|dx|se(p)|nz(p)|mf(p)] where:
a are total number of elements,
b are insertions,
c are deletions, thus losses, 
x are duplications,
e are singletons, 
p are the pseudogenes from this category
z are the elements that were excluded from the analysis.\n
f is the number of elements that were excluded from the analysis because of missing anchors

Tree:\n";
print $outt "$tstr\n";

#print diff numbers first here:
print $outt "functional genes:\n";
foreach my $dif (keys %numdiffs){
    print $outt "$dif $numdiffs{$dif}\n";
}

print $outt "pseudogenes:\n";
foreach my $pdif (keys %pnumdiffs){
    print $outt "$pdif $pnumdiffs{$pdif}\n";
}


close $outt;


#write summary file

print $outs "The following output shows the number of events that occur at the specified node in the tree.\n";
print $outs "Thus, summary for each event is a tab separated table with: tree_node number_of_events\n";
print $outs "Events are Insertion, Deletion, Duplication, Pseudogenization if chosen at the beginning.\n";
print $outs "\n\n";

##sort the entries in each entry for the hashes alphabetically


my %allins = ();
foreach my $in1 (sort keys %insertions) {
    if(exists($allins{$in1})){$allins{$in1} += $insertions{$in1};}
    else{$allins{$in1} = $insertions{$in1};}
}
foreach my $in2 (sort keys %plusnodes) {
    if(exists($allins{$in2})){$allins{$in2} += $plusnodes{$in2};}
    else{$allins{$in2} = $plusnodes{$in2};}
}

my %allpseins = ();
foreach my $in1 (sort keys %psein) {
    if(exists($allpseins{$in1})){$allpseins{$in1} += $psein{$in1};}
    else{$allpseins{$in1} = $psein{$in1};}
}
foreach my $in2 (sort keys %pseudos) {
    if(exists($allpseins{$in2})){$allpseins{$in2} += $pseudos{$in2};}
    else{$allpseins{$in2} = $pseudos{$in2};}
}


print $outs "EVENT: Insertion\n";
foreach my $in (sort keys %allins) {
    my $tmp="\(0\)";
    if(exists($allpseins{$in})){
	$tmp = "\($allpseins{$in}\)";
    }
    print $outs "$in\t$allins{$in}$tmp\n";
}
print $outs "\n\n";


print $outs "EVENT: Deletion\n";
foreach my $de (sort keys %minusnodes) {
    my $tmp="\(0\)";
    if(exists($psdels{$de})){
	$tmp = "\($psdels{$de}\)";
    }
    print $outs "$de\t$minusnodes{$de}$tmp\n";
}
print $outs "\n\n";


print $outs "EVENT: Duplication\n";
foreach my $du (sort keys %duplications) {
    print $outs "$du\t$duplications{$du}\n";
}
print $outs "\n\n";


print $outs "EVENT: Singletons\n";
foreach my $si (sort keys %singletons) {
    my $tmp="\(0\)";
    if(exists($pseusingles{$si})){
	$tmp = "\($pseusingles{$si}\)";
    }

    print $outs "$si\t$singletons{$si}$tmp\n";
}
print $outs "\n\n";


print $outs "EVENT: No Anchors Beginning\n";
foreach my $ni (sort keys %nonenums) {
    my $tmp="\(0\)";
    if(exists($pseunons{$ni})){
	$tmp = "\($pseunons{$ni}\)";
    }

    print $outs "$ni\t$nonenums{$ni}$tmp\n";
}
print $outs "\n\n";

print $outs "EVENT: No Anchors Inbetween\n";
foreach my $md (sort keys %missing_data) {
    my $tmp="\(0\)";
    if(exists($psemis{$md})){
	$tmp = "\($psemis{$md}\)";
    }

    print $outs "$md\t$missing_data{$md}$tmp\n";
}
print $outs "\n\n";


print $outs "CHECK: Differences between total numbers and sum of events\n";
foreach my $oi (sort keys %numdiffs) {
    my $tmp="\(0\)";
    if(exists($pnumdiffs{$oi})){
	$tmp = "\($pnumdiffs{$oi}\)";
    }

    print $outs "$oi\t$numdiffs{$oi}$tmp\n";
}
print $outs "\n\n";


##TODO add extra files for singletons and numdiffs?

####print output files for iTOL input
my $treeout2 = "$iTOLout\/F0tree.txt";
my $namesout = "$iTOLout\/F1tree_labels.txt";
my $totout = "$iTOLout\/F2tree_labels_tot.txt";
my $insout = "$iTOLout\/F3tree_labels_ins.txt";
my $sinout = "$iTOLout\/F4tree_labels_sin.txt";
my $delout = "$iTOLout\/F5tree_labels_del.txt";
my $dupout = "$iTOLout\/F6tree_labels_dup.txt";
#my $pseout=  "$iTOLout\/F7tree_labels_pse.txt";
my $nonout=  "$iTOLout\/F7tree_labels_non.txt";
my $misout=  "$iTOLout\/F8tree_labels_missing.txt";

open(my $outt2,">>",$treeout2);
open(my $outn,">>",$namesout);
open(my $outx,">>",$totout);
open(my $outi,">>",$insout);
open(my $outz,">>",$sinout);
open(my $outd,">>",$delout);
open(my $outu,">>",$dupout);
open(my $outo,">>",$nonout);
open(my $outm,">>",$misout);


print $outt2 $tree;

my $header_names = 
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL, Internal labels
COLOR,#000000\n";
    
my $header_tot = 
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL, Sum of elements & Legend
COLOR,#000000\n";

my $header_ins =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Insertions
COLOR,#00ff00\n";

my $header_sin =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Singletons
COLOR,#ffa500\n";

my $header_del =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Deletions
COLOR,#ff0000\n";

my $header_dup =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Duplications
COLOR,#0000ff\n";

my $header_non =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Excluded
COLOR,#551a8b\n";

my $header_mis =
"DATASET_TEXT
SEPARATOR COMMA
DATASET_LABEL,Missing
COLOR,#A9A9A9\n";


my $inbetweentext =
"#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#
 
#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#
    
#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap. Used only for text labels which are displayed on the outside
MARGIN,1.0
    
#applies to external text labels only; if set, text labels associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
SHOW_INTERNAL,0
   
#Rotate all labels by the specified angle
LABEL_ROTATION,0
    
#applies to external text labels only; If set to 1, labels will be displayed in arcs aligned to the tree (in circular mode) or vertically (in normal mode). All rotation parameters (global or individual) will be ignored.
ALIGN_TO_TREE,0
    
#font size factor; For external text labels, default font size will be slightly less than the available space between leaves, but you can set a multiplication factor here to increase/decrease it (values from 0 to 1 will decrease it, values above 1 will increase it)
SIZE_FACTOR,1
    
#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the \"DATA\" keyword              #
#=================================================================#
#the following fields are possible for each node:
#ID,label,position,color,style,size_factor,rotation
    
#position defines the position of the pie chart on the tree:
#  -1 = external label
#  a number between 0 and 1 = internal label positioned at the specified value along the node branch (for example, position 0 is exactly at the start of node branch, position 0.5 is in the middle, and position 1 is at the end)
#style can be 'normal',''bold','italic' or 'bold-italic'
#size factor will be multiplied with the standard font size\n";

my $legend =
"
LEGEND_TITLE,Element Counts
LEGEND_SHAPES,1,1,1,1,1,1,1
LEGEND_COLORS,#000000,#00ff00,#ffa500,#ff0000,#0000ff,#551a8b,#A9A9A9
LEGEND_LABELS,Total,Insertions,Singletons,Deletions,Duplications,Excluded,Missing_Data
";



print $outn "$header_names$inbetweentext\n";
print $outn "\n\nDATA\n";

print $outx "$header_tot$inbetweentext\n";
print $outx "$legend\n";
print $outx "\n\nDATA\n";


print $outi "$header_ins$inbetweentext";
print $outi "\n\nDATA\n";

print $outz "$header_sin$inbetweentext";
print $outz "\n\nDATA\n";


print $outd "$header_del$inbetweentext";
print $outd "\n\nDATA\n";

print $outu "$header_dup$inbetweentext";
print $outu "\n\nDATA\n";

print $outo "$header_non$inbetweentext";
print $outo "\n\nDATA\n";

print $outm "$header_mis$inbetweentext";
print $outm "\n\nDATA\n";



for(my $tt=0;$tt < scalar @T;$tt++){
    if($N[$tt] >= 0){##we do include the artificial root to debug
	my $n1 = $T[$tt];
	my $isleaf = 0;
	for(my $ll=0;$ll < scalar @L;$ll++){
	    if($L[$ll] eq $T[$tt]){
		$isleaf = 1;
		last;
	    }
	}
	if($isleaf){
	    ##leaf node, different line printing
	    my $nline = "$totnumbers{$T[$tt]}\($totpseudos{$T[$tt]}\)";
	    print $outx "$T[$tt],$nline,-1,#000000,normal,1,0\n";
	    my $sumi = $insertions{$T[$tt]} + $plusnodes{$T[$tt]};
	    my $psumi = $pseudos{$T[$tt]} + $psein{$T[$tt]};
	    print $outi "$T[$tt],$sumi\($psumi\),-1,#00ff00,normal,1,0\n";
	    print $outz "$T[$tt],$singletons{$T[$tt]}\($pseusingles{$T[$tt]}\),-1,#ffa500,normal,1,0\n";
	    print $outd "$T[$tt],$minusnodes{$T[$tt]}\($psdels{$T[$tt]}\),-1,#ff0000,normal,1,0\n";
	    print $outu "$T[$tt],$duplications{$T[$tt]},-1,#0000ff,normal,1,0\n";
	    print $outo "$T[$tt],$nonenums{$T[$tt]}\($pseunons{$T[$tt]}\),-1,#551a8b,normal,1,0\n";
	    print $outm "$T[$tt],$missing_data{$T[$tt]}\($psemis{$T[$tt]}\),-1,#A9A9A9,normal,1,0\n";		
	}
	else{
	    ##inner nodes, different line printing
	    print $outn "$T[$tt],$T[$tt],0,#000000,normal,1,0\n";
	    my $sumi2 = $insertions{$T[$tt]} + $plusnodes{$T[$tt]};
	    my $psumi2 = $pseudos{$T[$tt]} + $psein{$T[$tt]};
	    print $outi "$T[$tt],$sumi2\($psumi2\),0.6,#00ff00,normal,1,0\n";
	    print $outd "$T[$tt],$minusnodes{$T[$tt]}\($psdels{$T[$tt]}\),0.8,#ff0000,normal,1,0\n";
	    #duplications do not occur at inner nodes
	}
    }
}


##Find father of a node n in a newick tree:
##give each node in the newick string a positive number
##the number of opening gaps - closing gaps before the node
##the father of the current nodes is the next node in the newick
##string with a number one smaller than the one of the current node

##find lca: if both have the same number:
##simultaneously find the father and compare them
##if they have different numbers: start with the one with
##the higher number and check if the detected father is the node
##with the lower number. if the numbers become the same and the father is
##not found, go to the first case


sub findLCA{
    my @inp= @_;
    #tree and leaf string is separated by =
    my @T = split '=', $inp[0];
    my @L = split '=', $inp[1]; #names of the species
    my @Ltmp = split '=', $inp[1];
    
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

	if(scalar @pdiff == 0){
	    print STDERR "Error in countEvents, findLCA; pathes: $pstr1; $pstr2\n";
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
	print STDERR "Error in countEvents, findLCA: strange curLCA\n";
    }

    return $curlca;
}




##given a list of leaves and a tree, find the lca of the leaves s.t. there is NO child of the parent that is NOT in the given leaf list
##This function will be used to determine the deletions and the pseudogenes
##the return value is again a list of nodes in the tree
sub findParentOfAll{
    my @inp= @_;
    #tree and leaf string is separated by =
    my @T = split '=', $inp[0];
    my @L = split '=', $inp[1]; #names of the species
    my @Ltmp = split '=', $inp[1];
    
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
    #start with the leaf which is deepest in the tree
    #check if all its siblings are in the working array, too.
    #if yes, push father in current working array.
    #if not, push current node in output array
    while(scalar @Ltmp > 0){
	@L2tmp = ();
	@L2tmpids = ();
	my $maxdepth = 0;
	my $maxdid = -1;
	for(my $l=0;$l<scalar @Ltmp;$l++){
	    if($N[$Ltmpids[$l]] > $maxdepth){
		$maxdepth = $N[$Ltmpids[$l]];
		$maxdid = $l;
	    }
	}
	#get siblings to $L[$maxdid]
	#take all nodes from the tree that have equal N until there is a lower N in between
	my $goodup = 0; #set this to one if sibling are found to be in Ltmp or no siblings
	my $gooddown = 0;
	my $up= $Ltmpids[$maxdid] + 1;
	my $down = $Ltmpids[$maxdid]-1;
	my @goodsibs = (); #ids of sibs in T
	while($up < scalar @T){
	    if($N[$up]<0){$up++;next;}
	    if($N[$up]<$maxdepth)
	    {
		$goodup=1;
		last;
	    }
	    elsif($N[$up]>$maxdepth){$up++;next;}
	    else{
		my $found = 0;
		#sibling, check if it is containe in Ltmp
		for(my $lt = 0;$lt < scalar @Ltmp;$lt++){
		    if($T[$up] eq $Ltmp[$lt]){
			push @goodsibs, $up; 
			$found = 1;
			last;
		    }
		}
		$up++;
		if($found==0){last;}
	    }
	}
	while($down >= 0){
	    if($N[$down]<0 && $down>0){$down--;next;}
	    if($N[$down]<$maxdepth || $down==0)
	    {
		$gooddown=1;
		last;
	    }
	    elsif($N[$down]>$maxdepth){$down--;next;}
	    else{
		my $found = 0;
		#sibling, check if it is containe in Ltmp
		for(my $lt = 0;$lt < scalar @Ltmp;$lt++){
		    if($T[$down] eq $Ltmp[$lt]){
			push @goodsibs, $down; 
			$found = 1;
			last;
		    }
		}
		$down--;
		if($found==0){last;}
	    }
	}
	push @goodsibs, $Ltmpids[$maxdid];
	#update Ltmp to L2tmp without elems in goodsibs
	for(my $li = 0;$li < scalar @Ltmpids;$li++){
	    my $lexists = 0;
	    for(my $gs = 0; $gs < scalar @goodsibs;$gs++){
		if($goodsibs[$gs] == $Ltmpids[$li]){$lexists=1;last;}
	    }
	    if($lexists == 0){
		push @L2tmpids, $Ltmpids[$li];
		push @L2tmp, $Ltmp[$li];
	    }
	}
	
	#if good in both directions, push father in L2tmp, delete goodsibs from Ltmp
	if($goodup == 1 && $gooddown == 1 && scalar @Ltmp > 1){
	    my @path1 = ();
	    my @path2 = ();
	    my $num1 = $N[$goodsibs[0]];
	    my $num2 = $N[$goodsibs[1]];
	    for(my $t4=$goodsibs[0];$t4<scalar @T;$t4++){
		if($N[$t4] != $num1){next;}
		$num1--;
		if($num1 < 0){push @path1, $rootid;last;}
		push @path1, $t4;	
	    }
	    for(my $t5=$goodsibs[1];$t5<scalar @T;$t5++){
		if($N[$t5] != $num2){next;}
		$num2--;
		if($num2 < 0){push @path2, $rootid;last;}
		push @path2, $t5;	
	    }
	    ##look for smallest node that overlaps and push it in L2tmp
	    for(my $f1=0;$f1<scalar @path1;$f1++){
		my $ff = 0;
		for(my $f2=0;$f2<scalar @path2;$f2++){
		    if($path1[$f1] == $path2[$f2]){
			push @L2tmpids, $path2[$f2];
			push @L2tmp, $T[$path2[$f2]];
			$ff=1;
			last;
		    }
		}
		if($ff==1){last;}
	    }
	}
	#else push goodsibs and curnode in output, delete from Ltmp
	else{
	    for(my $gg=0;$gg<scalar @goodsibs;$gg++){
		push @output, $T[$goodsibs[$gg]];
	    }
	}
	
	#update Ltmp, empty L2tmp
	@Ltmp = @L2tmp;
	@Ltmpids = @L2tmpids;
    }

    return @output;
    
} 
