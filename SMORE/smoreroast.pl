#!/usr/bin/perl -w

# smoreroast is similar to toast but with clusterlist as input and bare
# event count files as output
# used to run smoremix and after roast summarize all bar output files
# deletions are collected in a separate file as they will be checked for missing data in the next step (smoreeat)

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
my $cluslistpath;
my $seqsim=0.8;
my $strucsim=0.8;
my $newicktree;
my $specieslist;
my $nocheck;
##options for verbose
my $printall;
my $printg;
my $printa;


GetOptions(
    #all modes
    'tool|t=s' => \$toolpath,
    'out|o=s' => \$outpath,
    'python=s' => \$pythonpath,
    'perl=s' => \$perlpath,
##options for toast and bake
    'in|i=s' => \$cluslistpath,
    'seqsim|s=f' => \$seqsim,
    'strucsim|p=f' => \$strucsim,
    'newick=s' => \$newicktree,
    'species=s' => \$specieslist,
    'nomiss' => \$nocheck,
##options for verbose
    'verbose' => \$printall,
    'graph' => \$printg,
    'aln' => \$printa
    ) or die "error in smoreroast";

my @C = split '/', $cluslistpath;
my $cluslist = $C[-1];

my $mode = 1;

my $err = "$outpath\/errors\_$cluslist";
my $errcmd = "touch $err";
readpipe("$errcmd");

my $del2check = "$outpath\/delCheck\_$cluslist";
my $pseudel2check = "$outpath\/pseudelCheck\_$cluslist";

my $remlist = "$outpath\/remoldings\_$cluslist";
my $inremlist = "$outpath\/inremoldings\_$cluslist";
my $remcmd = "touch $remlist";
my $inremcmd = "touch $inremlist";
readpipe("$remcmd");
readpipe("$inremcmd");

open(my $outsr,">>", $remlist);
open(my $outsi,">>", $inremlist);


my $cglist = "$outpath\/list\_cographs\_$cluslist";
my $ncglist = "$outpath\/list\_noncographs\_$cluslist";
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

my %blocks=();

#open cluslist and directly work with each cluster.
open FA,"<$cluslistpath" or die "can't open $cluslistpath\n";

my $leftanchor;
my $rightanchor;
my $cursize;
while(<FA>){
    chomp;
    my $line = $_;
    if($line=~/^>/){
	my @R=split ' ', $line;
	$leftanchor = $R[2];
	$rightanchor = $R[3];
	$cursize = $R[1];
	next;
    }
    if($cursize == 1){
	my @F = split '\t', $line;
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
    else{
	my $key = "$leftanchor\_$rightanchor";
	if(exists($blocks{$key})){
	    my $pf = $blocks{$key};
	    $blocks{$key} = "$pf\;$line";
	}
	else{$blocks{$key} = "$line";}
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

foreach my $k (keys %blocks){
    my @A = split '_', $k;
    my $leftanchor = $A[0];
    my $rightanchor = $A[1];
    #have a tmpfile for blocks string
    my $tmpfile0 = "$outpath\/tmp0$k";
    open(my $outtmp0, ">>", $tmpfile0);
    print $outtmp0 $blocks{$k};
    close $outtmp0;
    my $tmpfile1;
    if($printg || $printall){$tmpfile1 = "$graphfolder\/graph$k\.gr";}
    else{$tmpfile1 = "$outpath\/tmp1$k";}


    #build graph, tmpfile1 will contain the output
    my $graphcmd = "perl $toolpath\/buildEdgeList_fast.pl $tmpfile0 $mode $toolpath $strucsim $seqsim $tmpfile1 2>>$err";

    my @graphout = readpipe("$graphcmd");
    my $nodenum = $graphout[0];

    #check graph for cograph or not and edit slightly

    my $checkcmd = "perl $toolpath\/checkGraph_fast.pl $tmpfile1  $k $seqsim $strucsim $mode $cglist $ncglist 2>>$err";

    my @newoutgraph = readpipe("$checkcmd");
    my $newgraphstr = $newoutgraph[0];


    my $toprint = 0;
    if($printa || $printall){$toprint = 1;}
    my $filetoprint = "$alnfolder\/alignment$k\.aln";

    
    my $alncmd = "perl $toolpath\/createAlignments_fast2.pl $tmpfile1 $outpath $toolpath $seqsim $strucsim $mode 0 $newicktree $leftanchor $rightanchor $del2check $pseudel2check $toprint $filetoprint 2>>$err";
    #before countEvents, check for deletions with the files created here
    #create alignment
    #my $alncmd = "perl $scriptpath\/createAlignments_fast.pl $tmpfile1 $outpath $pathtonw $seqsim $strucsim $mode 0 $nwtree $inpath\/temp $leftanchor $rightanchor";

    my @outaln = readpipe("$alncmd"); #this array contains: dup mat insertion pseins pseudomatch deletion missinganchor missinanchorpseu deletionpseu
    my $rmcmd;
    if($printg || $printall){$rmcmd = "rm $tmpfile0";}
    else{$rmcmd = "rm $tmpfile1 $tmpfile0";}
    readpipe("$rmcmd");


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
}#end FA
close $outsi;
close $outsr;



##print singletons to file and create allsinglestr
my $singleout = "$outpath\/singletons\_$cluslist";
my $pseusingleout = "$outpath\/pseusingletons\_$cluslist";
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


#get output files which are the input to countEvents
my $matchout = "$outpath\/matches\_$cluslist";
my $duplout = "$outpath\/duplications\_$cluslist";
my $insout = "$outpath\/insertions\_$cluslist";
my $pseout = "$outpath\/pseudomatches\_$cluslist";
my $pseinsout = "$outpath\/pseudoinsertions\_$cluslist";

#dels and missing are empty as the check about missing anchors has to be done first

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

