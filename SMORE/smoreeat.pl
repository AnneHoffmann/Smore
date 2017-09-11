#!/usr/bin/perl -w

#smoreeat.pl will get the folder with all output files from smoreroast,
#summarize the counts and add them to the tree
#also get folder with output from smoremix for the counts of nones and totals
use IO::Uncompress::Gunzip qw($GunzipError);
use Data::Dumper;
use Getopt::Long qw(GetOptions);
use strict;
use warnings;


##options for all modes
my $toolpath;
my $outfolder;
my $pythonpath;
my $perlpath;
##options for toast and bake
my $pathtocam;
my $mixout;
my $roastout;
my $newicktree;
my $nocheck;

GetOptions(
    #all modes
    'tool|t=s' => \$toolpath,
    'out|o=s' => \$outfolder,
    'python=s' => \$pythonpath,
    'perl=s' => \$perlpath,
##options for toast and bake
    'prep=s' => \$pathtocam,
    'mix=s' => \$mixout,
    'roast=s' => \$roastout,
    'newick=s' => \$newicktree,
    'nomiss' => \$nocheck
    ) or die "error in smoreeat";

my %dupevents =();
my %matevents =();
my %insevents =();
my %misevents =(); #missing data=no anchors
my %delevents = ();#deletions
my %psevents = (); #insertions of pseudogenes
my %psedels = ();
my %psemis = ();
my %pseins = ();
my %singles = (); 
my %pseusingles = (); 
my %total = ();
my %pstotal = ();
my %nones= ();
my %psnones = ();

my $matchout = "$outfolder\/matches\.txt";
my $duplout = "$outfolder\/duplications\.txt";
my $insout = "$outfolder\/insertions\.txt";
my $pseout = "$outfolder\/pseudomatches\.txt";
my $psemisout = "$outfolder\/pseudomissing\.txt";
my $psedelout = "$outfolder\/pseudodeletions\.txt";
my $pseinsout = "$outfolder\/pseudoinsertions\.txt";
my $delout = "$outfolder\/deletions\.txt";
my $misout = "$outfolder\/missinganchors\.txt";



##deletions and missings
my $tmpfilelist = "$outfolder\/tmpfiles";
my $tmplscmd = "ls $pathtocam\/temp \> $tmpfilelist";
my @outlscmd = readpipe("$tmplscmd");

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

#do check for missing anchors
my $delchecklist = "$outfolder\/list\_delCheck";
my $delcheckcmd = "ls $roastout\/delCheck\* > $delchecklist";
readpipe("$delcheckcmd");
open FA,"<$delchecklist" or die "can't open $delchecklist\n";

while(<FA>){
    chomp;
    my $del2check = $_;
    ##now, we check if anchors exist with
    open DC,"<$roastout\/$del2check" or die "can't open $roastout\/$del2check \n";
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
}    

#do check for missing anchors
my $pdelchecklist = "$outfolder\/list\_pseudelCheck";
my $pdelcheckcmd = "ls $roastout\/pseudelCheck\* > $pdelchecklist";
readpipe("$pdelcheckcmd");
open FP,"<$pdelchecklist" or die "can't open $pdelchecklist\n";

while(<FP>){
    chomp;
    my $pseudel2check = $_;
    open PC,"<$roastout\/$pseudel2check" or die "can't open $roastout\/$pseudel2check \n";
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
}
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

##duplications
my $dupllist = "$outfolder\/list\_dupls";
my $duplcmd = "ls $roastout\/duplications\* > $dupllist";
readpipe("$duplcmd");
open DU,"<$dupllist" or die "can't open $dupllist\n";
while(<DU>){
    chomp;
    my $curdup = $_;
    open CDU,"<$roastout\/$curdup" or die "can't open $roastout\/$curdup\n";
    while(<CDU>){
	chomp;
	my @DE = split '\t', $_;
	if(exists($dupevents{$DE[0]})){$dupevents{$DE[0]}+=$DE[1];}
	else{$dupevents{$DE[0]}=$DE[1];}
    }
}
open(my $outd,">>",$duplout);

foreach my $du (sort keys %dupevents) {
    my @devs = split ',', $du;
    @devs = sort @devs;
    my $dustr = join(',',@devs);
    print $outd "$dustr\t$dupevents{$du}\n";
}
close $outd;

##matches
my $matlist = "$outfolder\/list\_mats";
my $matcmd = "ls $roastout\/matches\* > $matlist";
readpipe("$matcmd");
open MA,"<$matlist" or die "can't open $matlist\n";
while(<MA>){
    chomp;
    my $curmat = $_;
    open CMA,"<$roastout\/$curmat" or die "can't open $roastout\/$curmat\n";
    while(<CMA>){
	chomp;
	my @ME = split '\t', $_;
	if(exists($matevents{$ME[0]})){$matevents{$ME[0]}+=$ME[1];}
	else{$matevents{$ME[0]}=$ME[1];}
    }
}

open(my $outm,">>",$matchout);
foreach my $ma (sort keys %matevents) {
    my @mevs = split ',', $ma;
    @mevs = sort @mevs;
    my $mastr = join(',',@mevs);
    print $outm "$mastr\t$matevents{$ma}\n";
}
close $outm;

##pseumatches
my $pmatlist = "$outfolder\/list\_pmats";
my $pmatcmd = "ls $roastout\/pseudomatches\* > $pmatlist";
readpipe("$pmatcmd");
open PMA,"<$pmatlist" or die "can't open $pmatlist\n";
while(<PMA>){
    chomp;
    my $pcurmat = $_;
    open PCMA,"<$roastout\/$pcurmat" or die "can't open $roastout\/$pcurmat\n";
    while(<PCMA>){
	chomp;
	my @PME = split '\t', $_;
	if(exists($psevents{$PME[0]})){$psevents{$PME[0]}+=$PME[1];}
	else{$psevents{$PME[0]}=$PME[1];}
    }
}

open(my $outp,">>",$pseout);
foreach my $mi (sort keys %psevents) {
    my @mevs = split ',', $mi;
    @mevs = sort @mevs;
    my $mistr = join(',',@mevs);
    print $outp "$mistr\t$psevents{$mi}\n";
}
close $outp;

##insertions
my $inslist = "$outfolder\/list\_ins";
my $inscmd = "ls $roastout\/insertions\* > $inslist";
readpipe("$inscmd");
open IA,"<$inslist" or die "can't open $inslist\n";
while(<IA>){
    chomp;
    my $curins = $_;
    open CIA,"<$roastout\/$curins" or die "can't open $roastout\/$curins\n";
    while(<CIA>){
	chomp;
	my @IE = split '\t', $_;
	if(exists($insevents{$IE[0]})){$insevents{$IE[0]}+=$IE[1];}
	else{$insevents{$IE[0]}=$IE[1];}
    }
}
open(my $outi,">>",$insout);
foreach my $in (sort keys %insevents) {
    my @ievs = split ',', $in;
    @ievs = sort @ievs;
    my $instr = join(',',@ievs);
    print $outi "$instr\t$insevents{$in}\n";
}
close $outi;


##pseudoinsertions
my $pinslist = "$outfolder\/list\_pins";
my $pinscmd = "ls $roastout\/pseudoinsertions\* > $pinslist";
readpipe("$pinscmd");
open PIA,"<$pinslist" or die "can't open $pinslist\n";
while(<PIA>){
    chomp;
    my $curpins = $_;
    open CPIA,"<$roastout\/$curpins" or die "can't open $roastout\/$curpins\n";
    while(<CPIA>){
	chomp;
	my @PIE = split '\t', $_;
	if(exists($pseins{$PIE[0]})){$pseins{$PIE[0]}+=$PIE[1];}
	else{$pseins{$PIE[0]}=$PIE[1];}
    }
}
open(my $outpi,">>",$pseinsout);
foreach my $pin (sort keys %pseins) {
    my @pievs = split ',', $pin;
    @pievs = sort @pievs;
    my $pinstr = join(',',@pievs);
    print $outpi "$pinstr\t$pseins{$pin}\n";
}
close $outpi;


#parse none, singles, total
##singles
my $sinlist = "$outfolder\/list\_singles";
my $sincmd = "ls $roastout\/singletons\* > $sinlist";
readpipe("$sincmd");
open IS,"<$sinlist" or die "can't open $sinlist\n";
while(<IS>){
    chomp;
    my $cursin = $_;
    open CIS,"<$roastout\/$cursin" or die "can't open $roastout\/$cursin\n";
    while(<CIS>){
	chomp;
	my @SIE = split '\t', $_;
	my @Sspec = split '_', $SIE[0];
	if(exists($singles{$Sspec[0]})){$singles{$Sspec[0]}+=$SIE[1];}
	else{$singles{$Sspec[0]}=$SIE[1];}
    }
}
my $allsinglestr = "";
foreach my $sin (sort keys %singles) {
    $allsinglestr = "$allsinglestr\=$sin\-$singles{$sin}";
}
$allsinglestr = "$allsinglestr\!";
##pseusingles
my $psinlist = "$outfolder\/list\_psingles";
my $psincmd = "ls $roastout\/pseusingletons\* > $psinlist";
readpipe("$psincmd");
open IP,"<$psinlist" or die "can't open $psinlist\n";
while(<IP>){
    chomp;
    my $curpsin = $_;
    open CIP,"<$roastout\/$curpsin" or die "can't open $roastout\/$curpsin\n";
    while(<CIP>){
	chomp;
	my @SIP = split '\t', $_;
	my @Sspecp = split '_', $SIP[0];
	if(exists($pseusingles{$Sspecp[0]})){$pseusingles{$Sspecp[0]}+=$SIP[1];}
	else{$pseusingles{$Sspecp[0]}=$SIP[1];}
    }
}
foreach my $psin (sort keys %pseusingles) {
    $allsinglestr = "$allsinglestr\=$psin\-$pseusingles{$psin}";
}

my $allnonestr = "";
#nones
my $nonefile = "$mixout\/nones\.txt";
open NF,"<$nonefile" or die "can't open $nonefile\n";
while(<NF>){
    chomp;
    my @NL = split '\t', $_;
    my @Nspec = split '_', $NL[0];
    $allnonestr = "$allnonestr\=$Nspec[0]\-$NL[1]";
}
$allnonestr = "$allnonestr\!";
#pseudonones
my $pnonefile = "$mixout\/pseunones\.txt";
open NFP,"<$pnonefile" or die "can't open $pnonefile\n";
while(<NFP>){
    chomp;
    my @NLP = split '\t', $_;
    my @Nspecp = split '_', $NLP[0];
    $allnonestr = "$allnonestr\=$Nspecp[0]\-$NLP[1]";
}
#total
my $allspecstr = "";
my $totfile = "$mixout\/allTypes\.txt";
open TF,"<$totfile" or die "can't open $totfile\n";
while(<TF>){
    chomp;
    my @TL = split '\t', $_;
    my @Tspec = split '_', $TL[0];
    $allspecstr = "$allspecstr\=$Tspec[0]\-$TL[1]";
}
$allspecstr = "$allspecstr\!";
#pseudototals
my $ptotfile = "$mixout\/allPseudoTypes\.txt";
open TFP,"<$ptotfile" or die "can't open $ptotfile\n";
while(<TFP>){
    chomp;
    my @TLP = split '\t', $_;
    my @Tspecp = split '_', $TLP[0];
    $allspecstr = "$allspecstr\=$Tspecp[0]\-$TLP[1]";
}

#count events
##singletoncount = spec-num=spec-num...!spec-pnum=spec-pnum=...
##totelem and nonestr as singleton
my $err = "$outfolder\/errors";
my $cmderr = "touch $err";
readpipe("$cmderr");

##create summarypath and itolout
my $cmditol = "mkdir $outfolder\/data_iTOL";
readpipe("$cmditol");
my $cmdsum = "touch $outfolder\/geneticEvents\.txt";
readpipe("$cmdsum");
my $cmdtree = "touch $outfolder\/OutTree\.txt";
readpipe("$cmdtree");
my $countcmd = "perl $toolpath\/countEvents.pl $newicktree $allsinglestr $matchout $duplout $insout $pseout $psemisout $psedelout $pseinsout $delout $misout $outfolder\/OutTree\.txt $outfolder\/geneticEvents\.txt $allspecstr $allnonestr $outfolder\/data_iTOL 2>> $err";

readpipe("$countcmd");
