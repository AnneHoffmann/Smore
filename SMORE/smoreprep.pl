#!/usr/bin/perl -w

use Data::Dumper;
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use POSIX qw(strftime);
use File::Basename;
use File::Find;


##options for all modes
my $toolpath;
my $outpath;
my $pythonpath;
my $perlpath;
##options for prep
my $refspecies;
my $genomes;
my $mafs;
my $cmfile;
my $genelist;
my $locilist;
my $perc;
my $evalin;
my $bitvalin;
my $infernalpath;
my $pseudoscore;
my $errors;

GetOptions(
#all modes
    'tool|t=s' => \$toolpath,
    'out|o=s' => \$outpath,
    'python=s' => \$pythonpath,
    'perl=s' => \$perlpath,
##options for prep
    'ref|r=s' => \$refspecies,
    'genomes|g=s' => \$genomes,
    'maf|m=s' => \$mafs,
    'cm|c=s' => \$cmfile,
    'genes=s' => \$genelist,
    'loci=s' => \$locilist,
    'filter=f' => \$perc,
    'incE=f' => \$evalin,
    'incT=f' => \$bitvalin,
    'infernal=s' => \$infernalpath,
    'pseudo=f' => \$pseudoscore,
    'err=s' => \$errors
    ) or die "Some parameter for smore prep doesn't fit! \n";

##in theory, all parameters were checked in startsmore
##another check is probably not wrong

##Outpath
my $outpathstr = "";
if ($outpath){$outpathstr = "-o $outpath ";}
else{print STDERR "No output path given (option -o)!\n"; exit 1;}
##if output folder doesnt exist, create it!
if(-e $outpath){}
else{
    my $cmd42= "mkdir $outpath";
    my @out42 = readpipe("$cmd42");
}

##Refspecies
my $refspeciesstr = "";
if ($refspecies){$refspeciesstr="--ref $refspecies ";
		 #check if file exists
		 if($genomes){
		     my $n1 = "$refspecies\.fa";
		     my $n2 = "$refspecies\.fasta";
		     my $n3 = "$refspecies\.fa\.gz";
		     my $n4 = "$refspecies\.fasta\.gz";
		     if(-l "$genomes\/$n1"){}
		     elsif(-f "$genomes\/$n2"){}
		     elsif(-f "$genomes\/$n3"){}
		     elsif(-f "$genomes\/$n4"){}
		     elsif(-f "$genomes\/$refspecies"){}
		     else{print STDERR "No genome file for reference species $refspecies in $genomes (option -f)\n"; exit 1;}
		 }
}
else{print STDERR "No references species given! (option --ref)\n"; exit 1;}

my $pystr = "";
my $pestr = "";
##Python
if ($pythonpath){$pystr = "--python $pythonpath ";}
else{print STDERR "No path to python3 given! (option --python)\n"; exit 1;}

##Perl
if ($perlpath){$pestr = "--perl $perlpath ";}
else{print STDERR "No path to perl given! (option --perl)\n"; exit 1;}

my $mode = -1;

###Parameters for infernal mode

##CM
my $cmoption = "";
my $cmoptstr="";
if($cmfile){$cmoption = "-sg $cmfile"; $cmoptstr = "--cm $cmfile "; $mode = 0;}


##Genomes
my $genomesstr="";
if ($genomes){ 
    $genomesstr = "--genomes $genomes ";
    $cmoption = "$cmoption $genomes";
    #check if folder is not empty
    if(-e $genomes){} else{print STDERR "Genomes folder is empty! (option -g)\n"; exit 1;}
}

#infernal
my $infstr = "";
if ($infernalpath){
    $infstr = "--infernal $infernalpath ";
    $cmoption = "$cmoption $infernalpath";
}

##Maf
my $mafstr = "";
if ($mafs){$mafstr = "--maf $mafs ";}



##infernal parameter
my $inclopt="";
my $incloptstr = "";
if ($evalin && $cmfile){$inclopt = "-incE $evalin"; $incloptstr = "$incloptstr\--incE $evalin ";}
if ($bitvalin && $cmfile){$inclopt = "-incT $bitvalin";$incloptstr = "$incloptstr\--incT $bitvalin ";}


##Percentage of the lowest scoring blocks based on the MAF scores to be removed (between 0 and 100)
my $filterstr = "";
if ($perc){$filterstr = "--filter $perc ";}


if($cmfile){
    if($genomes && $mafs && $infernalpath){}
    else{print STDERR "Infernal mode but some required parameters are missing! See --help for more information! \n"; exit 1;}
}


###Genelist mode

if($genelist){$cmoption = "-og $genelist"; $cmoptstr = "--genes $genelist "; $mode = 1;}

##maf and filter should already be reported

if($genelist){
    if($mafs){}
    else{print STDERR "Genelist mode but --maf parameter is missing! See --help for more information! \n"; exit 1;}
}


my $prepcmd = "";
if(! $perc){
    $prepcmd = "$pythonpath\/python3 $toolpath\/main.py $cmoption $inclopt $mafs $outpath $toolpath $refspecies $pythonpath 2>> $errors";
}
else{
    $prepcmd = "$pythonpath\/python3 $toolpath\/main.py $perc $cmoption $inclopt $mafs $outpath $toolpath $refspecies $pythonpath 2>> $errors";
}

my @out = readpipe("$prepcmd");
print STDERR (join("",@out));

