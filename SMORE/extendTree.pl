#!/usr/bin/perl -w
#call: extendTree.pl treefile ids outtree
#return treefile with extended tree

#this program extends the tree file in case that inner node identifiers are
#missing. We need dummy identifier for the inner nodes in order to be
#able to draw the tree with iTOL.


use Data::Dumper;
use strict;
use warnings;

my $treefile = shift;
my $ids = shift;
my $outfile = shift;

open(my $outf, ">>", $outfile);


#we assume the id file to be a file translating taxnames into MAF ids
#header line: code    | taxid   |  primary taxid   | taxname |       multiZid
#in the tree file, there would be taxnames, we want to translate to multizid

my $dotrans = 0;
my %translation = ();
if($ids ne "="){
    open IF,"<$ids" or die "can't open $ids\n";
    while(<IF>){
	chomp;
	my $line = $_;
	if($line eq ""){next;}
	$line =~ s/[\s]+//g;
	my @L = split '\|', $line;
	$translation{$L[3]}=$L[4];
	$dotrans = 1;
    }
}




my $tree="";

open TF,"<$treefile" or die "can't open $treefile\n";
while(<TF>){
    chomp;
    $tree = $_;
    last;
}

if($tree eq ""){print "tree format doesn't fit!\n"; exit 1;}

my @T = split "", $tree;
my @Tnew = ();
my $cbrac = ")";
my $obrac = "(";
my $com = ",";
my $semco = ";";
my $space = " ";
my $count = 0;
my $newtree = "";
my $tmp = "";
if($T[-1] eq $semco){}
else{print "tree format doesn't fit! No \';\' at the end! \n"; exit 1;}

push @Tnew, $T[0];
for(my $i  = 1; $i < scalar @T; $i++){
    if($T[$i] eq $obrac || $T[$i] eq $com || $T[$i] eq $semco || $T[$i] eq $cbrac){ 
	if($tmp eq "" && $T[$i-1] eq $cbrac){
	    my $dumstr = "inode$count";		
	    push @Tnew, "$dumstr";
	    $count++;
	}
	if($tmp ne ""){
	    ##check here if translations contains the current node
	    my $tmp2="";
	    if($dotrans == 1){
		if(exists($translation{$tmp})){
		    $tmp2 = $translation{$tmp};
		}
		else{$tmp2 = $tmp;}
	    }
	    push @Tnew, $tmp2;
	    $tmp="";
	}
	push @Tnew, $T[$i];
    }
    else{
	$tmp = "$tmp$T[$i]";
    }
}


$newtree = join('',@Tnew);

print $outf $newtree;
