#!/usr/bin/perl -w
#call: checkTree.pl treefile specieslist
#return 1 or 0
#this program checks if the species identifier in specieslist also appear in
#the tree file. if this is not the case, the program stops and returns an error message

use Data::Dumper;
use strict;
use warnings;

my $treefile = shift;
my $specieslist = shift;

my $tree="";

open TF,"<$treefile" or die "can't open $treefile\n";
while(<TF>){
    chomp;
    $tree = $_;
    last;
}

if($tree eq ""){print "tree format doesn't fit!\n"; exit 1;}


open PS,"<$specieslist" or die "can't open $specieslist\n";

my $missingspecs = "";
my $missing = 0;

while(<PS>){
    chomp;
    my $line = $_;
    my @F = split '\/', $line;
    my @G = split '\.', $F[-1];
    my $spec = $G[0];
    if(index($tree, $spec) == -1){
	$missingspecs = "$missingspecs$spec,";
	$missing = 1;
    }
}


if($missing == 1){
    print "$missingspecs\n";
}
else{
    print "none\n";
}
