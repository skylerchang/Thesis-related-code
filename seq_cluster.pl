#!/usr/bin/perl

use strict; use warnings;
use Data::Dumper;

open COOR,"<$ARGV[0]";
open SEQ,"<$ARGV[1]";

my @coors=<COOR>;
my @seqs=<SEQ>;
my %results;
my %coors_index = map { $_ => $coors[$_] } reverse 0 .. @coors-1;
my %seqs_index = map { $_ => $seqs[$_] } reverse 0 .. @seqs-1;

foreach my $idx (keys %coors_index){
    push@{$results{$coors_index{$idx}}}, $seqs_index{$idx};
}

foreach my $coor (sort keys %results){

    print ">$coor@{$results{$coor}}";
}


close COOR;
close SEQ;






