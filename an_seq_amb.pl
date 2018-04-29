#!/usr/bin/perl

use strict; use warnings;
use List::Util qw(reduce);
use Carp;
use Data::Dumper;


open ANCHOR,"<$ARGV[0]";
open SEQ,"<$ARGV[1]";
my @anchors=<ANCHOR>;
my @seqs=<SEQ>;
my %results;
my @ratio1;
foreach my $seq (@seqs){
    my @HDs;
    foreach my $anchor (@anchors){
        my $HD = ( $seq ^ $anchor ) =~ tr/\0//c;
        push @HDs, $HD;
    }
    my $min_idx = get_min_idx(\@HDs);
    my $min_val = $HDs[$min_idx];
    $HDs[$min_idx] = 100;
    my $min_idx2 = get_min_idx(\@HDs);
    my $min_val2 = $HDs[$min_idx2];
    my $stability_idx = sprintf("%.3f",$min_val / $min_val2);
    if ($stability_idx == 1) {
        push @ratio1,$seq;
    }else{
        push @{$results{$anchors[$min_idx]}}, "$seq";
    }
}

foreach my $anchor ( sort keys %results){
    print ">$anchor@{$results{$anchor}}";
}

print ">ambiguous\n@ratio1";
close ANCHOR;
close SEQ;


#find index of lowest value in array (courtesy of Stack Overflow)
sub get_min_idx {
    my $ra = shift;
    croak "Sub expects array reference" if ref $ra ne 'ARRAY';
    return reduce { $ra->[$a] < $ra->[$b] ? $a : $b } 0..$#$ra;
}




