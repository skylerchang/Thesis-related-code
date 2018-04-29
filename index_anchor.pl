#!/usr/bin/perl

use strict; use warnings;
use List::Util qw(reduce);
use Carp;

open ANCHOR,"<$ARGV[0]";
open SEQ,"<$ARGV[1]";

my @anchors=<ANCHOR>;
my @seqs=<SEQ>;
my %results;

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
    push @{$results{$anchors[$min_idx]}->{bySeq}->{$seq}}, ($anchors[$min_idx2],$stability_idx);
    $results{$anchors[$min_idx]}->{byStability}->{$stability_idx}++;
}
foreach my $anchor (keys %results){
    foreach my $stability (keys %{$results{$anchor}->{byStability}}){
        print "$stability\t$results{$anchor}->{byStability}->{$stability}\n";
    }
  }

close ANCHOR;
close SEQ;


#find index of lowest value in array (courtesy of Stack Overflow)
sub get_min_idx {
    my $ra = shift;
    croak "Sub expects array reference" if ref $ra ne 'ARRAY';
    return reduce { $ra->[$a] < $ra->[$b] ? $a : $b } 0..$#$ra;
}




