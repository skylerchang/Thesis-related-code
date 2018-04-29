#!/usr/bin/perl

use strict;use warnings;
use List::Util qw(reduce);
use Carp;

open ANCHOR,"<$ARGV[0]";
open SEQ,"<$ARGV[1]";

my @anchors=<ANCHOR>;
my @seqs=<SEQ>;
my %results;
my %seqs_index = map { $seqs[$_] => $_ } reverse 0 .. @seqs-1;
foreach my $seq (keys %seqs_index){
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
    my $stability_idx = $min_val / $min_val2;
    if ($stability_idx == 1) {
        push @{$results{$seqs_index{$seq}}->{Coordinate}},"A";
    }else {
        push @{$results{$seqs_index{$seq}}->{Coordinate}}, $min_idx;
    }
}

foreach my $seq_index (sort {$a<=>$b}keys %results){
    foreach my $coordinate ( %{$results{$seq_index}}->{Coordinate}->[0]){
        print "$coordinate\n";
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




