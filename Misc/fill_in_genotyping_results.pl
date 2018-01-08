#!/usr/bin/perl -w
use strict;

# Single use script for Akiko's genotyping table
# to fill in locus ID on each line so the results can sorted by Locus ID using Excel


open(IN, "<$ARGV[0]");

my $locus;
my $line;

while (<IN>) {
    my @columns = split /\t/;
    if ($columns[1]) { 
        $locus = $columns[2];
        $line  = $columns[7];
    } else {
        $columns[2] = $locus;
        $columns[7] = $line;
    }
    print join "\t", @columns;
}