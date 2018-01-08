#!/usr/bin/perl -w

# Compare 2 sequence files. Prints out a table of sequence names that have
# identical sequences between the 2 files, and those with no matching sequences.
# If the same sequence is found multiple times in one file but not the other,
# only the first occurrence will be found.

use strict;
use feature qw/ say /;
use Bio::SeqIO;

my $usage = "Usage: diff_seqfiles seqFileA seqFileB";
die $usage unless scalar(@ARGV) == 2;

my $inFileA = $ARGV[0];
my $inFileB = $ARGV[1];

my $inA = Bio::SeqIO->newFh(-file => "<$inFileA");
my $inB = Bio::SeqIO->newFh(-file => "<$inFileB");

my %seqsA;
my %seqsB;
my %aToB;

while ( my $seq = <$inA> ) {
    my $seqName = $seq->display_id;
    if ( exists $seqsA{$seqName} ) {
        say STDERR "WARNING: duplicate sequence name $seqName in $inFileA";
        next;
    }
    $seqsA{$seqName} = $seq;
}

while ( my $seq = <$inB> ) {
    my $seqName = $seq->display_id;
    if ( exists $seqsB{$seqName} ) {
        say STDERR "WARNING: duplicate sequence name $seqName in $inFileB";
        next;
    }
    $seqsB{$seqName} = $seq;
}

foreach my $seqA ( values %seqsA ) {
    foreach my $seqB ( values %seqsB ) {
        if ( $seqA->seq eq $seqB->seq ) {
            $aToB{$seqA->display_id} = $seqB->display_id;
            delete $seqsA{$seqA->display_id};
            delete $seqsB{$seqB->display_id};
            last;
        }
    } 
}

say "\nIdentical sequences:\n";
foreach my $nameA ( sort keys %aToB ) {
    say sprintf "%-10s = %-s", $nameA, $aToB{$nameA};
}

say "\nNot found from $inFileA:\n";
foreach my $name ( sort keys %seqsA ) {
    say $name;
}

say "\nNot found from $inFileB:\n";
foreach my $name ( sort keys %seqsB ) {
    say $name;
}

say "\nDone.";


