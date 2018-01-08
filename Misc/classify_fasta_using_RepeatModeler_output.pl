#!/usr/bin/perl -w
use strict;
use feature qw /say/;


# Usage: <program.pl> Final_parsed.txt FASTA
# "Final_parsed.txt" is the summary output file of RepeatModeler.
# FASTA is the file containing the sequences to be labelled.
#
# This script reads the RepeatModeler classifications and adds them to the FASTA
# name using RepeatMasker format. It also removes anything else from the FASTA header
# line (i.e. the description). The new FASTA file is written to STDOUT.

my $classiFile = $ARGV[0];
my $fastaFile  = $ARGV[1];

my %classifications;

open IN1, $classiFile || die "Couldn't open $classiFile: $!";
while (<IN1>) {
    chomp;
    next if /^\s*$/; # ignore blank lines
    my @tuple = split(/\t/);
    my $name = $tuple[0];
    my $class = $tuple[12];
    my $sf = $tuple[21];
    next if $name eq 'RepName'; # skip header
    if ($class eq "dna transposon") {
        $class = "DNA";
    } elsif ($class eq "helitron") {
        #$class = "RC"; # RC does not appear to be properly recognized by RM in the output summary table (counted as Unclassified)
        $class = "DNA";
        $sf = "Helitron";
    } elsif ($class eq "ltr retrotransposon") {
        $class = "LTR";
    } elsif ($class eq "non-ltr retrotransposon") {
        if ($sf eq "sine") {
            $class = "SINE";
            $sf = " ";
        } else {
            $class = "LINE";
        }
    } elsif ($class eq " ") {
            $class = "Unknown";
    } else {
        die "ERROR! Unknown class: $class";
    }
    $sf =~ s/\//-/g;
    $sf = ucfirst($sf);
    $classifications{$name} = [$class, $sf];    # need ref instead?
}

open IN2, $fastaFile || die "Couldn't open $fastaFile: $!";
while (<IN2>) {
    chomp;
    if (/^>([\w-]*)/) {
        my $name = $1;
        unless (exists $classifications{$name}) {
            die "Couldn't find $name";
        }
        my ($class, $sf) = @{$classifications{$name}};
        if ($sf eq " ") {
            say ">$name#$class";
        } else {
            say ">$name#$class/$sf";
        }
    } else {
        say;
    }
}


1;
