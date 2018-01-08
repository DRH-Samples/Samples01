#!/usr/bin/perl -w
use strict;
use File::Copy 'move';

# Rename each file based on the name given in the fasta header on the first line.
# NOTE: Uses Perl glob and only examines first command-line argument, so must be
# used with double quotes, e.g.:
# rename_files_to_fasta_name.pl "*.fasta"

foreach my $infile (glob $ARGV[0]) {
    open IN, $infile || die "$infile open failed: $!";
    my $first = <IN>;
    close IN;
    $first =~ />(\w+)/ || die "Unable to parse FASTA name from $infile";
    my $outfile = "$1.fasta";
    -e $outfile && die "$outfile already exists";
    move $infile, $outfile || die "move $infile $outfile failed: $!"
}