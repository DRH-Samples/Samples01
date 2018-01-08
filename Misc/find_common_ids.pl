#!/usr/bin/perl -w

# Input: a list of file names
# Output: the common IDs (lines) in each pair of files.
# This script prints common IDs between pairs of files. It is designed to be used 
# with columns that are exported from spreadsheet, so it assumes that each line
# contains one ID. (In other words, it actually just looks for identical lines.)
# Note that it does not check for uniqueness: if an ID, which occurs in file1 at 
# least once, occurs more than once in file2, it will be printed more than once.

use strict;

my $usage = "Usage: find_common_ids.pl file1 file2 [file*]";
die $usage unless scalar @ARGV > 1;

# Compare each unique pair of files
while (my $file1 = shift @ARGV) {
    foreach my $file2 (@ARGV) {
        intersect( $file1, $file2 );
    }
}

sub intersect {
    my ($file1, $file2) = @_;
    
    # Load IDs
    open (FILE1, $file1) or die "Could not open file $file1: $!";
    my @set1 = <FILE1>;
    close FILE1;
    open (FILE2, $file2) or die "Could not open file $file2: $!";
    my @set2 = <FILE2>;
    close FILE2;
    
    # Print intersection
    my %set1Hash;
    foreach my $id1 (@set1) {
        $set1Hash{$id1} = 1;
    }
    print "Intersection between $file1 and $file2:\n";
    foreach my $id2 (sort @set2) {
        print $id2 if $set1Hash{$id2};
    }
    print "\n\n";
}
