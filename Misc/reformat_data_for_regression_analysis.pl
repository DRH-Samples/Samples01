#!/usr/bin/perl -w

# Usage:    cat input_CSV_file.csv | perl reformat_data_for_regression.pl > output_CSV_file.csv
# Input:    CSV file of format: "line,Survive 1st,Survive 2nd,Survive 3rd,Dead 1st,Dead 2nd,Dead 3rd"
#           Also, each field except the first must be a positive decimal or the line is ignored,
#           except that '-' or empty fields are treated as '0'.
# Output:   CSV file reformatted for regression analysis, i.e., "line, replicate, [1/0]".

use strict;

print "SALK_ID,Block,Survival\n";
my @wt;
LINE: while (<>) {
    s/\r//g;    # trim carriage returns (chomp doesn't)
    chomp;
    my ($id, $live1, $live2, $live3, $dead1, $dead2, $dead3) = split( /,/ );#, $_, 7 ); # use limit syntax to handle trailing null fields
    if ($id eq 'A') {   
        @wt = ($id, $live1, $live2, $live3, $dead1, $dead2, $dead3);
    } else {
        my $first;
        if (length $id > 8) {
            $first = substr $id, 0, 1;
            $first = 'Z' unless $first =~ /[A-Za-z]/;
            $id = substr $id, -8;
        }
        $first = 'Z' unless $first;
        unless (( substr $id, 0, 1 ) =~ /[A-Za-z]/) {
            $id = $first . substr $id, -7;
        }
        unless ( print_reformatted($id, $live1, $live2, $live3, $dead1, $dead2, $dead3) ) {
            next LINE;
        }
    }
}
die "No w.t." unless @wt;
print_reformatted(@wt);         # print w.t. last, otherwise SAS f-s up
print STDERR "Done.\n";


sub print_reformatted {
    my ($id, $live1, $live2, $live3, $dead1, $dead2, $dead3) = @_;
    foreach my $number ($live1, $live2, $live3, $dead1, $dead2, $dead3) {
        #if (!defined($number) or $number eq '' or $number eq '-') {
        if ($number eq '' or $number eq '-') {
            $number = 0;
        }
        unless ($number =~ /^\d+$/) {
            print STDERR "Invalid line; $number is not a number; skipping: $_\n";
            return 0;
        }
    }
    for (1..$live1) { print "$id,1,1\n"; }
    for (1..$dead1) { print "$id,1,0\n"; }
    for (1..$live2) { print "$id,2,1\n"; }
    for (1..$dead2) { print "$id,2,0\n"; }
    for (1..$live3) { print "$id,3,1\n"; }
    for (1..$dead3) { print "$id,3,0\n"; }
    return 1;
}


