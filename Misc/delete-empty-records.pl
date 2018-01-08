#!/usr/bin/perl -w
use strict;
use Pod::Usage;


main();
exit;

sub main {
    if ( scalar(@ARGV) != 1 ) {
        pod2usage( "missing filename\n" );
        pod2usage(2);
    }
    my $file = $ARGV[0];
    open IN, "<$file"
        or die "Can't open file '$file': $!\n";
        
        
    # print header
    while ( <IN> ) {
        print;
        last unless /^#/;
    }
        
    # print non-empty records
    my $record = "";
    my $lastLine = "";

    while ( <IN> ) {
        $record .= $_;
        if ( /^\/\// ) {    # must include '^' b/c some descriptions contain "//", e.g. "http://"
            unless ( $record =~  /No hits detected that satisfy reporting thresholds/ ) {
                print $record;
            }
            $record = "";
        }
    }
}


__END__


=head1 NAME

delete-empty-records.pl - Reads in a hmmscan output file and outputs only the
non-empty records contained therein.


=head1 SYNOPSIS

delete-empty-records.pl hmmscan_output_file


=head1 DESCRIPTION

B<delete-empty-records.pl>

See above.

=cut
