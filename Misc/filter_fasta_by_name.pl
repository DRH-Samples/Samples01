#!/usr/local/bin/perl -w
use strict;
use Pod::Usage;
use Bio::SeqIO;

main();
exit;

sub main {
    if ( scalar(@ARGV) != 2 ) {
        pod2usage(2);
    }
    my $in = new Bio::SeqIO( -file => "<$ARGV[0]" );
    my $out = new Bio::SeqIO( -format => 'fasta' );
    my $filterString = $ARGV[1];
    my @filters = split(/,/, $filterString);
    while ( my $seq = $in->next_seq() ) {
        my $id = $seq->id;
        foreach my $filter (@filters) {
            if ( $id =~ /$filter/i ) {
                $out->write_seq( $seq );
            }
        }
    }
}

__END__


=head1 NAME

filter_fasta_by_name.pl - Reads in a FASTA sequence file and prints out
only the sequences containing the specified string as part of their name
(case insensitive). If more than one substring is specified as a comma-separated
list, returns sequences matching any of the names.


E.g. To make a FASTA file containing AT3G04605 from one containing AT3G04605, AT3G04605_1, and
AT3G04605_2, use something like:

$ filter_fasta_by_name.pl all_genes.fasta "AT3G04605" >AT3G04605.fasta

=head1 SYNOPSIS

filter_fasta_by_header_text.pl fastaFile "Text to search for"


=head1 DESCRIPTION

B<filter_fasta_by_name.pl>

See above.

=cut