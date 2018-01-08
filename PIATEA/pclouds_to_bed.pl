#!/usr/bin/perl -w

# Converts a PClouds ".region" file (which contains genome positions) into a
# BED file. Apparently, when given a multi-FASTA file for input, PClouds
# concatenates all the sequences and reports positions on the concatenated
# sequence. This program maps the positions reported in the ".region" file to
# their positions on the original sequences and reports the results in BED format.

use strict;
use Bio::SeqIO;
#use feature qw (say);

my $usage = "Usage: pclouds_to_bed <PClouds_region_file> <original_FASTA_file>";
die $usage unless 2 == scalar @ARGV;

convert_region_file( $ARGV[0], read_fasta($ARGV[1]) );


# Get the sequence names and lengths from the FASTA file
# Returns array ref.
sub read_fasta {
    my ($fasta_file) = @_;
    my $in = Bio::SeqIO->new(-file => $fasta_file);
    my @seq_data; 
    while (my $seq = $in->next_seq) {
        push @seq_data, [$seq->display_name, $seq->length];                     # can't use hash because need to maintain original order of sequences
    }
    return \@seq_data;
}

sub convert_region_file {
    my ($regions_file, $seq_data_ref) = @_;
    my @seq_data = @$seq_data_ref;
    open IN, $regions_file || die "Couldn't open regions file: $!";
    my $out_file = (split /\./, $regions_file)[0] . ".bed";
    open OUT, ">$out_file" || die "Couldn't open outfile $out_file: $!";
    my $current_seq_name;
    my $current_seq_length;     # length of current sequence
    my $previous_seq_end;       # end of previous sequence on concatenated coordinates
    my $current_seq_end = 0;    # end of current sequence on concatenated coordinates
    my ($start, $end);
    while (<IN>) {
        chomp;
        my ($field1, $field2) = split;
        # It's normal for the last line to be a non-location, but all other lines should be locations
        unless ($field2) {
            last if eof IN;
            die "ERROR: The following line is not a location: $_"; 
        }
        ($start, $end) = ($field1 - 1, $field2);                                # BED locations are 0-based, so need to subtract 1
        unless ($current_seq_name && $start <= $current_seq_end) {
            $previous_seq_end = $current_seq_end;
            ($current_seq_name, $current_seq_length) = @{shift @seq_data};
            $current_seq_end += $current_seq_length;
            #say "processing $current_seq_name";
            print "processing $current_seq_name\n";
        }
        if ($end <= $current_seq_end) {                                         # normal
            #say OUT join "\t", ($current_seq_name, $start - $previous_seq_end, $end - $previous_seq_end);
            print OUT join "\t", ($current_seq_name, $start - $previous_seq_end, $end - $previous_seq_end), "\n"; 
        } else {                                                                # exception: location spans 2 sequences, ignore
            #say "WARNING: $_ goes beyond the end of $current_seq_name ($current_seq_end). Ignoring.";
            print "WARNING: $_ goes beyond the end of $current_seq_name ($current_seq_end). Ignoring.\n";
        }
    }
    close IN;
    close OUT;
    $end <= $current_seq_end ||
        die "ERROR: end of last location ($end) is beyond total length of sequences ($current_seq_end)!";
}