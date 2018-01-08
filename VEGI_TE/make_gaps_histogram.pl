#!/usr/bin/perl -w

# Input:    BED file representing the chain depth track.
# Outputs:  1. Histograms of gap sizes (tab-separated values);
#           2. BED files with gap locations.

use strict;
use feature qw (say);
use List::Util qw (max min);
use Bio::FeatureIO;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Annotated;
#use Data::Dumper;

my $in_file = $ARGV[0];
my @perl_sucks_ass = split /\./, $in_file;
my $file_base_name = $perl_sucks_ass[0];
open IN, "<$in_file" || die "open failed with $!";

# d->lengths, 0<=d<=8, d is the max Chain Depth
# d=0 is gaps of depth 0; d=1 is gaps with max 1, etc.
my %all_gaps;
my @gap_depths = (0..8);

# d->lengths, 1<=d<=9, b is min Chain Depth
# d=9 is blocks of depth 9, d=8 is blocks of depth 9 or 8, etc.
#my %all_blocks;
#my @block_depths = (9..1);

# main
say "Starting...";
do_input();
my %all_gap_lengths = get_lengths(%all_gaps);
#my %all_block_lengths = get_lengths(%all_blocks);
do_output();
say "Done.";
# end

# subs
sub do_input {
    my %gap_starts;
    my %gap_ends;
    my %gap_seq_ids;
    my $current_seq_id = '';
    while (<IN>) {
        chomp;
        my ($seq_id, $start, $end, $depth) = split;
        if ($seq_id ne $current_seq_id) {
            say "Reading $seq_id...";
            $current_seq_id = $seq_id;
        }
        foreach my $limit (@gap_depths) {
            if ( $gap_seq_ids{$limit} and $gap_seq_ids{$limit} ne $seq_id ) {   # exception: end of seq                   
                complete_gap( $limit,
                              $gap_starts{$limit},
                              $gap_ends{$limit},
                              $gap_seq_ids{$limit} );
                $gap_starts{$limit}  = undef;
                $gap_ends{$limit}    = undef;
                $gap_seq_ids{$limit} = undef;
            }
            if ($gap_starts{$limit}) {
                if ($depth <= $limit) {                                         # continue existing gap
                    $gap_ends{$limit} = $end;
                } else {                                                        # complete existing gap
                    complete_gap( $limit,
                                  $gap_starts{$limit},
                                  $gap_ends{$limit},
                                  $gap_seq_ids{$limit} );
                    $gap_starts{$limit}  = undef;
                    $gap_ends{$limit}    = undef;
                    $gap_seq_ids{$limit} = undef;
                }
            } else {                                                            # start new gap
                if ($depth <= $limit) {
                    $gap_starts{$limit} = $start;                                  
                    $gap_ends{$limit} = $end;
                    $gap_seq_ids{$limit} = $seq_id;
                }
            }
        }
    }
}

sub complete_gap {
    my ($limit, $start, $end, $seq_id) = @_;
    my @gap_data = ($seq_id, $start, $end);                                     # data container (gap_data)
    push @{$all_gaps{$limit}}, \@gap_data;
}

# in:   depth_limit->[gap_data]
# out:  depth_limit->[sorted lengths]
sub get_lengths {
    my %all_gap_data = @_;
    say "Sorting lengths...";
    my %all_lengths;
    foreach my $limit (keys %all_gap_data) {
        my @lengths;
        foreach my $gap_data_ar (@{$all_gap_data{$limit}}) {
            my @gap_data = @$gap_data_ar;
            push @lengths, $gap_data[2] - $gap_data[1];                         # no need to add one because input BED is 0-based
        }
        $all_lengths{$limit} = [sort {$a<=>$b} @lengths];
    }
    return %all_lengths;
}

sub do_output {
    output_histograms();
    output_gff();    
}

sub output_histograms {
    my ($raw_histograms_ref, $sliding_histograms_ref) = calculate_histograms();
    say "Writing output...";
    my %raw_histograms = %$raw_histograms_ref;
    my %sliding_histograms = %$sliding_histograms_ref;
    my @output;
    my $max_length = 0;
    foreach my $limit (@gap_depths) {
        $max_length = max($max_length, scalar(@{$raw_histograms{$limit}}));
        push @output, $raw_histograms{$limit};
        push @output, $sliding_histograms{$limit};
    }
    open( OUT, ">$file_base_name-histograms.tsv" );
    say OUT join "\t", "Length", "Depth";
    print OUT "\t";
    for my $d (@gap_depths) {
        print OUT join "\t", $d, "", "";
    }
    say OUT "";
    print OUT "\t";
    for my $d (@gap_depths) {
        print OUT join "\t", "Raw", "Sliding", "";
    }
    say OUT "";
    foreach my $length (1..$max_length) {
        #next if $length%2;                                                      # ignore odd-number rows (because input data is only on even rows)
        my $empty = 1;
        my @row;
        push @row, $length;
        foreach my $c (0..$#output) {
            my $col = $output[$c];
            push @row, defined $col->[$length] ? $col->[$length] : '';
            $empty = 0 if ( $col->[$length] and !($c%2) );                            
        }
        say OUT join "\t", @row unless $empty;
    }
    say OUT "";
    close(OUT);
}

sub calculate_histograms {
    say "Calculating histograms...";
    my $window_size = 50;
    my $sliding_window_max = 50000;                                             # to limit file size and compute time
    my %all_raw_histograms;
    my %all_sliding_histograms;
    foreach my $limit (@gap_depths) {
        my @gaps = @{$all_gap_lengths{$limit}};
        my @raw_histogram = ((0) x ($gaps[$#gaps]+1));
        foreach my $length (@gaps) {
            ++$raw_histogram[$length];
        }
        @{$all_raw_histograms{$limit}} = @raw_histogram;
        #say Dumper(@raw_histogram);
        #say "";
        my $slide_len = min($#raw_histogram, $sliding_window_max);
        my @sliding_histogram = ((0) x $slide_len);
        foreach my $start (0 .. $slide_len) {
            foreach my $position_in_window ($start .. $start + $window_size - 1) {  # window starts at current position
                my $val = $raw_histogram[$position_in_window];
                $sliding_histogram[$start] += $val if ($val);
            }
        }
        @{$all_sliding_histograms{$limit}} = @sliding_histogram;
        #say Dumper(@sliding_histogram);
        #say "";
    }
    return (\%all_raw_histograms, \%all_sliding_histograms);
}

sub output_gff {
    foreach my $limit (@gap_depths) {
        my $out_file = "$file_base_name-gap_depth_$limit.bed";
        open OUT, ">$out_file" || die "Couldn't open $out_file for writing: $!";
        foreach my $gap_data (@{$all_gaps{$limit}}) {
            say OUT join "\t", @$gap_data;
        }
        close OUT;
    }
}
