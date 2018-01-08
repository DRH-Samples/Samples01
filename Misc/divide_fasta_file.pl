#!/usr/bin/perl -w

#
#   Splits a single multi-FASTA file containing sequences for all
#   lyrata scaffolds into 9 files, one each for scaffolds 1-8 and one
#   for scaffolds 9 and above.
#

use strict;
use feature 'say';
use Bio::SeqIO;

my $in_file = $ARGV[0];

my $seq_in = Bio::SeqIO->new(
                            # -fh     => \*STDIN,
                            -file   => $in_file,
                            -format => 'fasta',
                            );

my %seq_table;

while (my $in_seq = $seq_in->next_seq) {
    my $seq_name = $in_seq->display_id;
    $seq_name =~ /.*scaffold(\d+).*/;
    my $scaffold_num = $1;
    unless (defined $scaffold_num) {
        die "unable to parse scaffold number from $seq_name";
    }
    unless ($seq_table{$scaffold_num}) {
        $seq_table{$scaffold_num} = [];
    }
    push @{$seq_table{$scaffold_num}}, $in_seq;
}

foreach my $scaffold_num (sort {$a <=> $b} keys %seq_table) {
    if ($scaffold_num < 9) {
        my $file_name = "AL_scaffold$scaffold_num.fasta";
        my $seq_out = Bio::SeqIO->new(
                            -file   => ">$file_name",
                            -format => 'fasta',
                            );
        say "Writing scaffold $scaffold_num sequences to $file_name";
        foreach my $seq (@{$seq_table{$scaffold_num}}) {
            $seq_out->write_seq($seq)
        }
    }
}

my $file_name = "AL_scaffolds9plus.fasta";
my $seq_out = Bio::SeqIO->new(
                            -file   => ">$file_name",
                            -format => 'fasta',
                            );
foreach my $scaffold_num (sort {$a <=> $b} keys %seq_table) {
    if ($scaffold_num >= 9) {
        say "Writing scaffold $scaffold_num sequences to $file_name";
        foreach my $seq (@{$seq_table{$scaffold_num}}) {
            $seq_out->write_seq($seq)
        }
    }
}
