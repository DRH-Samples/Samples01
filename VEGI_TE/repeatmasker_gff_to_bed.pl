#!/usr/bin/perl -w

# Uses The first "Target" value as the name, for RepeatMasker GFF conversion.
# I tried using Bio::Tools::GFF, but it fails with RepeatMasker GFFs because
# the scores are right-justified using leading spaces, which causes
# Bio::SeqFeature::Generic to error out.

use strict;
use feature qw(say switch);

while (<>) {
    next if (/^#/);
    chomp;
    my ($seq_id, $source, $feature, $start, $end, $score, $strand, $frame, $group)
        = split /\t/;
        
    # BED locations are 0-based
    --$start;
        
    # trim spaces from score, just in case
    $score =~ s/\s+//;
        
    # this is specific to RepeatMasker output
    $group =~ /Target "Motif:(\S+)"/ || die "Couldn't parse repeat name from: $group";
    my $name = $1; 
        
    say join "\t", $seq_id, $start, $end, $name, $score, $strand;
}