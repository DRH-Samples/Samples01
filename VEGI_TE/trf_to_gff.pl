#!/usr/bin/perl -w

# Converts one or more Tandem Repeat Finder output data files
# (-d option from trf execution) to GFF format.

use strict;
use Bio::Tools::TandemRepeatsFinder;
use Bio::FeatureIO;
use Bio::SeqFeature::Annotated;

foreach my $infile (@ARGV) {
    my $parser = Bio::Tools::TandemRepeatsFinder->new(-file => $infile);
    my $seq_id;
    my $out;
    while( my $feature = $parser->next_result ) {
        if ($seq_id) {
            $feature->seq_id() eq $seq_id || die "$infile has results for more than one sequence!";
        } else {
            $seq_id = $feature->seq_id();
            $out = Bio::FeatureIO->new( -format     => 'gff',
                                        -version    => 3,
                                        -file       => ">trf-$seq_id.gff"
                                        );
        }
        my $annotated = Bio::SeqFeature::Annotated->new();
        $annotated->seq_id(         $seq_id);
        $annotated->source_tag(     'trf');
        $annotated->primary_tag(    'tandem_repeat'); # from SOFA, see http://www.sequenceontology.org/resources/intro.html
        $annotated->start(          $feature->start);
        $annotated->end(            $feature->end);
        $annotated->score(          $feature->score);
        $annotated->strand(         0);                 # needed to prevent warnings from gff.pm: "Use of uninitialized value in numeric eq (==)"
        $out->write_feature($annotated);
    }
}
