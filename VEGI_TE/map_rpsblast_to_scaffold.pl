#!/usr/local/bin/perl -w
use strict;
use feature qw( say );
use List::Util qw( min max);
use Bio::SeqIO;
use Bio::SearchIO;
#use Bio::Tools::GFF;
#use DRH::MD::CD::Hmmpfam::TemporaryHmmscanParser;

my $theScaffold;                                                                # designed to process one scaffold at a time

main();
exit;

sub main {
    my $usage = "Usage: map_rpsblast_to_scaffold.pl <rpsblast_file> <gff_file> max_evalue track_name track_description";
    2 <= scalar(@ARGV) || die $usage;
    my $proteinLocationMap = read_gff( $ARGV[1] );  # proteinID -> CDS location
    my @features = remap_domains( $ARGV[0], $proteinLocationMap, $ARGV[2] );
    #$ARGV[0] =~ /(\S+)\.\S*/ || die "Error parsing file name";
    write_bed( $ARGV[3], $ARGV[4], \@features );
}

# GFF input must look like this:
#AL_scaffold1	FGENESH	gene	5	1147	.	-	.	ID=AL_scaffold1_1
#AL_scaffold1	FGENESH	CDS	5	61	.	-	.	ID=AL_scaffold1_1
#AL_scaffold1	FGENESH	CDS	172	252	.	-	.	ID=AL_scaffold1_1
#AL_scaffold1	FGENESH	CDS	407	802	.	-	.	ID=AL_scaffold1_1
#AL_scaffold1	FGENESH	gene	1248	2534	.	-	.	ID=AL_scaffold1_2
#AL_scaffold1	FGENESH	CDS	1419	1642	.	-	.	ID=AL_scaffold1_2
#AL_scaffold1	FGENESH	CDS	1803	2035	.	-	.	ID=AL_scaffold1_2
#AL_scaffold1	FGENESH	CDS	2124	2365	.	-	.	ID=AL_scaffold1_2

sub read_gff {
    my ($gffFile) = @_;
    my %modelLocations;
    my $currentModel = '';
    open GFF_IN, "<$gffFile" || die "Couldn't open $gffFile for reading: $!";
    while(<GFF_IN>) {
        chomp;
        my @fields = split /\t/;
        my ($scaffold, $primary_tag, $start, $end, $strand) = @fields[0, 2..4, 6];
        $theScaffold || ($theScaffold = $scaffold);
        $scaffold eq $theScaffold || die "Expected $theScaffold: $_";
        $fields[8] =~ /ID=$scaffold[_](\d+)/ || die $fields[8] . " does not match 'ID=[x]'.";
        if ('+' eq $strand) {
            $strand = 1;
        } elsif ('-' eq $strand) {
            $strand = -1;
        } else {
            die "Invalid strand: $_";
        }
        my $model = $1;
        if ('gene' eq $primary_tag) {
            $modelLocations{$model} && die "model $model already processed: $_";
            $currentModel = $model;
            my $location = new Bio::Location::Split;
            $location->seq_id($scaffold);
            $modelLocations{$model} = $location;
        } elsif ('CDS' eq $primary_tag) {
            $model eq $currentModel || die "expecting model $currentModel: $_";
            $modelLocations{$model}->add_sub_Location(
                new Bio::Location::Simple(
                        -seq_id => $scaffold,
                        -start => $start,
                        -end => $end,
                        -strand => $strand ));
        } else {
            die "Primary tag must be 'gene' or 'CDS', not $primary_tag.";
        }
    }
    close GFF_IN;
    return \%modelLocations;
}

sub remap_domains {
    my ( $rpsBlastFile, $proteinLocations, $maxEvalue ) = @_;
    $maxEvalue || ($maxEvalue = 10);
    my $in = new Bio::SearchIO( -file => $rpsBlastFile, -format => 'blast');
    my @features;
    my $saidErr = 0;
    while ( my $result = $in->next_result ) {
        unless ($result->query_name =~ /Model_(\d+)--$theScaffold$/) {
            unless ($saidErr) {
                say STDERR
                    "Not on $theScaffold, ignoring this and all models on other scaffolds: "
                    . $result->query_name;
                $saidErr = 1;
            }
            next;
        }
        my $protein = $1;
        my $cdsLocation = $proteinLocations->{$protein} || die "No CDS location for protein model $protein.";
        my @coveredLocations;                                                   # array of: [start, end, count] 
        while ( my $hit = $result->next_hit ) {
            HSP: while ( my $hsp = $hit->next_hsp ) {
                next if $hsp->evalue > $maxEvalue;
                my $hadOverlap = 0;
                foreach my $covered (@coveredLocations) {
                    my $overlap = get_overlap([$hsp->start, $hsp->end], $covered);
                    if ($overlap) {
                        if ($overlap > 80 && $covered->[2] >= 5) {              # magic numbers: keep up to 5 domains with overlap >80%
                            next HSP;
                        } else {
                            $covered->[0] = min($covered->[0], $hsp->start);
                            $covered->[1] = max($covered->[1], $hsp->end);
                            ++$covered->[2];
                            $hadOverlap = 1;
                            last;
                        }
                    }
                }
                unless ($hadOverlap) {
                    push @coveredLocations, [$hsp->start, $hsp->end, 1];
                }
                push @features, transform_domain_location(
                                    $cdsLocation, $protein, $hit, $hsp );
            }
        }
    }
    return @features;
}

# Locations are: (start, end). Returns the % of the query that overlaps the subject. 
sub get_overlap {
    my ($queryLoc, $subjectLoc) = @_;
    my $overlapStart = max($queryLoc->[0], $subjectLoc->[0]);
    my $overlapEnd   = min($queryLoc->[1], $subjectLoc->[1]);
    if ($overlapStart > $overlapEnd) {
        return 0;
    } else {
        return 100 *
            ($overlapEnd - $overlapStart + 1) /
            ($queryLoc->[1] - $queryLoc->[0] + 1);
    }
}

# see load_hmmer3_searches_into_db.pl->new_gene_domain_features
sub transform_domain_location {
    my ( $cdsLocation, $proteinID, $hit, $hsp ) = @_;
    my @features;
    my $strand = $cdsLocation->strand();
    die "Strand $strand not allowed" unless ( $strand == 1 || $strand == -1 );
    my $dnaLoc;
    my $isMultiExon = 0;
    $dnaLoc = Bio::Location::Split->new( -seq_id => $cdsLocation->seq_id );
    my $totalDistToDomain5 = $hsp->start*3 - 2;                 # total bases to the 5' end of the domain start (from the 5' end of the gene)
    my $totalDistToDomain3 = $hsp->end*3;                       # total bases to the 3' end of the domain start
    my $distToDomain5 = $totalDistToDomain5;                    # remaining bases to the 5' end of the domain start (from the 3' end of the current exon)
    my $distToDomain3 = $totalDistToDomain3;                    # remaining bases to the 3' end of the domain start
    my $cumIntronLength = 0;                                    # cumulative bases in introns up to the current exon
    my ($prevExon3, $domain5, $domain3);
    if ( $strand == 1 ) {           
        my $gene5 = $cdsLocation->start;
        foreach my $currExon ($cdsLocation->sub_Location(1)) {
            my $currExon5 = $currExon->start;
            my $currExon3 = $currExon->end;                     
            if ( defined $prevExon3 ) {
                $cumIntronLength += $currExon5 - $prevExon3 - 1;
            }
            my $subloc = Bio::Location::Simple->new(
                                -seq_id => $cdsLocation->seq_id,
                                -start  => $currExon5,
                                -end    => $currExon3,
                                -strand => $strand );
            $distToDomain5 -= $currExon->length;                
            if ( !defined $domain5 && $distToDomain5 <= 0 ) {
                $domain5 = $gene5 + $totalDistToDomain5 + $cumIntronLength - 1;
                $subloc->start($domain5);
            }
            $distToDomain3 -= $currExon->length;                
            if ( $distToDomain3 <= 0 ) {
                $domain3 = $gene5 + $totalDistToDomain3 + $cumIntronLength - 1;
                $subloc->end($domain3);
                $dnaLoc->add_sub_Location( $subloc );
                last;
            }
            if ( defined $domain5 ) {
                $dnaLoc->add_sub_Location( $subloc );
            }
            $prevExon3 = $currExon->end;                                
        }
    } else {    # strand == -1       
        # >sub_Location(-1): reverse sort order: exons in minus-strand order 
        # (same as gene order); however, exon positions remain on plus strand, start < end
        # differences to strand == 1 method are indicated with "# delta"
        my $gene5 = $cdsLocation->end;                                            
        foreach my $currExon ($cdsLocation->sub_Location(-1)) {   
            my $currExon5 = $currExon->end;                                 
            my $currExon3 = $currExon->start;                               
            if ( defined $prevExon3 ) {
                $cumIntronLength += $prevExon3 - $currExon5 - 1;            
            }                
            my $subloc = Bio::Location::Simple->new(
                                -seq_id => $cdsLocation->seq_id,
                                -start  => $currExon3,                      
                                -end    => $currExon5,                      
                                -strand => $strand );
            $distToDomain5 -= $currExon->length;                
            if ( !defined $domain5 && $distToDomain5 <= 0 ) {
                $domain5 = $gene5 - $totalDistToDomain5 - $cumIntronLength + 1; 
                $subloc->end($domain5);                                     
            }
            $distToDomain3 -= $currExon->length;                
            if ( $distToDomain3 <= 0 ) {
                $domain3 = $gene5 - $totalDistToDomain3 - $cumIntronLength + 1; 
                $subloc->start($domain3);                                   
                $dnaLoc->add_sub_Location( $subloc );
                last;
            }
            if ( defined $domain5 ) {
                $dnaLoc->add_sub_Location( $subloc );
            }
            $prevExon3 = $currExon->start;                                  
        }
    }
    my $feat = Bio::SeqFeature::Generic->new (
        -seq_id         => $cdsLocation->seq_id,
        -start          => $dnaLoc->start,
        -end            => $dnaLoc->end,
        -strand         => $dnaLoc->strand,
        -score          => $hsp->score,
    );
    $feat->add_tag_value( 'model',           $proteinID );
    $feat->add_tag_value( 'cdd_name',        $hit->name );
    $feat->add_tag_value( 'desc',            $hit->description );
    $feat->add_tag_value( 'length',          $hsp->length );
    $feat->add_tag_value( 'evalue',          $hsp->evalue );
    $feat->add_tag_value( 'frac_identical',  $hsp->frac_identical );
    $feat->add_tag_value( 'frac_conserved',  $hsp->frac_conserved );
    $feat->add_tag_value( 'query_string',    $hsp->query_string );
    $feat->add_tag_value( 'subject_string',  $hsp->hit_string );
    $feat->add_tag_value( 'homology_string', $hsp->homology_string );
    $feat->add_tag_value( 'query_start',     $hsp->start('query') );
    $feat->add_tag_value( 'query_end',       $hsp->end(  'query') );
    $feat->add_tag_value( 'subject_start',   $hsp->start ('subject' ) );
    $feat->add_tag_value( 'subject_end',     $hsp->end( 'subject' ) );
    
    push @features, $feat;
    foreach my $subLoc ( $dnaLoc->sub_Location ) {
        my $subFeat = Bio::SeqFeature::Generic->new (
            -seq_id         => $cdsLocation->seq_id,
            -start          => $subLoc->start,
            -end            => $subLoc->end,
            -strand         => $subLoc->strand,
        );
        $feat->add_SeqFeature( $subFeat ); 
    }
    return @features;
}

sub write_bed {
    my ( $trackName, $trackDescription, $features ) = @_;
    my $scaffold_prefix = "AL_scaffold";                                        # good enough for now
    say "track name=$trackName type=bedDetail description=\"$trackDescription\"";               
    foreach my $feat (@$features) {
        $feat->seq_id =~/$scaffold_prefix(\d+)/ || die "Error parsing scaffold from" . $feat->seq_id;
        my $scaffold = "scaffold_$1";
        my $strand = $feat->strand == 1 ? '+' : '-';
        my $score = $feat->score;
        $score = min ($score, 1000);
        my $desc = get_tag_value($feat, 'desc');
        $desc =~ /(\S+), (\S+), ([^.]+)\. ?(.*)/ || die "Couldn't parse description: $desc";
        my ($name1, $name2, $name3, $description) = ($1, $2, $3, $4);
        say join "\t",
            $scaffold,                                                          # chrom
            $feat->start-1,                                                     # chromStart
            $feat->end,                                                         # chromEnd
            $name1 eq $name2 ? $name1 : "$name1:$name2",                         # name
            $score,                                                             # score
            $strand,                                                            # strand
            $feat->start-1,                                                     # thickStart
            $feat->end,                                                         # thickEnd
            "0,255,0",                                                          # itemRGB (green)
            getBlockCountSizesAndStarts($feat),                                 # blockCount, blockSizes, blockStarts
            $name1,
            get_details($feat, $name1, $name2, $name3, $description);
    }
}

sub get_details {
    my ($feat, $name1, $name2, $name3, $description) = @_;
    
    my $details = sprintf "<BR><a href=\"http://www.ncbi.nlm.nih.gov/cdd?term=%s\">%s</a>, ",
                    $name1,
                    get_tag_value($feat, 'cdd_name');
    if ($name1 =~ /pfam(\d+)/) {
        $details .= sprintf "<a href=\"http://pfam.sanger.ac.uk/family/PF%s\">%s</a>, ",
                    $1,
                    $name1;
    } else {
        $details .= sprintf "%s, ",
                    $name1;
    }
    $details .= sprintf "%s, %s.<BR><BR> %s<BR>",
                    $name2,
                    $name3,
                    $description;
    $details .= "<CODE><PRE>";
    $details .= sprintf "Expect = %s, Identities = %.1f%%, Positives = %.1f%%, Length = %s <BR><BR>",
                    get_tag_value($feat, 'evalue'),
                    get_tag_value($feat, 'frac_identical')*100,
                    get_tag_value($feat, 'frac_conserved')*100,
                    get_tag_value($feat, 'length');
    $details .= sprintf "Gene    %3d  %s  %3d <BR>",
                    get_tag_value($feat, 'query_start'),
                    get_tag_value($feat, 'query_string'),
                    get_tag_value($feat, 'query_end');
    $details .= sprintf "             %s <BR>",
                    get_tag_value($feat, 'homology_string');
    $details .= sprintf "Domain  %3d  %s  %3d <BR>",
                    get_tag_value($feat, 'subject_start'),
                    get_tag_value($feat, 'subject_string'),
                    get_tag_value($feat, 'subject_end');
    $details .= "</PRE></CODE>";
    return $details;
}

sub getBlockCountSizesAndStarts {
    my ($feat) = @_;
    my $count = 0;
    my @sizes;
    my $size = 0;
    my @genomic_starts;
    foreach my $subFeat ( $feat->get_SeqFeatures ) {                            
        ++$count;
        $size = $subFeat->length;
        push @sizes, $size;
        push @genomic_starts, $subFeat->start;
    }
    if (-1 == $feat->strand) {                                                  # reverse order for minus strand
        @sizes = reverse(@sizes);
        @genomic_starts = reverse(@genomic_starts);
    }
    my @block_starts = map {$_ - $feat->start} @genomic_starts;
    return join "\t",
                $count,                                                         # blockCount
                (join ",", @sizes),                                             # blockSizes
                (join ",", @block_starts);                                      # blockStarts
}

sub get_tag_value {
    my ( $feat, $tag ) = @_;
    my @arr = $feat->get_tag_values($tag);
    return $arr[0];
}