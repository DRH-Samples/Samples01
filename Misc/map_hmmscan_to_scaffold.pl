#!/usr/local/bin/perl -w
use strict;
use feature qw/ say /;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Tools::GFF;
use DRH::MD::CD::Hmmpfam::TemporaryHmmscanParser;

main();
exit;

sub main {
    if ( scalar(@ARGV) != 3 ) {     # hmmpfam file; CDS GFF file; out file base name
        pod2usage(2);
    }
    my $proteinLocationMap = read_cds( $ARGV[1] );  # proteinID -> CDS location
    my @features = remap_domains( $ARGV[0], $proteinLocationMap );
    my $outFileBaseName = $ARGV[2];
    write_gff(      "$outFileBaseName.gff",     \@features );
    write_gbrowse(  "$outFileBaseName.gbrowse", \@features );
}

sub write_gff {
    my ( $outFileName, $features ) = @_;
    open ( GFF, ">$outFileName" ) or die "Couldn't open $outFileName for writing: $!";
    foreach my $feat (@$features) {
        say GFF $feat->gff_string;
        # stupid $feat->gff_string ignores subfeatures, so do it manually
        foreach my $subFeat ( $feat->get_SeqFeatures ) {
            say GFF $subFeat->gff_string;   
        }
    }
    close GFF;
}

sub write_gbrowse {
    my ( $outFileName, $features ) = @_;
    open ( GBROWSE, ">$outFileName" ) or die "Couldn't open $outFileName for writing: $!";
    
    say GBROWSE '[domain]';
    say GBROWSE 'glyph     = segments';
    say GBROWSE 'connector = dashed';
    say GBROWSE 'bgcolor   = wheat';
    say GBROWSE 'key       = Conserved Domains ("Protein Coding" Genes)';
    say GBROWSE '';
    
    foreach my $feat (@$features) {
        my $model           = get_tag_value( $feat, 'model' );
        my $name            = get_tag_value( $feat, 'name');
        my $desc            = get_tag_value( $feat, 'desc');
        my $interProID      = get_tag_value( $feat, 'interproID');
        my $interProName    = get_tag_value( $feat, 'interproName');
        my $eValue;
        if ( $feat->score =~ /^\D.*\D$/ ) {   # not numeric
            $eValue = $feat->score;
        } else {
            $eValue = sprintf "%1.1G", $feat->score;
        }
        my $source = $feat->source_tag;
        my $topLabel;
        if ( $interProID =~ /NULL/ ) {
            $topLabel = join( ' ', $model, $name, $desc, $eValue, $source );
        } else {
            $topLabel = join( ' ', $model, $name, $desc, $interProID, $eValue, $source );
        }
        $topLabel = "\"$topLabel\"";
        my $bottomLabel = '';
        unless ( $interProName =~ /NULL/ ) {
            my @parts = split /;/, $interProName;
            $bottomLabel = "Note=\"" . $parts[0] . "\"";
        }
        my $loc = to_location_string($feat);
        say GBROWSE join( "\t", 'domain', $topLabel, $loc, $bottomLabel );
        #say GBROWSE;
    }
    close GBROWSE;
}

sub get_tag_value {
    my ( $feat, $tag ) = @_;
    my @arr = $feat->get_tag_values($tag);
    return $arr[0];
}

sub to_location_string {
    my ( $feat ) = @_;
    my $str = $feat->seq_id . ':';
    my @subFeatures = $feat->get_SeqFeatures;
    unless ( $feat->strand ) {
        say STDERR "No strand for feature: " . $feat->gff_string;
    }
    if ( $feat->strand == 1 ) {
        while ( my $eachSubFeature = shift @subFeatures ) {
            $str .= sprintf "%d..%d,", $eachSubFeature->start, $eachSubFeature->end;
        }
    } else {
        while ( my $eachSubFeature = pop @subFeatures ) {
            $str .= sprintf "%d..%d,", $eachSubFeature->start, $eachSubFeature->end;
        }
    }
    my $chopped = chop $str;
    unless ( $chopped eq ',' ) {
        say STDERR "Unexpected last character ($chopped) for feature: " . $feat->gff_string;
        $str .= $chopped;
    }
    return $str;   
}

sub read_cds {
    my ( $gffFile ) = @_;
    my $parser = Bio::Tools::GFF->new( -file => $gffFile, -gff_version => 3 );
    my %proteinLocationMap;
    my $currentProteinID = '';
    my $currentLocation;
    my $feature;
    while( $feature = $parser->next_feature() ) {
        # assume all CDS for each protein are given in a single block
        if ( $feature->primary_tag =~ /CDS/ ) {
            unless ( $feature->has_tag('Parent') ) {
                die "no Parent for: " .  $feature->gff_format;
            }
            my @parents = $feature->get_tag_values('Parent');
            my $proteinID = $parents[0];
            unless ( $proteinID =~ /AT.G\d{5}\.\d/ ) {
                die "unexpected ID $proteinID for: " . $feature->gff_format;
            }
            if ( $proteinID eq $currentProteinID ) {
                $feature->location->seq_id( $feature->seq_id );
                $currentLocation->add_sub_Location( $feature->location );
            } else {
                $proteinLocationMap{ $currentProteinID } = $currentLocation;
                $currentProteinID = $proteinID;
                $currentLocation = new Bio::Location::Split;
                $currentLocation->add_sub_Location( $feature->location );
                $currentLocation->seq_id( $feature->seq_id );
            }
        } else {
            my $foo;
        }
    }
    $parser->close();
    return \%proteinLocationMap;
}

sub remap_domains {
    my ( $hmmscanFileName, $proteinLocationMap ) = @_;
    my $in = new TemporaryHmmscanParser( $hmmscanFileName );
    my @features;
    while ( my $result = $in->next_result ) {
        my $proteinID = $result->query_name;
        unless ( $proteinID =~ /(AT(\d)G\d{5}\.\d+)/ ) {
            say STDERR "unexpected query name: $proteinID";
        }
        my $cdsLocation = $proteinLocationMap->{$proteinID};
        unless ($cdsLocation) {
            say STDERR "No CDS location for $proteinID, skipping.";
            next;
        }
        my $start           = $result->TODO;
        my $end             = $result->TODO;
        my $locationOnProtein = new Bio::Location::Simple( -start => $start, -end => $end );
        
        my $domainName      = $result->TODO;
        my $domainDesc      = $result->TODO;
        my $score           = $result->TODO;
        push @features, transform_domain_location(
            $locationOnProtein, $cdsLocation, $proteinID,
            "hmmscan", $domainName, $domainDesc, $score );
    }
    return @features;
}

# see load_hmmer3_searches_into_db.pl->new_gene_domain_features
sub transform_domain_location {
    my ( $locationOnProtein, $cdsLocation, $proteinID,
         $source, $domainName, $domainDesc, $score ) = @_;
    my @features;
    my $strand = $cdsLocation->strand();
    die "Strand $strand not allowed" unless ( $strand == 1 || $strand == -1 );
    my $dnaLoc;
    my $isMultiExon = 0;
    $dnaLoc = Bio::Location::Split->new();
    my $totalDistToDomain5 = $locationOnProtein->start*3 - 2;   # total bases to the 5' end of the domain start (from the 5' end of the gene)
    my $totalDistToDomain3 = $locationOnProtein->end*3;         # total bases to the 3' end of the domain start
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
        #-frame          => $frame,
#        -display_name   => $domainName,        
        -primary_tag    => 'domain',
        -source         => $source,
        -score          => $score,
    );
    $feat->add_tag_value( 'model',        $proteinID );
    $feat->add_tag_value( 'name',         $domainName );
    $feat->add_tag_value( 'desc',         $domainDesc );
    
    push @features, $feat;
    foreach my $subLoc ( $dnaLoc->sub_Location ) {
        my $subFeat = Bio::SeqFeature::Generic->new (
            -seq_id         => $cdsLocation->seq_id,
            -start          => $subLoc->start,
            -end            => $subLoc->end,
#            -strand         => $subLoc->strand,
#            -display_name   => $domainName,
            -primary_tag    => 'domain',
            -source         => $source,
        );
#        $subFeat->add_tag_value('ID', $domainName.':'.$dnaLoc->start.':'.$dnaLoc->end.':part'); # same ID allowed iff only the location (and score) varies; allows grouping in GBrowse
        $feat->add_SeqFeature( $subFeat ); 
    }

    

    ## Copy any tags from $hsp. Currently this is just 'truncation'.
    #foreach my $tag ( $hsp->get_all_tags ) {
    #    foreach my $value ( $hsp->get_tag_values($tag) ) {
    #        $feat->add_tag_value( $tag, $value);
    #    }
    #}
    
    # Use lowercase names for ad-hoc tags, uppercase are predefined or reserved;
    # See http://www.sequenceontology.org/gff3.shtml .
    #$feat->add_tag_value('Target', $domainName);
    #$feat->add_tag_value('Target', $hsp->hit->start);
    #$feat->add_tag_value('Target', $hsp->hit->end);
    
    #$feat->add_tag_value('Name',   $feat->display_name);   # display_name not printed by Bio::Tools::GFF but is stored in DB
    #$feat->add_tag_value('bits', $hsp->score);
    
    return @features;
}

sub aa_loc_to_dna_loc {
    my ( $aaLoc, $offset ) = @_;
    my $strand = $aaLoc->strand;
    die "Strand $strand not allowed" unless ( $strand == 1 || $strand == -1 );
    my $aaStart   = $aaLoc->start;
    my $aaEnd     = $aaLoc->end;
    my $length    = 3*( $aaEnd - $aaStart + 1 );
    if ($strand == 1) {        
        my $start = $offset + ( 3*$aaStart - 2 );        
        my $end =    $start + ( $length - 1 );
        return Bio::Location::Simple->new(-start => $start, -end => $end, -strand => $strand);
    } else {                
        my $end = $offset - ( 3*$aaStart - 2 );
        my $start  = $end - ( $length - 1 );
        return Bio::Location::Simple->new(-start => $start, -end => $end, -strand => $strand);
    }
}

sub calc_gene_offset {
    my ($loc) = @_;
    if ( $loc->strand == 1 ) {
        return ($loc->start-1);
    } else {
        return ($loc->end+1);
    }    
}

__END__


=head1 NAME

map_domains_to_scaffolds.pl

=head1 SYNOPSIS

map_domains_to_scaffolds.pl Given the locations of domains on CDS, and the location
of CDS on scaffolds, maps the location of the domains on the scaffolds.


=head1 DESCRIPTION

B<map_domains_to_scaffold.pl>

Inputs:
1. The location of domains on CDS, specified by either of two methods:
    a. A TAIR domains file, or
    b. An hmmscan output file.
2. An GFF file containing CDS locations on scaffolds.

Outputs:
1. A gbrowse annotation file that can be uploaded and displayed on gbrowse.
2. A GFF file with the same data (TO DO if needed).

If not all of the genes included in input #1 are included in input #2, the missing
genes are listed to stdout.

E.g. To map the domains of annotated "protein coding" genes on Chr1, input #1
should be the TAIR domains file (which includes all chromosomes) and input #2 should
be the Chr1 NCBI feature table (both of which can be downloaded from arabidopsis.org):

    $ map_domains_to_scaffolds.pl TAIR10_all_domains NCBI_Chr1.tbl 1>chr1_domains.txt 2>err.log

E.g. To map the domains of annotated "transposable element" genes on Chr1, based
on an hmmscan of those genes:
    
    $ map_domains_to_scaffolds.pl chr1_TE_genes-PfamAB-nonempty.hmmscan NCBI_Chr1.tbl 1>chr1_domains.txt 2>err.log


=cut