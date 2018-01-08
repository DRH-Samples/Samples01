#!/usr/local/bin/perl -w
use strict;
use feature qw(say);
use Log::Log4perl qw(:easy);
use Bio::Location::Simple;
use Bio::Factory::FTLocationFactory;
use Bio::DB::SeqFeature::Store;
use Bio::Tools::GFF;
use Bio::Perl;
use Data::Dumper;
use Cwd;
use DRH::MD::CD::Hmmpfam::TemporaryHmmscanParser;

Log::Log4perl->easy_init($DEBUG);

my $gffio = Bio::Tools::GFF->new( -gff_version => 3 );

my $chrAttributeTable = {
    1 =>    ['NC_003070', 30427671], 
    2 =>    ['NC_003071', 19698289],
    3 =>    ['NC_003074', 23459830],
    4 =>    ['NC_003075', 18585056],
    5 =>    ['NC_003076', 26975502],
};

my $RFA_TABLE =
{
    1 => { 0 => undef, 1 => undef, 2 => undef},
    2 => { 0 => undef, 1 => undef, 2 => undef},
    3 => { 0 => undef, 1 => undef, 2 => undef},
    4 => { 0 => undef, 1 => undef, 2 => undef},
    5 => { 0 => undef, 1 => undef, 2 => undef},                       
};

my $DB= Bio::DB::SeqFeature::Store->new(
                                        -adaptor => "DBI::mysql",
                                        -dsn     => 'gbrowse_test2',
                                        -user    => 'nobody',
                                    );

main();
exit;

sub main {
    #my $inFileGlob = "/Users/dhoen/Work/Data/Arabidopsis/TAIR9/Custom/Conserved_Domains/TE_Genes/chr1_TE_genes-PfamAB-nonempty.hmmscan";
    my $inFileGlob = "/Users/dhoen/Work/Data/Arabidopsis/TAIR9/Custom/Conserved_Domains/TE_Genes/*-nonempty.hmmscan";
    my @featuresByChr = load_all_genes( $inFileGlob );
    foreach my $features (@featuresByChr) {
        write_gbrowse($features);
    }
}

sub write_gbrowse {
    my ( $features ) = @_;
    my $chr = @$features[0]->seq_id;
    my $outFileName = $chr . "_te_gene_domains.gbrowse";
    open ( GBROWSE, ">$outFileName" ) or die "Couldn't open $outFileName for writing: $!";
    say GBROWSE '[domain]';
    say GBROWSE 'glyph     = segments';
    say GBROWSE 'connector = dashed';
    say GBROWSE 'bgcolor   = wheat';
    say GBROWSE 'key       = Conserved Domains ("Transposable Element" Genes)';
    say GBROWSE '';
    
    foreach my $feat (@$features) {
        my $locus           = get_tag_value( $feat, 'locus' );
        my $name            = get_tag_value( $feat, 'Target');
        my $truncation      = get_tag_value( $feat, 'truncation');
        my $eValue;
        if ( $feat->score =~ /^\D.*\D$/ ) {   # not numeric
            $eValue = $feat->score;
        } else {
            $eValue = sprintf "%1.1G", $feat->score;
        }
        my $topLabel = join( ' ', $name, $eValue, $truncation );
        $topLabel = "\"$topLabel\"";
        my $bottomLabel = "Note=\"" . $locus . "\"";
        my $loc = to_location_string($feat);
        #my $loc = sprintf "%s:%d..%d", $feat->seq_id, $feat->start, $feat->end;
        say GBROWSE join( "\t", 'domain', $topLabel, $loc, $bottomLabel );
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

sub load_all_genes {
    my ( $inFileGlob ) = @_;
    my @featuresByChr;
    foreach my $file ( glob( $inFileGlob ) ) {
        say "Processing file $file";
        my @features = load_genes($file);
        push @featuresByChr, \@features;
    }
    return @featuresByChr;
}

# E.g. Query:       AT2G01022_2  [L=794]
sub load_genes {
    my $file = shift;
    my $in = new TemporaryHmmscanParser( $file );
    my $theChrNum;
    my @features;
    while ( my $result = $in->next_result ) {
        my ($locus, $chrNum, $frame) = parse_query_name($result->query_name);    # $frame not defined if == 0
        if (defined $theChrNum) {
            die 'one chromosome at a time' unless $chrNum == $theChrNum;
        } else {
            $theChrNum = $chrNum;
        }
        my $exons = parse_query_description($result->query_description);
        while ( my $hit = $result->next_hit() ) {
            while ( my $hsp = $hit->next_hsp() ) {
                my $feat = new_gene_domain_features($hsp, $hit->name, $locus, $exons, $chrNum, $frame);                
                push @features, $feat;
            }
        }
    }
    my $outFileName = "te_gene_domains-chr$theChrNum.gff3";
    store_features($outFileName, @features);
    return @features;
}


# if it's a:
# 1) chromosome: returns chrNum, strand, frame
# 2) ORF: returns chrNum, strand, frame, start, end
# 3) gene: returns locus, chrNum, (frame)
#
# frame is only for single-exon pseudogenes which were translated into 3
# frames because it isn't clear if they should contain
# multiple exons but were overlooked b/c of their status (pseudo which may simply
# mean TE-related)
sub parse_query_name {
    my ($queryName) = @_;
    # chr or ORF:
    if ( $queryName =~ /Chr([1-5])_([1-6])(?::(\d+)\.\.(\d+))?/ ) {     # '?:' is to make the '()' a non-grouping group
        # $1=chrNum; $2=frame/strand token; $3=orf_start; 4=$orf_end
        my ($strand, $frame) = parse_strand_and_frame($2);
        return ( $1, $strand, $frame, $3, $4 );
    # gene:
    } elsif ( $queryName =~ /(AT(\d)G\d{5})(?:_([1-2]))?/ ) {
        return ($1,$2,$3);
    } else {
        die 'Error parsing query name ' . $queryName;
    }    
}

# location may be a Simple location , if there is one exon, or a Split location
# if there are multiple exons; note that exons may be a single base, e.g. AT1G79740:
# complement(join(30004367..30004767,30005074..30006627,30006715))[description...]
sub parse_query_description {
    my ($queryDescription) = @_;
    if ( $queryDescription =~
        /^((?:(?:join|complement|complement\(join)\()?(?:\d+(?:\.\.\d+)?,?)+\)?\)?)/ ) {
            my $loc = Bio::Factory::FTLocationFactory->from_string($1);
            return $loc;
    } else {
        die "Can't parse location from $queryDescription"
    }
}


sub store_features {
    my ($fileName, @features) = @_;
    #if ( -e $fileName && !$OVERWRITE ) {
    #    die getcwd() . "/$fileName exists";
    #}
    open( my $out, '>', $fileName ) or die "Couldn't open $fileName for writing: $!";
    say $out '##gff-version 3';
    foreach my $feat (@features) {
        say $out $gffio->gff_string($feat);
    }
    say "wrote ", scalar(@features), " features to $fileName.";    
}

# these are the tokens I used in the files to indicate strand + frame (an unfortunate choice)
# token 1,2,3 -> strand +1, frame 0,1,2
# token 4,5,6 => strand -1, frame 0,1,2
sub parse_strand_and_frame {
    my $token = shift;
    die "token must be in range 1-6" if $token < 1 or $token > 6;
    if ( $token < 4 ) {
        return ( 1, $token-1 );
    } else {
        return ( -1, $token-4 );
    }
}

# domain feature from gene-based search
# frame is optional
# Use the location info stored in query name (passed in as $exons) to
# determine exon positions because
# there are multiple gene models and I (arbitrarily) chose the ".1" model
# (which is the longest one) to use in my searches; it won't be easy to
# tease this back out of the database afaik.
#
# For now, just use this lookup to verify the position; specifically the lookup
# gives the gene position, which should always contain the exon positions
# (and will be identical in the case of pseudogenes).
sub new_gene_domain_features {
    my ($hsp, $domainName, $locus, $exons, $chrNum, $frame) = @_;
    $frame = 0 unless defined $frame;
    my $chrAccession = $chrAttributeTable->{$chrNum}[0];
    my @lookup = $DB->features(
                        -type => ['gene', 'pseudogene'],
                        -seq_id => $chrAccession,
                        -attribute => {locus_tag => $locus} 
                     );
    unless ( @lookup == 1 ) {
        die "Expected exactly 1 feature in database for $locus but found " . scalar(@lookup);
    }
    my $geneFeat = $lookup[0];
    unless ( $geneFeat->contains($exons)) {
        die "$locus is not in the expected position"
    }
    my $strand = $exons->strand();
    die "Strand $strand not allowed" unless ( $strand == 1 || $strand == -1 );
    my $dnaLoc;
    my $isMultiExon = 0;
    $dnaLoc = Bio::Location::Split->new();
    # note: hsp locations are oriented with the gene not the chromosome and
    # are in aa-coordinates not dna-coordinates
    my $totalDistToDomain5 = $hsp->start*3 - 2; # total bases to the 5' end of the domain start (from the 5' end of the gene)
    my $totalDistToDomain3 = $hsp->end*3;       # total bases to the 3' end of the domain start
    my $distToDomain5 = $totalDistToDomain5;    # remaining bases to the 5' end of the domain start (from the 3' end of the current exon)
    my $distToDomain3 = $totalDistToDomain3;    # remaining bases to the 3' end of the domain start
    my $cumIntronLength = 0;                    # cumulative bases in introns up to the current exon
    my ($prevExon3, $domain5, $domain3);
    my @locations;
    if ( $exons->isa('Bio::Location::Simple') ) {
        push @locations, $exons;
    } else {
        @locations = $exons->sub_Location(1);
    }
    if ( $strand == 1 ) {           
        my $gene5 = $exons->start;
        # does not appear necessary to specify sort order '1' with current data,
        # but just in case (e.g. for different searches or parsers), since if not
        # specified the order is apparently just the order that they were added
        foreach my $currExon ( @locations ) {
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
        my $gene5 = $exons->end;                                            # delta
        foreach my $currExon ( @locations ) {   
            my $currExon5 = $currExon->end;                                 # delta
            my $currExon3 = $currExon->start;                               # delta
            if ( defined $prevExon3 ) {
                $cumIntronLength += $prevExon3 - $currExon5 - 1;            # delta
            }                
            my $subloc = Bio::Location::Simple->new(
                                -start  => $currExon3,                      # delta
                                -end    => $currExon5,                      # delta
                                -strand => $strand );
            $distToDomain5 -= $currExon->length;                
            if ( !defined $domain5 && $distToDomain5 <= 0 ) {
                $domain5 = $gene5 - $totalDistToDomain5 - $cumIntronLength + 1; # delta
                $subloc->end($domain5);                                     # delta
            }
            $distToDomain3 -= $currExon->length;                
            if ( $distToDomain3 <= 0 ) {
                $domain3 = $gene5 - $totalDistToDomain3 - $cumIntronLength + 1; # delta
                $subloc->start($domain3);                                   # delta
                $dnaLoc->add_sub_Location( $subloc );
                last;
            }
            if ( defined $domain5 ) {
                $dnaLoc->add_sub_Location( $subloc );
            }
            $prevExon3 = $currExon->start;                                  # delta
        }
    }
    my $aaSeq = $hsp->query_string;
    #unless ( verify_location($dnaLoc, $aaSeq, $chrNum) ) {
    #    die "Location calculation error";
    #}
    
    my $feat = Bio::SeqFeature::Generic->new (
        -seq_id         => "Chr$chrNum",
        -start          => $dnaLoc->start,
        -end            => $dnaLoc->end,
        -strand         => $dnaLoc->strand,
        -frame          => $frame,
        -display_name   => $domainName,        #"$hitName:$chrAccession:" . $loc->start . ':' . $loc->end,
        -primary_tag    => 'domain',
        -source         => 'hmmpfam',
        -score          => $hsp->significance,  # use this instead of the bit score because it's more meaningful
    );
    $feat->add_tag_value( 'locus',  $locus);
    $feat->add_tag_value( 'name',   $domainName );
    $feat->add_tag_value( 'bits',   $hsp->score);
    my $query = uc $hsp->query_string;
    $query =~ s/-//g;
    $feat->add_tag_value('query', $query);
    
    $feat->add_tag_value('Target', $domainName);
    $feat->add_tag_value('Target', $hsp->hit->start);
    $feat->add_tag_value('Target', $hsp->hit->end);
    
    # Copy any tags from $hsp. Currently this is just 'truncation'.
    foreach my $tag ( $hsp->get_all_tags ) {
        foreach my $value ( $hsp->get_tag_values($tag) ) {
            $feat->add_tag_value( $tag, $value);
        }
    }  
    foreach my $subLoc ( $dnaLoc->sub_Location ) {
        my $subFeat = Bio::SeqFeature::Generic->new (
            -seq_id         => "Chr$chrNum",
            -start          => $subLoc->start,
            -end            => $subLoc->end,
            -primary_tag    => 'domain',
            -source         => 'hmmpfam',
        );
        $feat->add_SeqFeature( $subFeat ); 
    }
    return $feat;
}

# hsp must be from a translated chromosome search
# only needed if these values have not yet been computed (e.g. from a previous run)
sub record_reverse_frame_adjustment {
    my ( $chrNum, $frame, $hsp ) = @_;
    print 'Calculating reverse frame adjustment... ';
    die "\nAlready calculated" if $RFA_TABLE->{$chrNum}->{$frame};            
    my $aaLoc = Bio::Location::Simple->new( -start => $hsp->start, -end => $hsp->end, -strand => -1 );
    my $aaSeq = $hsp->query_string; 
    for ( my $adjust = -1; $adjust < 4; ++$adjust ) {
        my $offset = calc_chr_offset($chrNum, -1);  # minus strand
        my $dnaLoc = transform_to_dna_loc( $aaLoc, $frame, $offset, $adjust );
        if ( verify_location( $dnaLoc, $aaSeq, $chrNum ) ) {
            $RFA_TABLE->{$chrNum}->{$frame} = $adjust;
            say $adjust;
            return;
        }
    }
    die "\nError in location calculation.";  # locAdjust should not exceed 4 (afaik)    
}

# The offset is the distance from the start of the chromosome to the 'beginning' of
# the feature, where 'beginning' is the start position if it's on the plus strand
# and the end position if it's on the minus strand. This is to allow locations to
# be calculated just knowing the AA-coordinates of a sub-feature (e.g. a domain). 
# The offset is 0 if the start of the feature is at the first base of the chromosome 
# on the plus strand, or at the last base of the chromosome on the minus strand.
#
# Offset calculations are different for orf-based and gene-based searches because,
# for the minus strand, the start of the ORF is measured from the end of the chromosome,
# since it is the result of translated complement of the chromosome DNA sequence,
# whereas that of the gene is measured from the beginning. (Also, ORF-based positions
# must be converted from AA-coordinates to DNA-coordinates.)
sub calc_orf_offset {
    my ($chrNum, $loc, $frame) = @_;
    if ( $loc->strand == 1 ) {
        return 3*($loc->start-1);
    } else {
        my $chrLength = $chrAttributeTable->{$chrNum}[1];        
        my $offset = $chrLength - 3*($loc->start-1);
        return $offset;
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

# used for chromosome-based searches; needed because the 'beginning' of features
# on the minus strand is the end of the chromosome
sub calc_chr_offset {
    my ($chrNum, $strand) = @_;
    if ( $strand == 1 ) {
        return 0;
    } else {
        my $chrLength = $chrAttributeTable->{$chrNum}[1];        
        return $chrLength;
    }    
}

# adjust: reverse-strand frame adjustment
sub transform_to_dna_loc {
    my ( $aaLoc, $frame, $offset, $adjust ) = @_;
    $adjust = 0 unless defined $adjust;
    my $strand = $aaLoc->strand;
    die "Strand $strand not allowed" unless ( $strand == 1 || $strand == -1 );
    my $aaStart   = $aaLoc->start;
    my $aaEnd     = $aaLoc->end;
    my $length    = 3*( $aaEnd - $aaStart + 1 );
    if ($strand == 1) {        
        my $start = $offset + ( 3*$aaStart - 2 + $frame );        
        my $end =    $start + ( $length - 1 );
        return Bio::Location::Simple->new(-start => $start, -end => $end, -strand => $strand);
    } else {                
        #$end    = $chrLength - ( $offset + 3*$aaStart - 2 + $frame ) + $adjust;      # reverse strand but still set start < end
        my $end = $offset - ( 3*$aaStart - 2 + $frame ) + $adjust;
        my $start  = $end - ( $length - 1 );
        return Bio::Location::Simple->new(-start => $start, -end => $end, -strand => $strand);
    }
}


sub verify_location {
    my ($loc, $aaSeq, $chrNum) = @_;
    my $accession = $chrAttributeTable->{$chrNum}[0];

    
    # TODO: make this work if $loc is a SplitLocation:
    my $dna = '';
    if ( $loc->isa( 'Bio::Location::Simple' ) ) {
        $dna = $DB->fetch_sequence($accession, $loc->start, $loc->end);
    } else {
        # should be a Bio::Location::Split
        # sub_Location(1): guarantee ordering with increasing position on + strand 
        foreach my $subLoc ( $loc->sub_Location(1) ) {
            $dna .= $DB->fetch_sequence($accession, $subLoc->start, $subLoc->end);
        }
    } 
    
    
    if ( $loc->strand == -1 ) {
        $dna = Bio::Perl::reverse_complement_as_string($dna);
    }
    if ( $dna eq '' ) {
        die "Unable to fetch sequence $accession:" . $loc->start . '-' . $loc->end .
        ". Are you sure a sequence exists in the database with ID $accession" .
        " and that it is long enough?";
    }
    my $aa = Bio::Perl::translate_as_string($dna);
    my $expected = uc $aaSeq;
    $expected =~ s/-//g;
    unless ( $aa eq $expected ) {
        my $debugBreakpointHere;
    }
    return ($aa eq $expected);
}

1;
__END__