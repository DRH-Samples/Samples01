#!/usr/bin/perl -w

# Run FGENESH on one or more sequences.
# Usage:    run_fgenesh seqFile parFile
# Input:    (1) FASTA file with DNA sequence(s), genomic or CDS, and
#           (2) ParFile (e.g. "/usr/local/fgenesh/Dicots").
# Output:   (1) FGENESH output file(s), one per input sequence (xyz.fgenesh),
#           (2) a FASTA file with all predicted aa sequences (xyz.aa), and
#           (3) a FASTA file with all predicted mRNA (CDS) sequences (xyz.cds).


use strict;
use feature qw/ say /;
use File::Basename;
use Bio::SeqIO;
use Bio::Tools::Fgenesh;

my $usage = "Usage: run_fgenesh seqFile parFile";
die $usage unless scalar(@ARGV) == 2;

my $inFile  = $ARGV[0];
my $parFile = $ARGV[1];

my $in = Bio::SeqIO->newFh(-file => "<$inFile", '-format' => 'fasta');
my($inFilename, $inDirectories, $inSuffix) = fileparse($inFile, qr/\.[^.]*/);   # see http://perldoc.perl.org/File/Basename.html
my $aaOut  = makeFastaOut($inFilename, "aa" );
my $cdsOut = makeFastaOut($inFilename, "cds");

while ( my $seq = <$in> ) {
    my $fgeneshOut = runIt($seq);
    parseIt($fgeneshOut, $seq->display_id);
}

# done


###############
# subroutines #
###############

sub makeFastaOut {
    my ($prefix, $suffix) = @_;
    my $outFile = "$prefix.$suffix";
    die "File already exists named $outFile" if (-e $outFile);
    return Bio::SeqIO->new(-file => ">$outFile", '-format' => 'fasta');
}

sub runIt {
    my ($seq) = @_;
    my $seqName = $seq->display_id;
    my $tmpFileName = "$seqName.run_fgenesh.tmp";
    my $seqOut = Bio::SeqIO->new( -file => ">$tmpFileName", -format => 'fasta');
    $seqOut->write_seq($seq);
    $seqOut->close;
    my $fgeneshOut = "$seqName.fgenesh";
    if (-e $fgeneshOut) {
        unlink($tmpFileName)
            || say STDERR "Unable to delete temporary file $tmpFileName";
        die "File already exists named $fgeneshOut";
    }
    my $command = "fgenesh $parFile $tmpFileName -pmrna";
    my $result = `$command`;
    unlink($tmpFileName) || say STDERR "Unable to delete temporary file $tmpFileName";
    open (FGENESH_OUT, ">$fgeneshOut") || die "Couldn't open $fgeneshOut for writing: $!";
    print FGENESH_OUT $result;
    close FGENESH_OUT;
    return $fgeneshOut;
}

sub parseIt {
    my ($fgeneshFile, $seqName) = @_;
    my $fgenesh = Bio::Tools::Fgenesh->new(-file => $fgeneshFile);
    my $gene = $fgenesh->next_prediction();
    die "No predictions found in $fgeneshFile" unless $gene;
    say STDERR "WARNING: More than one prediction found in $fgeneshFile" if $fgenesh->next_prediction();
    my $aa  = $gene->{_predicted_aa};
    my $cds = $gene->{_predicted_cds};
    die "No predicted aa sequence in $fgeneshFile" unless $aa;
    die "No predicted cds sequence in $fgeneshFile" unless $cds;
    $aa->display_id($seqName);
    $cds->display_id($seqName);
    $aaOut->write_seq($aa);
    $cdsOut->write_seq($cds);
}