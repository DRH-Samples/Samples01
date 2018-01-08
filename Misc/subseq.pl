#!/usr/local/bin/perl -w
use strict;
use feature qw[say];
use Getopt::Long;
use Pod::Usage;
use List::Util qw[min max];
use Bio::SeqIO;
use Bio::Seq;

my $file;
my $seqName;
my $start;
my $end;
my $flank = 0;
my $reverseComplement = 0;

parseCommandLine();
my $in = new Bio::SeqIO(-file => $file);
my $seq;
while (my $eachSeq = $in->next_seq) {
    if ($eachSeq->id eq $seqName) {
        $seq = $eachSeq;
        last;
    }
}
die "$seqName not found" unless $seq;
$start = max( $start - $flank, 1 );
$end   = min( $end + $flank, $seq->length );
my $subSeq;
unless ($reverseComplement) {
    $subSeq = Bio::Seq->new( -display_id => "$seqName:$start..$end",
                             -seq        => $seq->subseq($start,$end) );
} else {
    $subSeq = Bio::Seq->new( -display_id => "$seqName:$end..$start",
                             -seq        => $seq->subseq($start,$end) );
    $subSeq = $subSeq->revcom;
}
Bio::SeqIO->new(-format => 'fasta')->write_seq($subSeq);

sub parseCommandLine {
    my $help = 0;
    GetOptions(     'help'          => \$help,
                    'input=s'       => \$file,
                    'name=s'        => \$seqName,
                    'start=i'       => \$start,
                    'end=i'         => \$end,   
                    'flank=i'       => \$flank,
                    'reverse'       => \$reverseComplement,
              )
      || pod2usage(2);
    pod2usage(2) if $help;
    unless ( defined $file      &&
             defined $seqName   &&
             defined $start     &&
             defined $end   
            )
    {
        pod2usage( "missing required parameter\n" );
        pod2usage(2);
    }
}


__END__


=head1 NAME

subseq.pl - Prints the specified subsequence plus flanking sequence (optional).

=head1 SYNOPSIS

subseq.pl -i seq_file -n seq_name -s start -e end [-f flank]

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exit.

=item B<-man>

Print the manual page and exit.

=item S<-input>

File containing the original sequence(s).

=item S<-name>

Sequence name (in case there is more than one sequence in the file).

=item I<-start>

Subsequence start.

=item I<-end>

Subsequence start.

=item I<-flank>

Flanking sequence length (optional).

=back

=cut