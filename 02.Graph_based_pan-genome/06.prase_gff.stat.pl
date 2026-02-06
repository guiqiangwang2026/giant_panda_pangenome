#!/usr/bin/env perl


use strict;
use warnings;
use utf8;
no warnings qw( experimental::smartmatch );
#use v5.24;
use Carp qw/confess carp/; # carp=warn;confess=die

use Data::Dumper;


sub open_in_fh {
    my $file = shift;
    if ($file =~ /\.gz$/) {
        open my $fh, "gzip -dc $file |" or die "Cannot open gzip file $file: $!";
        return $fh;
    } else {
        open my $fh, '<', $file or die "Cannot open file $file: $!";
        return $fh;
    }
}

sub open_out_fh {
    my $file = shift;
    if ($file =~ /\.gz$/) {
        open my $fh, "| gzip > $file" or die "Cannot open gzip file $file for writing: $!";
        return $fh;
    } else {
        open my $fh, '>', $file or die "Cannot open file $file for writing: $!";
        return $fh;
    }
}

my ($in, $out) = @ARGV;

$in //= 'Ref.gff.gene.gz';
$out //= 'Ref.gff.gene.gz.stat';

my $I = open_in_fh($in);
my $O = open_out_fh($out);

my %hash;
while(<$I>){
    chomp;
    next if /^#/;
    my @F = split /\t/;
    my $type = $F[2];
    my $gene = $F[8];
    my $len = $F[4] - $F[3];
    #my $len = $F[4] - $F[3] + 1;
    $hash{$gene}{$type} += $len;
}

foreach my $gene (sort keys %hash){
    foreach my $type (sort keys %{$hash{$gene}}){
        say $O join("\t", $gene, $type, $hash{$gene}{$type});
    }
}