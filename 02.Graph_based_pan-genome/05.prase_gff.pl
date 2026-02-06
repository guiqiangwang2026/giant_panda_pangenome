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

$in //= 'Ref.gff.gz';
$out //= 'Ref.gff.gene.gz';

my $I = open_in_fh($in);
my $O = open_out_fh($out);

my $gene_now = "";
while(<$I>){
    chomp;
    next if /^#/;
    my @F = split /\t/;
    if ($F[2] eq 'gene'){
        $F[8] =~ /ID=([^;]+)(;|$)/ or die;
        my $gene_id = $1;
        $gene_now = $gene_id;
        $F[8] = "GENE=$gene_id";
    } else {
        die unless defined $gene_now;
        $F[8] = "GENE=$gene_now";
    }
    print $O join("\t", @F), "\n";
}