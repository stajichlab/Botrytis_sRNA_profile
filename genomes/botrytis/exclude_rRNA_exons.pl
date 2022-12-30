#!/usr/bin/env perl
use strict;
use warnings;

my $in = shift || 'BcinereaB05-10.features.sort.wRMrRNA.gff3';

open(my $fh => $in) || die "Cannot open $in. $!";
my %rRNAgenes;

while(<$fh>) {
	next if /^#/;
	my @row = split(/\t/,$_);
	if ($row[2] eq 'rRNA') {
		my %group = map { split(/=/,$_) } split(';',pop @row);
		$rRNAgenes{$group{'ID'}}++;
	}
}
# reopen unless we just call seek instead
open($fh => $in) || die "Cannot open $in. $!";
while(<$fh>) {
	if ( /^#/ ) {
		print;
		next
	}
        my @row = split(/\t/,$_);
        if ($row[2] eq 'exon') {
                my %group = map { split(/=/,$_) } split(';',pop @row);
		next if (exists $rRNAgenes{$group{'Parent'}} );
        }
	print $_
}
