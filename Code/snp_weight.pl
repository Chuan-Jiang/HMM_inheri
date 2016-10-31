#!/usr/bin/perl
use warnings; use strict;
### snp weight in windows 
my($snp_f,$block_f) = @ARGV;
my $win = 1000000;

open SNP, $snp_f or die $!;

my %ha;
while(<SNP>){
	chomp;
	my @ar = split ;
	my $chr = $ar[0];
	my $p = $ar[1];
	my @ty = @ar[3,4,5,6];
	my $k = int($p/$s);
	for (my $i = 0 ; $i <4;$i ++){
		my $w;
		if ($ty[$i] < 0.45){
			$w = 0;
		}elsif($ty[$i] >= 0.45 and $ty[$i] < 1){
			$w = 0.5;
		}else{
			$w = 1;
		}
		$ha{$chr}{$k}[$i] += $w;
	}
}

while (my($k,$v)= each %ha){
	my %sh = %$v;
	foreach my $p (sort {$a <=>$b} keys %sh){
		my @a = @{$sh{$p}};
		my $na = join "\t",@a;
		print "$k\t$p\t$na\n";
	}
}


