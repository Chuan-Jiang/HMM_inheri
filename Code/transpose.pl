#!/usr/bin/perl 
use warnings; use strict;

### this script used to combine hmm states

open IN, shift @ARGV or die $!; 
my @rc = (-1,-1,-1,-1);
my %win;
my $w = 1000000; 
my @col = ("red","green","yellow","blue","grey");
open BL, ">state_block.txt" or die $!;
while(<IN>){
	chomp;
	my @ar = split;
	if ($ar[0] ne $rc[0] or $ar[2] ne $rc[2]){
		if($rc[0] ne "-1"){
			print BL "$rc[1]\t$col[$rc[2]]\n";
		}
		if ($ar[0] eq $rc[0] and $ar[2] ne $rc[2]){
			print BL "$ar[0]\t$rc[1]\t$ar[1]\t$col[4]\n";
		}
		print BL "$ar[0]\t$ar[1]\t";
	}
	print BL "$ar[1]\t$col[$ar[2]]\n" if (eof(IN));	
	@rc = @ar;

	### window; 
	my $p = (int($ar[1] / $w)) * $w ;
	my @g = split ',', $ar[3];
	#print "$ar[0]\t$p\t$ar[1]\n";
	for (my $i = 0; $i < @g; $i ++){
		$win{$ar[0]}{$p}{$i} += $g[$i];
		$win{$ar[0]}{$p}{s} = $ar[2];
		#print "HH\t${$win{$ar[0]}{$p}}[$i]\t$g[$i]\n";
	}
}

while (my($k,$v) = each %win){
	my %pos_ha = %$v;
	for my $s ( sort {$a <=> $b} keys %pos_ha ){
		my %f = %{$pos_ha{$s}};
		#print "$k\t$s\t@f\n";
		my @sc;
		my @cl;
		my $be = $f{s};
		push @sc,$f{$be};
		push @cl,$be;
		delete $f{$be};
		delete $f{s};

		foreach my $i (sort {$b <=> $a} (keys %f)){
			push @sc,$f{$i};
			push @cl,$i;
		}
		#my $tmp_j = join "\t",@sc,@cl;
		my $ss = $s + $w -1;	
		my $sta = 0;
		for (my $i = 0; $i < @sc; $i ++){
			my $new_sta = $sta + $sc[$i];
			print "$k\t$s\t$sta\t$ss\t$new_sta\t$col[$cl[$i]]\n";
			$sta = $new_sta;
		}
	}
}
