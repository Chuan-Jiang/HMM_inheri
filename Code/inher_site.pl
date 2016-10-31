#!/usr/bin/perl 
use warnings; use strict;


######################  PARA  ######################
####################################################

my ($jie,$mei,$out)= @ARGV;
open IN,"bedtools intersect -a $jie -b $mei -wo | " or die $!;
open MO,">$out" or die $!;
#open SN,">$out.bed" or die $!;

my %ge2inh = (
	ababaaaa => [0.85,0.05,0.05,0.05], #"i",
	ababaaab => [0.05,0.45,0.45,0.05], #"m:p",
	abababab => [0.45,0.05,0.05,0.45], #"i:n",
	ababaabb => [0.05,0.05,0.05,0.85], #"n",
	aaabaaaa => [0.45,0.45,0.05,0.05], #"i:m",
	aaababab => [0.45,0.45,0.05,0.05], #"i:m",
	aaabaaab => [0.05,0.05,0.45,0.45], #"p:n",
	abaaaaaa => [0.45,0.05,0.45,0.05], #"i:p",
	abaaabab => [0.45,0.05,0.45,0.05], #"i:p",
	abaaaaab => [0.05,0.45,0.05,0.45], #"m:n",
);

my %dec = (
	i=>0,m=>1,p=>2,n=>3
);

my @P = (0.25,0.25,0.25,0.25);   # priori probabilty 
my %psi;      # backpointer matix; 
my %d;         # each probability
my %ip;   # each individual probability
################################################

my %haps;
my %gene_ha;
my @poss = (9,10,11,21,22,23);  # informative positons in bedtools output file
my @rcds = (0,0,0);  # record the last informatio. which will provided for viterbi
while(<IN>){
	chomp;
	my $line = $_;
	my $fb = filter_vcf($line);   # filter low quality snps
	next if ($fb == 0);
	
	my @ar = split "\t",$line;
	
	my $dot = 0;
	my $sla = 0;
	my %tmpha;
	foreach my $i (@poss ){
		if ($ar[$i] =~ /((\S)([\/\|])(\S))/){
			#print "$1\t";
			if ($2 eq "." or $4 eq "."){
				$dot = 1;
			}else{
				my @a = sort {$a <=> $b} ($2,$4);
				$ar[$i] = "$a[0]$a[1]";
				$tmpha{$i}{1} = $2;
				$tmpha{$i}{2} = $4;
				if($3 eq "/"){
					$sla = 1;
				}
			}
		}else{
			die "why!!!$ar[$i]\n";
		}
	}
	if (! $dot and !$sla){
		foreach my $i (@poss){
			$haps{$i}{1} .= $tmpha{$i}{1};
			$haps{$i}{2} .= $tmpha{$i}{2};
		}
	}

	if ($dot){
		#print "\n";
		next;
	}
	my ($chr,$site,$s1,$f1,$m1,$s2,$f2,$m2) = @ar[0,1,9,10,11,21,22,23];
	if ($f1 ne $f2 or $m1 ne $m2){
		#print "\n";
		next;
	}
	
	if ( ($f1+$m1) < 2 or (($f1 + $m1) ==2 and ($s1+$s2) < 2)){
		for my $i ($s1,$s2,$f1,$m1){
			$i =~ tr/01/ab/;
		}
	}elsif( ($f1+$m1) > 2 or (($f1 + $m1) == 2 and  ($s1 + $s2) >= 2)){
		for my $i ($s1,$s2,$f1,$m1){
			$i =~ tr/01/ba/;
			$i = reverse $i;
		}
	}
	my @sibs = sort ($s1,$s2);
	my $g_p = join "",($f1,$m1,@sibs);
	#print "GT:$g_p\t";
	if ($ge2inh{$g_p}){
		my $inh_p = $ge2inh{$g_p};	
		$ip{$chr}{$site} = join ",", @$inh_p;
		#print SN "$chr\t$site\t$site\t@$inh_p\n";
		viterbi(\@rcds,[$chr,$site,$inh_p]);
		@rcds = ($chr,$site); 
	}else{
		#print "\n";
		#print STDERR "NO exist\n";
	}
}

#####################
#### subroutines ####
#####################
sub filter_vcf {
	my $line = shift @_;
	my @ar = split /\t/,$line;
	my $boo = 1;
	my @pedi = @ar[9,21,10,11];
	foreach my $hu (@pedi){
		my ($ge,$va) = (split /:/ , $hu)[0,1];
		#print "WHY $line\n" unless ($va);
		next if ($hu =~ /\./);
		my ($ref,$alt) = split /,/ ,$va;
		if ($ge =~ /1/ and $alt < 3){
			$boo = 0;
			last;
		}
	}

	if ($ar[5] < 100 or $ar[17] < 100){
		$boo = 0;
	}
	return $boo;
}

sub cal_trans {  # defaul average distance : 1cM/Mbs
	my ($dis) = shift @_;
	my $cm = ($dis/1000000) * 0.01 * 4;
	$cm = 0.5 if ($cm > 0.5);
	my $slf = 1-$cm;
	my $cha = $cm/2;
	my $er = 0.005 * $cha;;
	my @a = (
		[$slf,$cha,$cha,$er],
		[$cha,$slf,$er,$cha],
		[$cha,$er,$slf,$cha],
		[$er,$cha,$cha,$slf]
	);
	#print "TT: $slf\t$cha\t$er\n";
	return \@a;
}

sub viterbi {
	my ($pre,$cur) = @_;
	my($pchr,$ppos) = @$pre;
	my($cchr,$cpos,$cpar) = @$cur;
	#print "VIT:$pchr\t$ppos\t$cchr\t$cpos\n";
	if ($pchr ne $cchr){
		$pchr = $cchr;
		$ppos = 0 ;
		$psi{$cchr}{0}= [0,0,0,0];
		$d{$cchr}{0} = \@P;
	}

	my $dis = $cpos - $ppos + 1;
	my $t = cal_trans($dis);
		
	my $low = 0;
	for (my $j = 0; $j < 4; $j ++){
		my @tmparr;
		for (my $i = 0 ; $i < 4; $i ++){
			#print "$d{$cchr}{$ppos}[$i] * $$t[$i][$j] * $$cpar[$j]\n";
			$tmparr[$i] = $d{$cchr}{$ppos}[$i] * $$t[$i][$j] * $$cpar[$j];
		}
		$d{$cchr}{$cpos}[$j] = -1;  # -1 were assigned by chance
		for (my $x = 0; $x < 4; $x ++){
			if ($tmparr[$x] > $d{$cchr}{$cpos}[$j]){
				$d{$cchr}{$cpos}[$j] = $tmparr[$x];
				$psi{$cchr}{$cpos}[$j] = $x;
			}elsif ($tmparr[$x] == $d{$cchr}{$cpos}[$j]){
				#print "Ambiguous path\n";
			}
		}
		if($d{$cchr}{$cpos}[$j] < 1e-200){
			$low = 1;
		}
	}
	if($low == 1){
		#print "NEED to blow up\n";
		for (my $i = 0; $i < 4; $i ++){
			$d{$cchr}{$cpos}[$i] = $d{$cchr}{$cpos}[$i] * 1e200;
		}
	}
=head
	for my $i (@$t){
		print "\t@$i\n";
	}
	
	print "\tEM\t@$cpar\n";
	print "@{$psi{$cchr}{$cpos}}\n";
	print "@{$d{$cchr}{$cpos}}\n";
=cut
}
############## end of subroutines  ##############
#################################################

=head
for my $c (sort keys %d){
	my %dc = %{$d{$c}};
	my %psic = %{$psi{$c}};
	for my $k ( sort {$a<=>$b}  (keys %dc)){
		print "$c\t$k\t@{$dc{$k}}\n";
		print "\t$c\t$k\t@{$psic{$k}}\n";
	}
}
=cut
foreach my $c (sort keys %d){
	my %cinf = %{$d{$c}};
	my @allpos = sort {$b <=> $a} keys %cinf;
	my $last_p = shift @allpos;
    
	my $pstar = $cinf{$last_p}[0];
	my $qstar = 0;
	my $psi_r = $psi{$c}{$last_p}[0];
	for (my $x = 1; $x < 4; $x ++){
		if ( $cinf{$last_p}[$x] > $pstar){
			$pstar = $cinf{$last_p}[$x];
			$qstar = $x;
			$psi_r = $psi{$c}{$last_p}[$x]; 
		}
	} 
	#print MO "$c\t$last_p\t$pstar\t$qstar\n";
	print MO "$c\t$last_p\t$qstar\t$ip{$c}{$last_p}\n";
	
	my $pos_r = $last_p;	
	foreach my $p ( @allpos ){
		my $pro = $d{$c}{$p}[$psi_r];
		#print MO "$c\t$p\t$pro\t$psi_r\n";
		next if ($p == 0);
		print MO "$c\t$p\t$psi_r\t$ip{$c}{$p}\n";
		$pos_r = $p;
		$psi_r = $psi{$c}{$p}[$psi_r];
	}
}


