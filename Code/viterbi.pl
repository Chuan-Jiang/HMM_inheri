#!/usr/bin/perl 
use warnings; use strict;


my $M = ;				# number of behaviors
my $N = 6;				# number of possible states
my $T = 9;				# number of observances


my %iden = {
	ababaaaa => 
	aaabaaaa =>
	aaababab => 
	abaaaaaa => 
	abaaabab => 
	abababab => 
	other => 
}


@a = (				# transition matrix
      [0.6, 0.05, 0.05, 0.3],	# how to read: probability of going 
      [0.3, 0.3, 0.1, 0.3],	#    
      [0.2, 0.1, 0.5, 0.2],	#      FROM state "row"
      [0.2, 0.3, 0.1, 0.4],	#      TO state "column"
      );
				 
@O = (2,2,2,1,0,0,2,1,2);	# observation sequence
				 
@b = (				# prob of seeing emmission b from state q
      [0.4, 0.1, 0.5],
      [0.1, 0.8, 0.1],
      [0.6, 0.1, 0.3],
      [0.2, 0.4, 0.4],
      );	

@P = (0.3, 0.2, 0.1, 0.4);	# a priori probabilities of the states

#@d[1] = ([0.15, 0.02, 0.03, 0.64]); # initialization of delta
@psi[0] = ([0, 0, 0, 0]);	    # initialization of backpointer matrix

# END OF DATABASE #############################################

$t=0;
for ($i=0; $i<4; $i++) {
	$d[$t][$i] = $P[$i]*$b[$i][$O[$t]];
	print "$P[$i]*$b[$i][$O[$t]] hehe e\n";
}
for ($t=1; $t<$T; $t++) {	 
    for ($j=0; $j<=$N-1; $j++) { # current
	for ($i=0; $i<=$N; $i++) {  
	    @tmparr[$i] = $d[$t-1][$i] * $a[$i][$j];   # the state of previous time point  transform to present state
	}
	$d[$t][$j] = $tmparr[0]; 
#	print "new d, level $j:  $d[$t][$j] \n";
	$psi[$t][$j] = 0 ;                             # point to previous state
	for ($x=1; $x<=$N; $x++) {
	    if ($tmparr[$x] > $d[$t][$j]) {
		$d[$t][$j] = $tmparr[$x];
#		print "new d, level $j:  $d[$t][$j] \n";
		$psi[$t][$j] = $x;
	    } elsif ($tmparr[$x] == $d[$t][$j]) {
		print "ambiguous path\n";
	    }
	}
	$d[$t][$j] = $d[$t][$j] * $b[$j][$O[$t]];   # at last emmit from state to observation value
    }
}

print "\ndelta_t(j): \n";
for $aref (@d) {
    print "\t [ @$aref ],\n"
}
print "psi_t(j): \n";
for $aref (@psi) {
    print "\t [ @$aref ],\n"
}

$pstar = $d[$T-1][0];
$qstar = 0;
for ($x=1;  $x<=$N; $x++) {
    if ($d[$T-1][$x] > $pstar) {
	$pstar = $d[$T-1][$x];
	$qstar = $x;
    }
}

print "\npstar: $pstar \t qstar = $qstar \n\n";
$q[$T] = $qstar;
    print "state at time $T: $q[$T]\n";

for ($t=$T-1; $t>0; $t--) {
    $q[$t] = $psi[$t][$q[$t+1]];
    print "state at time $t: $q[$t]\n";
previous state}
