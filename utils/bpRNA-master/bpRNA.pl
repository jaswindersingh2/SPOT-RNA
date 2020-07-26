#!/usr/bin/perl -w
use Graph;
use strict;

#####################################################################
# 2018  Padideh Danaee,Michelle Wiley, Mason Rouches, David Hendrix #
# http://hendrixlab.cgrb.oregonstate.edu                            #
# ----------------------------------------------------------------- #              
# MAIN                                                              #
#####################################################################

my $DEBUG = 0;
my $usage = "Usage:\n$0 <bpseq file or dot-bracket file> \n";

my $inputFile = $ARGV[0] or die $usage;

my($id,$bp,$seq);

if(isBPSeqFile($inputFile)) {
    ($id) = $inputFile =~ /([^\/]*?)\.bpseq/ or die "Could not parse BPSEQ filename: $inputFile\nExpected file extension \".bpseq\"";
    ($bp,$seq) = readBPSeqFile($inputFile);
} elsif(isDotBracketFile($inputFile)) {
    ($id) = $inputFile =~ /([^\/]*?)\.dbn/ or die "Could not parse dotbacket filename: $inputFile\nExpected file extension \".dbn\"";
    ($bp,$seq) = readDotBracketFile($inputFile);    
} else {
    die "Could not determine file type. Expecting BPSEQ or DBN (dot-bracket) file formats.";
}

if(keys %{$bp}) {
    my $allSegments = getSegments($bp);
    for(my $i=0;$i<@{$allSegments};$i++) {
	my $w = @{$allSegments->[$i]};
	print "$i $w\n" if($DEBUG);
    }
    my($segments,$knots,$warnings) = separateSegments($allSegments);
    $bp = filterBasePairs($bp,$knots);
    $segments = getSegments($bp);
    my($dotbracket,$s,$k,$structureTypes,$pageNumber) = buildStructureMap($segments,$knots,$bp,$seq);
    printStructureTypes($id,$seq,$dotbracket,$s,$k,$structureTypes,$pageNumber,$warnings);
} else {
    my $dotbracket = "." x length($seq);
    my $s = "E" x length($seq);
    my $k = "N" x length($seq);
    my $stFile = $id . ".st";
    open(STF,">$stFile") or die "Could not open $stFile for writing.\n";
    print STF "#Name: $id\n";
    print STF "#Length: ", length($seq), "\n";
    print STF "#PageNumber: 0\n";
    print STF "$seq\n$dotbracket\n$s\n$k\n";
    close(STF);
}

###############
# SUBROUTINES #
###############

sub printStructureTypes {
    my($id,$seq,$dotbracket,$s,$k,$structureTypes,$pageNumber,$warnings) = @_;
    $s = join("",@{$s});
    $k = join("",@{$k});
    my $stFile = $id . ".st";
    open(STF,">$stFile") or die "Could not open $stFile for writing.\n";
    print STF "#Name: $id\n";
    print STF "#Length: ", length($seq), "\n";
    print STF "#PageNumber: $pageNumber\n";
    print STF "$warnings"; # warning contains label and newline
    print STF "$seq\n$dotbracket\n$s\n$k\n";
    for my $item (@{$structureTypes->{"S"}}) {
	print STF "$item";
    }    
    for my $item (@{$structureTypes->{"H"}}) {
	print STF "$item";
    }
    for my $item (@{$structureTypes->{"B"}}) {
	print STF "$item";
    }
    for my $i (sort {$a <=> $b} keys %{$structureTypes->{"I"}}) {
	for my $item (@{$structureTypes->{"I"}{$i}}) {
	    print STF "$item";
	}
    }    
    for my $m (sort {$a <=> $b} keys %{$structureTypes->{"M"}}) {
	for my $item (@{$structureTypes->{"M"}{$m}}) {
	    print STF "$item";
	}
    }
    for my $item (@{$structureTypes->{"X"}}) {
	print STF "$item";
    }    
    for my $item (@{$structureTypes->{"E"}}) {
	print STF "$item";
    }    
    for my $item (@{$structureTypes->{"PK"}}) {
	print STF "$item";
    }    
    for my $item (@{$structureTypes->{"PKBP"}}) {
	print STF "$item";
    }    
    for my $item (@{$structureTypes->{"NCBP"}}) {
	print STF "$item";
    }    
    for my $item (@{$structureTypes->{"SEGMENTS"}}) {
	print STF "$item";
    }
    close(STF);
}

sub buildStructureMap {
    my($segments,$knots,$bp,$seq) = @_;
    my($G,$edges) = buildSegmentGraph($seq,$bp,$segments,$knots);
    my($dotbracket,$pageNumber) = computeDotBracket($segments,$knots,$seq);
    my($s,$pk) = computeStructureArray($dotbracket,$bp,$seq);
    my %edges;
    my %pkLoops;
    my %regions;
    my %structureTypes;
    foreach my $edge (@{$edges}) {
	# $u and $v are the segment IDs
	my($u,$v,$s1_pos,$s2_pos,$label) = @{$edge};
	# the start and stop of edge is adjacent to segment:
	my $lStart = $s1_pos+1;
	my $lStop = $s2_pos-1;
	my $type = "M";
	$type = "H" if($u == $v);
	#print "$type\n";
	push(@{$edges{$type}},[$lStart,$lStop]);
	#print "Loaded edge ",$#{$edges{"M"}}," from $u to $v type=$label from $lStart to $lStop\n" if $type eq "M";
    }
    # First, handle Multiloops
    if($edges{"M"}) {
	my $mG = Graph->new; # a directed graph
	for(my $m=0;$m<@{$edges{"M"}};$m++) {
	    # this is important for X regions not adjacent to any other X regions.
	    $mG->add_vertex($m);
	    # 1-based coords
	    my($mStart,$mStop) = @{$edges{"M"}[$m]};
	    my $mSeq = substr($seq,$mStart-1,$mStop-$mStart+1);
	    #print "$m $mStart..$mStop $mSeq\n";
	    for(my $n=0;$n<@{$edges{"M"}};$n++) {
		# 1-based coords
		my($nStart,$nStop) = @{$edges{"M"}[$n]};
		my $nSeq = substr($seq,$nStart-1,$nStop-$nStart+1);
		if(loopLinked($mStart,$mStop,$nStart,$nStop,$bp)) {
		    #print "linking $m to $n\n";
		    $mG->add_edge($m,$n);
		}
	    }
	}
	my $mUG = $mG->undirected_copy();
	my @mCC = $mUG->connected_components();
	my @multiLoops;
	my @externalLoops;
	# create a sorted list of multiloops.
	for my $c (@mCC) {
	    #print "found CC of size ", scalar(@{$c}), "\n";
	    if(isMultiLoop($c,$mG)) {
		# sort the branches of this multiloop by position.
		my @c = sort {$edges{"M"}[$a]->[0] <=> $edges{"M"}[$b]->[0]} @{$c};
		push(@multiLoops,\@c);
	    } else {
		my @c = sort {$edges{"M"}[$a]->[0] <=> $edges{"M"}[$b]->[0]} @{$c};
		push(@externalLoops,\@c);
		#print "found a X region\n";
	    }
	}
	@multiLoops = sort {$edges{"M"}[$a->[0]]->[0] <=> $edges{"M"}[$b->[0]]->[0]} @multiLoops;
	my $mCount = 0;
	for my $c (@multiLoops) {
	    $mCount++;
	    my $mp = 0;
	    for my $m (@{$c}) {
		$mp++;
		my($mStart,$mStop) = @{$edges{"M"}[$m]};
		my $mLength = $mStop - $mStart + 1;
		my $mSeq = substr($seq,$mStart-1,$mLength);
		# define adjacet positions in 1-based
		# retrieve nucleotide in 0-based
		my $bp5_pos1 = $mStart-1;
		my $nuc5_1 = substr($seq,$bp5_pos1-1,1);
		my $bp5_pos2 = $bp->{$bp5_pos1};
		my $nuc5_2 = substr($seq,$bp5_pos2-1,1);
		my $bp3_pos1 = $mStop+1;
		my $nuc3_1 = substr($seq,$bp3_pos1-1,1);
		my $bp3_pos2 = $bp->{$bp3_pos1};
		my $nuc3_2 = substr($seq,$bp3_pos2-1,1);
		my @mKnots = includesKnot($mStart,$mStop,$knots);
		my $PK = @mKnots ? "PK{".join(',',@mKnots)."}" : "";
		push(@{$structureTypes{"M"}{$mCount}},"M$mCount.$mp $mStart..$mStop \"$mSeq\" ($bp5_pos1,$bp5_pos2) $nuc5_1:$nuc5_2 ($bp3_pos1,$bp3_pos2) $nuc3_1:$nuc3_2 $PK\n");
		for my $k (@mKnots) {
		    push(@{$pkLoops{$k}},["M$mCount.$mp",$mStart,$mStop]);
		}
		########################################
		# updating multiloops in %regions here #
		########################################
		push(@{$regions{"M"}},[$mStart,$mStop]);
		# if this inequality doesn't hold, it's a branch of length 0.
		if($mStart <= $mStop) {
		    # in a multiloop part
		    for(my $i=$mStart;$i<=$mStop;$i++) {
			# convert to 0-based
			$s->[$i-1] = "M";
		    }
		}
	    }
	}
	for my $c (@externalLoops) {
	    # these need to be classified as "X" for unpaired
	    my $xCount = 0;
	    my @c = sort {$edges{"M"}[$a]->[0] <=> $edges{"M"}[$b]->[0]} @{$c};
	    for my $x (@c) {
		$xCount++;
		my($xStart,$xStop) = @{$edges{"M"}[$x]};
		my $xSeq = substr($seq,$xStart-1,$xStop-$xStart+1);
		my $bp5_pos1 = $xStart-1;
		my $nuc5_1 = substr($seq,$bp5_pos1-1,1);
		my $bp5_pos2 = $bp->{$bp5_pos1};
		my $nuc5_2 = substr($seq,$bp5_pos2-1,1);
		my $bp3_pos1 = $xStop+1;
		my $nuc3_1 = substr($seq,$bp3_pos1-1,1);
		my $bp3_pos2 = $bp->{$bp3_pos1};
		my $nuc3_2 = substr($seq,$bp3_pos2-1,1);
		#print "checking $u  at $uStart $uStop\n";
		my @xKnots = includesKnot($xStart,$xStop,$knots);
		my $PK = @xKnots ? "PK{".join(',',@xKnots)."}" : "";
		push(@{$structureTypes{"X"}},"X$xCount $xStart..$xStop \"$xSeq\" ($bp5_pos1,$bp5_pos2) $nuc5_1:$nuc5_2 ($bp3_pos1,$bp3_pos2) $nuc3_1:$nuc3_2 $PK\n") if($xStart < $xStop);
		for my $k (@xKnots) {
		    #print "visited $k uKnot";
		    push(@{$pkLoops{$k}},["X$xCount",$xStart,$xStop]);
		}
		# Check if in a multiloop part or something other than U/X
		for(my $i=$xStart;$i<=$xStop;$i++) {
		    # convert to 0-based
		    die "incorrect value at $i $s->[$i-1]\n" unless($s->[$i-1] =~ /[uX]/);
		}
	    }
	}
    }
    # hairpins 
    my $hCount = 0;
    if($edges{"H"}) {
	#print "found H edge\n";
	for(my $h=0;$h<@{$edges{"H"}};$h++) {
	    $hCount++;
	    my($hStart,$hStop) = @{$edges{"H"}[$h]};
	    # -1 for zero-based
	    my $hSeq = substr($seq,$hStart-1,$hStop-$hStart+1);	    
	    # -1 again for flanking
	    my $pos5 = $hStart-1;
	    #  +1 for flanking
	    my $pos3 = $hStop+1;
	    # subtrackt 1 to get zero-based to get nucleotides flanking (closing base pair)
	    my $nuc1 = substr($seq,$pos5-1,1);
	    my $nuc2 = substr($seq,$pos3-1,1);
	    push(@{$regions{"H"}},[$hStart,$hStop]);
	    my @hKnots = includesKnot($hStart,$hStop,$knots);
	    my $PK = @hKnots ? "PK{".join(',',@hKnots)."}" : "";
	    push(@{$structureTypes{"H"}},"H$hCount $hStart..$hStop \"$hSeq\" ($pos5,$pos3) $nuc1:$nuc2 $PK\n");
	    #print "stored H into structureTypes\n";
	    for my $k (@hKnots) {
		#print "visited $k hknot\n";
		push(@{$pkLoops{$k}},["H$hCount",$hStart,$hStop]);
	    }
	}
    }
    ###################
    # Extract regions #
    ###################
    my $prevChar = "";
    my $thisStart = 0;
    my $thisStop = 0;
    my $FIRST = 1;
    ####################################################
    # Stems, internal loops, multiloops, hairpin loops #
    ####################################################
    for(my $i=0;$i<@{$s};$i++) {
	unless($s->[$i]) {
	    die "no value in \@s for $i\n"
	}
	if($s->[$i] ne $prevChar) {
	    # in a new region.
	    # store previous region info:	    
	    # ignore multiloops and hairpins, they are already added above
	    unless($prevChar =~ /[MH]/) {
		# store in 1-based to be consistent with M and H above
		push(@{$regions{$prevChar}},[$thisStart+1,$thisStop+1]) unless($FIRST);
	    }
	    # reset the start/stop:
	    $thisStart = $i;
	    $thisStop = $i;
	    $prevChar = $s->[$i];
	    $FIRST = 0;
	} elsif(($prevChar eq "S")&&($bp->{$i+1}+1 != $bp->{$i})) {
	    # in a break in a stem into another stem. store in 1-based
	    push(@{$regions{$prevChar}},[$thisStart+1,$thisStop+1]);
	    $thisStart = $i;
            $thisStop = $i;
            $prevChar = $s->[$i];
	    $FIRST = 0;
	} else {
	    # continue/extend the same region
	    $thisStop = $i;
	}
    }
    # store the last region.
    push(@{$regions{$prevChar}},[$thisStart+1,$thisStop+1]) unless($FIRST);    
    #bulges
    if($regions{"B"}) {
	my $bCount = 0;
	foreach my $bulge (@{$regions{"B"}}) {
	    # 1-based
	    my($bStart,$bStop) = @{$bulge};
	    # all substr calls need 0-based, hence -1 for all:
	    my $bSeq = substr($seq,$bStart-1,$bStop-$bStart+1);	    
	    # 1-based positions of flanking:
	    my $bp5_pos1 = $bStart-1;
	    my $bp3_pos1 = $bStop+1;
	    my $bp5_nt1 = substr($seq,$bp5_pos1-1,1);
	    my $bp3_nt1 = substr($seq,$bp3_pos1-1,1);
	    my $bp5_pos2 = $bp->{$bp5_pos1};
	    my $bp3_pos2 = $bp->{$bp3_pos1};
	    my $bp5_nt2 = substr($seq,$bp5_pos2-1,1);
	    my $bp3_nt2 = substr($seq,$bp3_pos2-1,1);
	    $bCount++;
	    my @bKnots = includesKnot($bStart,$bStop,$knots);
	    my $PK = @bKnots ? "PK{".join(',',@bKnots)."}" : "";
	    push(@{$structureTypes{"B"}},"B$bCount $bStart..$bStop \"$bSeq\" ($bp5_pos1,$bp5_pos2) $bp5_nt1:$bp5_nt2 ($bp3_pos1,$bp3_pos2) $bp3_nt1:$bp3_nt2 $PK\n");
	    for my $k (@bKnots) {
		push(@{$pkLoops{$k}},["B$bCount",$bStart,$bStop]);
	    }
	}
    }
    if($regions{"I"}) {
	my $iG = Graph::Undirected->new;
	for(my $i=0;$i<@{$regions{"I"}};$i++) {
	    my($iStart,$iStop) = @{$regions{"I"}[$i]};
	    for(my $j=$i+1;$j<@{$regions{"I"}};$j++) {
		my($jStart,$jStop) = @{$regions{"I"}[$j]};
		if(loopLinked($iStart,$iStop,$jStart,$jStop,$bp)) {
		    $iG->add_edge($i,$j);
		}
	    }
	}
	my @iCC = $iG->connected_components();
	my @internalLoops;
	for my $c (@iCC) {
	    my @c = sort {$regions{"I"}[$a]->[0] <=> $regions{"I"}[$b]->[0]} @{$c};
	    push(@internalLoops,\@c);
	}
	my $iCount = 0;
	@internalLoops = sort {$regions{"I"}[$a->[0]]->[0] <=> $regions{"I"}[$b->[0]]->[0]} @internalLoops;
	for my $c (@internalLoops) {
	    $iCount++;
	    my $componentSize = scalar(@{$c});
	    my $ip = 0;
	    for my $v (@{$c}) {
		$ip++;		
		# 1-based positions
		my($iStart,$iStop) = @{$regions{"I"}[$v]};
		my $iLength = $iStop - $iStart + 1;
		# subtract 1 to convert to 0-based
		my $iSeq = substr($seq,$iStart-1,$iLength);
		# positions in 1-based here:
		my $bp5_pos1 = $iStart-1;
		my $bp5_pos2 = $bp->{$bp5_pos1};
		my $nuc5_1 = substr($seq,$bp5_pos1-1,1);
		my $nuc5_2 = substr($seq,$bp5_pos2-1,1);
		my @iKnots = includesKnot($iStart,$iStop,$knots);
		my $PK = @iKnots ? "PK{".join(',',@iKnots)."}" : "";
		push(@{$structureTypes{"I"}{$iCount}},"I$iCount.$ip $iStart..$iStop \"$iSeq\" ($bp5_pos1,$bp5_pos2) $nuc5_1:$nuc5_2 $PK\n");
		for my $k (@iKnots) {
		    push(@{$pkLoops{$k}},["I$iCount.$ip",$iStart,$iStop]);
		}
	    }
	}
    }
    if($regions{"E"}) {
	my $eCount = 0;
	for(my $e=0;$e<@{$regions{"E"}};$e++) {
	    my($eStart,$eStop) = @{$regions{"E"}[$e]};
	    my $eSeq = substr($seq,$eStart-1,$eStop-$eStart+1);
	    $eCount++;
	    my @eKnots = includesKnot($eStart,$eStop,$knots);
	    my $PK = @eKnots ? "PK{".join(',',@eKnots)."}" : "";
	    push(@{$structureTypes{"E"}},"E$eCount $eStart..$eStop \"$eSeq\" $PK\n");
	    for my $k (@eKnots) {
		push(@{$pkLoops{$k}},["E$eCount",$eStart,$eStop]);
	    }
	}
    }
    #stems
    my @stemList;
    my %VISITED; # a hash to keep track of regions in stems that have been collected already
    if($regions{"S"}) {
	my $sCount = 0;
	foreach my $stem (sort {$a->[0] <=> $b->[0]} @{$regions{"S"}}) {
	    my($sStart1,$sStop1) = @{$stem};
	    unless($VISITED{$sStart1}{$sStop1}) {
		$sCount++;
		# 1-based
		my $sStart2 = $bp->{$sStop1};
		my $sStop2 = $bp->{$sStart1};
		# all substr calls need 0-based, hence -1 for all:
		my $sSeq1 = substr($seq,$sStart1-1,$sStop1-$sStart1+1);
		my $sSeq2 = substr($seq,$sStart2-1,$sStop2-$sStart2+1);
		push(@{$structureTypes{"S"}},"S$sCount $sStart1..$sStop1 \"$sSeq1\" $sStart2..$sStop2 \"$sSeq2\"\n");
		$VISITED{$sStart2}{$sStop2} = 1;
		push(@stemList,[$sStart1,$sStop1,"S$sCount"]);
	    }
	}
    }
    my $ncCount=0;
    for my $i (keys %{$bp}) {
	my $j = $bp->{$i};
	if($i < $j) {
	    my $b1 = substr($seq,$i-1,1);
	    my $b2 = substr($seq,$j-1,1);
	    if(nonCanonical($b1,$b2)) {	    
		$ncCount++;
		my $thisLabel = "";
		foreach my $stem (@stemList) {
		    my($start,$stop,$label) = @{$stem};
		    if(($start<=$i)&&($i<=$stop)) {
			$thisLabel = $label;
		    }
		}
		die "No label found for $i,$j\n" unless($thisLabel);
		push(@{$structureTypes{"NCBP"}},"NCBP$ncCount $i $b1 $j $b2 $thisLabel\n");
	    }
	}
    }
    if(@{$knots}) {
	for(my $k=0;$k<@{$knots};$k++) {
	    my $knotID = $k+1;
	    my @knot = @{$knots->[$k]};
	    my $knotSize = @knot;
	    my $first = shift(@knot);
	    my($k_5pStart,$k_3pStart) = @{$first};
	    my $last = @knot ? pop(@knot) : $first;
	    my($k_3pStop,$k_5pStop) = @{$last};	    
	    my $linkedLoops = "";
	    if(@{$pkLoops{$knotID}} == 2) {		
		@{$pkLoops{$knotID}} = sort {$a->[1] <=> $b->[1]} @{$pkLoops{$knotID}};
		my($lType1,$lStart1,$lStop1) = @{$pkLoops{$knotID}[0]};
		my($lType2,$lStart2,$lStop2) = @{$pkLoops{$knotID}[1]};
		$linkedLoops = "$lType1 $lStart1..$lStop1 $lType2 $lStart2..$lStop2";
	    } else {
		die "Expected two loops linked for PK$knotID\n";
	    }
	    push(@{$structureTypes{"PK"}},"PK$knotID ${knotSize}bp $k_5pStart..$k_3pStop $k_5pStop..$k_3pStart $linkedLoops\n");
	    push(@stemList,[$k_5pStart,$k_3pStop,"PK$knotID"]);
	    my $n=0;
	    for my $pair (@{$knots->[$k]}) {
		my($k_5p,$k_3p) = @{$pair};
		#print "knot $k: $k_5p $k_3p\n";
		# positions in 1-based 
		my $b_5p = substr($seq,$k_5p-1,1);
		my $b_3p = substr($seq,$k_3p-1,1);
		$n++;
		push(@{$structureTypes{"PKBP"}},"PK$knotID.$n $k_5p $b_5p $k_3p $b_3p\n");
		# check if PKBP is non-canonical.
		if(nonCanonical($b_5p,$b_3p)) {
		    $ncCount++;
		    push(@{$structureTypes{"NCBP"}},"NCBP$ncCount $k_5p $b_5p $k_3p $b_3p PK$knotID.$n\n");
		}
	    }
	}
    }    
    for(my $i=0;$i<@{$segments};$i++) {
	my $segmentID = $i+1;
	my @segment = @{$segments->[$i]};
	my $segSize = @segment;
	my $firstPair = shift(@segment);	    
	my($seg5pStart, $seg3pStart)= @{$firstPair};
	my $lastPair =  @segment ? pop(@segment) : $firstPair;
	my($seg3pStop, $seg5pStop) = @{$lastPair};
	my $segSeq1 = substr($seq,$seg5pStart-1,$seg3pStop-$seg5pStart+1);
	my $segSeq2 = substr($seq,$seg5pStop-1,$seg3pStart-$seg5pStop+1);
 	push(@{$structureTypes{"SEGMENTS"}},"segment$segmentID ${segSize}bp $seg5pStart..$seg3pStop $segSeq1 $seg5pStop..$seg3pStart $segSeq2\n");
    }
    my @k = split(//,"N" x length($dotbracket));
    for(my $i=0;$i<@{$s};$i++) {
	$k[$i] = "K" if($pk->[$i]);
    }
    return($dotbracket,$s,\@k,\%structureTypes,$pageNumber);
}

sub buildSegmentGraph {
    my($seq,$bp,$segments,$knots) = @_;
    # 
    # 5'start   3'stop
    #       ACGUA
    # start ||||| stop
    #       UGCAU
    # 3'start   5'stop
    #
    my($firstPos,$lastPos) = getExtremePositions($bp);
    my $G = Graph->new(multiedged => 1);
    my @edges;
    for(my $i=0;$i<@{$knots};$i++) {
	my @knot1 = @{$knots->[$i]};
	my $firstPair1 = shift(@knot1);	    
	my($k1_5pStart, $k1_3pStart)= @{$firstPair1};
	my $lastPair1 =  @knot1 ? pop(@knot1) : $firstPair1;
	my($k1_3pStop, $k1_5pStop) = @{$lastPair1};
	my $k1Seq1 = substr($seq,$k1_5pStart-1,$k1_3pStop-$k1_5pStart+1);
	my $k1Seq2 = substr($seq,$k1_5pStop-1,$k1_3pStart-$k1_5pStop+1);
    }
    if(@{$segments} > 0) {
	for(my $i=0;$i<@{$segments};$i++) {
	    $G->add_vertex($i);
	    my @segment1 = @{$segments->[$i]};
	    my $firstPair1 = shift(@segment1);
	    my($s1_5pStart, $s1_3pStart)= @{$firstPair1};
	    my $lastPair1 =  @segment1 ? pop(@segment1) : $firstPair1;
	    my($s1_3pStop, $s1_5pStop) = @{$lastPair1};
	    my $s1Seq1 = substr($seq,$s1_5pStart-1,$s1_3pStop-$s1_5pStart+1);
	    my $s1Seq2 = substr($seq,$s1_5pStop-1,$s1_3pStart-$s1_5pStop+1);
	    # get the next paired bases in each direction (see diagram above)
	    my $p1_5pStart = getPrevPair($s1_5pStart,$bp,$firstPos,$knots); # 1
	    my $n1_3pStop = getNextPair($s1_3pStop,$bp,$lastPos,$knots);    # 2
	    my $p1_5pStop = getPrevPair($s1_5pStop,$bp,$firstPos,$knots);   # 3
	    my $n1_3pStart = getNextPair($s1_3pStart,$bp,$lastPos,$knots);  # 4	    
	    for(my $j=0;$j<@{$segments};$j++) {
		my @segment2 = @{$segments->[$j]};
		my $firstPair2 = shift(@segment2);
		my($s2_5pStart, $s2_3pStart)= @{$firstPair2};
		my $lastPair2 =  @segment2 ? pop(@segment2) : $firstPair2;
		my($s2_3pStop, $s2_5pStop) = @{$lastPair2};
		if($n1_3pStart == $s2_5pStart) {
		    $G->add_edge($i,$j);
		    # extra label "1" is for debugging.
		    push(@edges,[$i,$j,$s1_3pStart,$s2_5pStart,"1"]);
		    #print "Connected $i to $j\n";
		}
		if($n1_3pStart == $s2_5pStop) {
		    $G->add_edge($i,$j);
		    push(@edges,[$i,$j,$s1_3pStart,$s2_5pStop,"2"]);
		    #print "Connected $i to $j\n";
		}
		if($n1_3pStop == $s2_5pStart) {
		    $G->add_edge($i,$j);
		    push(@edges,[$i,$j,$s1_3pStop,$s2_5pStart,"3"]);
		    #print "Connected $i to $j\n";
		}
		if($n1_3pStop == $s2_5pStop) {
		    $G->add_edge($i,$j);
		    push(@edges,[$i,$j,$s1_3pStop,$s2_5pStop,"4"]);
		    #print "Connected $i to $j\n";
		}
	    }
	}
    }
    return($G,\@edges);
}

sub computeStructureArray {
    my($dotbracket,$bp,$seq) = @_;
    # Takes bpseq file and adds column listing structure type of each base:
    # paired "Stem"     S
    # Multiloop 	M
    # Internal loop     I
    # Bulge 	        B
    # Hairpin loop      H
    # pseudoKnot	K
    # dangling End	E
    # eXternal loop     X
    my $knotBracket = getKnotBrackets();
    my @s;
    for(my $i = 0; $i < length($dotbracket); $i++) {
	my $x = substr($dotbracket, $i, 1); 
	my $loopStructure = "";
	# stem
	if($x eq "(" || $x eq ")") {
	    $loopStructure = "S";
	} elsif(($x eq ".")||($knotBracket->{$x})) {
	    # loops
	    my ($fwd,$fwdIndex) = fwdFinder($i,$dotbracket,$x,$knotBracket);
	    my ($bwd,$bwdIndex) = bwdFinder($i,$dotbracket,$x,$knotBracket);
	    my $fwdIndexPair = $bp->{$fwdIndex+1}; # returns position on RNA (1-based)
	    my $bwdIndexPair = $bp->{$bwdIndex+1};
	    my $length = $fwdIndex-$bwdIndex-1;
	    if($bwd eq "("){
		if ($fwd eq "("){
		    if ($fwdIndexPair == $bwdIndexPair-1) {
			$loopStructure = "B";
		    } else {
			$loopStructure = between($fwdIndexPair-1,$bwdIndexPair-1,$dotbracket);
		    }
		} elsif($fwd eq ")") {
		    $loopStructure = "H";
		} elsif($fwd eq "") {
		    $loopStructure = "E";
		} else {
		    $loopStructure = "What is this fwd1: $fwd";
		}
	    } elsif ($bwd eq ")") {
		if ($fwd eq "(") {
		    $loopStructure = "X";
		} elsif ($fwd eq ")") {
		    if ($fwdIndexPair == $bwdIndexPair-1) {
			$loopStructure = "B";
		    } else {
			$loopStructure = between($fwdIndexPair-1,$bwdIndexPair-1,$dotbracket);
		    }
		} elsif ($fwd eq "") {
		    $loopStructure = "E";
		} else {
		    $loopStructure = "What is this fwd2: $fwd";
		}
	    } elsif ($bwd eq "" || $bwd eq "."){
		$loopStructure = "E";
	    } else {
		$loopStructure = "What is this bwd1: $bwd";
	    }
	}
	$s[$i] = $loopStructure;
    }
    my @pk = split(//,"0" x length($dotbracket));
    for(my $i = 0; $i < length($dotbracket); $i++) {
	my $x = substr($dotbracket, $i, 1); 
	if($knotBracket->{$x}) {
	    $pk[$i] = 1;
	}
    }
    return(\@s,\@pk);
}

sub printStructureData {
    my($regions) = @_;
    # collect lines for output file
    my @list;
    foreach my $type (keys %{$regions}) {
	if(($type)&&($type ne "N")) {
	    my $count = 0;
	    foreach my $region (@{$regions->{$type}}) {
		$count++;
		my($start,$stop) = @{$region};
		push(@list,[$start,"$id\tbpRNA\t$type\t$start\t$stop\t.\t+\t.\tID=$type$count"]);
	    }
	}
    }
    my $gffFile = $id."_structure.gff";
    open(GFF,">$gffFile") or die "Could not open $gffFile for writing.\n";
    foreach my $region (sort {$a->[0] <=> $b->[0]} @list) {
	my($start,$line) = @{$region};
	print GFF "$line\n";
    }
    close(GFF);
}

sub isMultiLoop() {
    my($c,$mG) = @_;    
    $,=" ";
    # make copy of components
    my @c = @{$c};
    if(scalar(@c) == 1) {
	return 0;
    }
    my $f = shift(@c);
    my $v = $f;
    while(@c) {
	if(my @successors = $mG->successors($v)) {
	    # it should only have one successor. Reality check here.
	    die "Fatal error: Too many successors of multiloop part.\n$v is connectd to @successors\n" unless(@successors == 1);
	    my $w = shift(@successors);
	    my $index = -1;	
	    for(my $i=0;$i<@c;$i++) {
		if($c[$i] == $w) {
		    $index = $i;
		}
	    }
	    if($index >= 0) {
		splice(@c, $index, 1);
	    } else {
		return 0;
	    }
	    $v = $w;	
	} else {
	    return 0;
	}
    }
    # if here, all vertices in @c were removed.
    # check if successor of last node is first node. 
    return 1;
}

sub computeDotBracket {
    my($segments,$knots,$seq) = @_;
    my $dotbracket = "." x length($seq);
    ## For the segments, set all pairs to parentheses
    for(my $i=0;$i<@{$segments};$i++) {    
	foreach my $pair (@{$segments->[$i]}) {
	    my($l,$r) = @{$pair};
	    substr($dotbracket,$l-1,1) = "(";
	    substr($dotbracket,$r-1,1) = ")";	    
	}
    }
    # loop through each base pair, update bracket
    my %page;    
    my $n = 1;
    my $unlabeledKnots = 1;
    my $firstI = 0;
    if(0) {
	for(my $i=0;$i<@{$knots};$i++) {
	    #print "knot $i has ", scalar(@{$knots->[$i]}), "bp\t";
	    my @knot1 = @{$knots->[$i]};
	    my $first1 = shift(@knot1);
	    my($k1_5pStart,$k1_3pStart) = @{$first1};
	    my $last1 = @knot1 ? pop(@knot1) : $first1;
	    my($k1_3pStop,$k1_5pStop) = @{$last1};
	    #print "$k1_5pStart..$k1_3pStart\n";
	}
    }
    while($unlabeledKnots) {
	# assume no more left
	$unlabeledKnots = 0;
	for(my $i=0;$i<@{$knots};$i++) {	
	    unless($page{$i}) {
		$firstI = $i;
		last;
	    }	
	}
	$page{$firstI} = $n;
	# the most 5' knot with undefined page is what everything is compared to 
	# start an array contaning members of this page
	my @checkList = ();
	push(@checkList,$firstI);
	# start at next knot with undefined page
	for(my $i=$firstI+1;$i<@{$knots};$i++) {
	    unless($page{$i}) {
		# for now consider it not labeled
		my $CROSSING  = 0;
		# if it doesn't cross, then use $n, otherwise label later.
		for my $checkI (@checkList) {
		    if(knotsCross($knots->[$checkI],$knots->[$i])) {
			# there is an unlabeled knot left.
			$unlabeledKnots = 1;
			$CROSSING = 1;
		    }
		}
		unless($CROSSING) {
		    $page{$i} = $n;
		    push(@checkList,$i);
		}
	    }
	}
	$n++;
    }
    my $pageNumber = 0;
    for(my $i=0;$i<@{$knots};$i++) {
	my $n = $page{$i}; # page of PK. For all PKs, n > 0, where n=0 is "("
	if($n > $pageNumber) {
	    $pageNumber = $n;
	}
	foreach my $pair (@{$knots->[$i]}) {
	    my($l,$r) = @{$pair};
	    my($lB,$rB) = getBrackets($n);
	    substr($dotbracket,$l-1,1) = $lB;
	    substr($dotbracket,$r-1,1) = $rB;
	}	
    }
    $pageNumber += 1; # convert to 1-based
    return($dotbracket,$pageNumber);
}

sub knotsOverlap {
    my($knot1,$knot2) = @_;
    # treat $knot like a segment
    #
    # 5'start   3'stop
    #       ACGUA
    #       |||||
    #       UGCAU
    # 3'start   5'stop
    #	
    my @knot1 = @{$knot1};
    my $first1 = shift(@knot1);
    my($k1_5pStart,$k1_3pStart) = @{$first1};
    my $last1 = @knot1 ? pop(@knot1) : $first1;
    my($k1_3pStop,$k1_5pStop) = @{$last1};
    my @knot2 = @{$knot2};
    my $first2 = shift(@knot2);
    my($k2_5pStart,$k2_3pStart) = @{$first2};
    my $last2 = @knot2 ? pop(@knot2) : $first2;
    my($k2_3pStop,$k2_5pStop) = @{$last2};
    if(($k1_5pStart <= $k2_5pStart)&&($k2_5pStart <= $k1_3pStart)||
       ($k1_5pStart <= $k2_3pStart)&&($k2_3pStart <= $k1_3pStart)) {
	return 1;
    } 
    return 0;
}

sub knotsCross {
    my($knot1,$knot2) = @_;
    # treat $knot like a segment
    #
    # 5'start   3'stop
    #       ACGUA
    #       |||||
    #       UGCAU
    # 3'start   5'stop
    #	
    my @knot1 = @{$knot1};
    my $first1 = shift(@knot1);
    my($k1_5pStart,$k1_3pStart) = @{$first1};
    my @knot2 = @{$knot2};
    my $first2 = shift(@knot2);
    my($k2_5pStart,$k2_3pStart) = @{$first2};
    if(pkQuartet($k1_5pStart,$k1_3pStart,$k2_5pStart,$k2_3pStart)) {
	return 1;
    }
    return 0;
}

sub loopLinked {
    # 0123456789X
    # ((.))((.))...
    my($iStart1,$iStop1,$iStart2,$iStop2,$bp) = @_;
    # here it is forced to be in order 5' to 3'
    if($iStop1+1 == $bp->{$iStart2-1}) {
	return 1;
    }
    return 0;
}

sub getKnotBrackets {
    my $knots = "[]{}<>" . join("",("a".."z")) . join("",("A".."Z"));
    my %knotBracket;
    foreach my $c (split(//,$knots)) {
	$knotBracket{$c} = 1;
    }
    return \%knotBracket;
}

sub getBrackets {
    my($n) = @_;
    die "Fatal error: too many (n>29) PKs to represent in dotbracket! $id\n" if($n >= 30);
    my $left = "([{<" . join("",("A".."Z"));
    my $right = ")]}>" . join("",("a".."z"));
    return (substr($left,$n,1),substr($right,$n,1));
}

# Find index of the next paired base
sub fwdFinder {
    my ($i,$dotbracket,$x,$knotBracket) = @_;
    my $B = substr($dotbracket,$i,1);
    while (($B eq ".")||($knotBracket->{$B})){
	$i++;
	$B = substr($dotbracket,$i,1);
	if ($i >= length($dotbracket)) {
	    return("",$i);
	}
    } 
    return($B,$i);
}

# Find index of the previous paired base
sub bwdFinder {
    my ($i,$dotbracket,$x,$knotBracket) = @_; 
    my $B = substr($dotbracket,$i,1);
    while (($B eq ".")||($knotBracket->{$B})){
	$i--;
	$B = substr($dotbracket,$i,1);
	if ($i < 0){
	    return("",$i);
	}
    }
    return($B,$i);
}

# Determine if multiloop or internal loop
# Usage:			# get to zero base
# $type = between($fwdPairIndex-1,$bwdPairIndex-1,$dotbracket);
sub between {
    my ($fwdIP,$bwdIP,$dotbracket) = @_;
    $fwdIP++;
    while ($fwdIP < $bwdIP-1){
	my $fwdP = substr($dotbracket,$fwdIP,1);
	if ($fwdP eq "(" || $fwdP eq ")"){
	    return("X");
	}
	$fwdIP++; 
    }
    return("I");
}

sub pkQuartet {
    my($i,$j,$k,$l) = @_;
    # assumption: i < j  and k < l
    if((($i < $k)&&($k < $j)&&($j < $l))||
       (($i < $l)&&($l < $j)&&($k < $i))) {
        return 1;
    }
    return 0;
}

sub getSegments {
    my($bp) = @_;
    #
    # 5'start   3'stop
    #       ACGUA
    #       |||||
    #       UGCAU
    # 3'start   5'stop
    #
    my @allSegments = ();
    my($firstPos,$lastPos) = getExtremePositions($bp);
    # initialize the i,j
    my $i = getNextPair(0,$bp,$lastPos);
    if($i) {
	my $j = $bp->{$i};
	# Build the first segment starting with the first basepair
	# while still examining paired positions		
	while(($firstPos <= $i) && ($i <= $lastPos)) {
	    my $INSEGMENT;
	    my @segment = ();
	    if($i<$j){
		push(@segment,[$i,$j]);
		$INSEGMENT=1;
	    }
	    # grow the segment to include more base pairs
	    while($INSEGMENT) {
		#print "getting nextPair for $i and prevPair of $j\n";
		my $n = getNextPair($i,$bp,$lastPos);
		my $p = getPrevPair($j,$bp,$firstPos);
		if($n && $p) {
		    if($bp->{$n} == $p && $n<$p) {
			# in the segment. Store these basepairs
			push(@segment,[$n,$p]);
			$i = $n;
			$j = $p;	    
		    } else {		
			$INSEGMENT = 0;
			push(@allSegments,\@segment);		
		    }
		} else {
		    $INSEGMENT = 0;
		    push(@allSegments,\@segment);
		}
	    }
	    # go to first base pair of NEXT segment
	    $i = getNextPair($i,$bp,$lastPos);
	    $j = $bp->{$i};
	}
    }
    # by construction, the segments should be ordered 
    # by the most 5' position of each segment
    return \@allSegments;
}

sub filterBasePairs {
    my($bp,$knots) = @_;
    my %bp;
    for my $i ( keys %{$bp} ) {
	my $j = $bp->{$i};
	unless(inKnot($i,$knots)) {
	    unless(inKnot($j,$knots)) {
		$bp{$i} = $j;
	    } else {
		die "filterBasePairs: $j in PK, but $i is not.\n";
	    }
	}
    }
    return \%bp;
}

sub separateSegments {
    my($segments) = @_;
    my $warnings = "";
    my %knot;
    my $G = Graph::Undirected->new;    
    if(@{$segments} > 1) {
	for(my $i=0;$i<@{$segments}-1;$i++) {
	    my @segment1 = @{$segments->[$i]};
	    my $firstPair1 = shift(@segment1);
	    my($s1_5pStart, $s1_3pStart)= @{$firstPair1};
	    #print "$i ($s1_5pStart, $s1_3pStart)\n";
	    for(my $j=$i+1;$j<@{$segments};$j++) {
		my @segment2 = @{$segments->[$j]};
		my $firstPair2 = shift(@segment2);
		my($s2_5pStart, $s2_3pStart)= @{$firstPair2};
		if(pkQuartet($s1_5pStart,$s1_3pStart,$s2_5pStart,$s2_3pStart)) {
		    $G->add_edge($i,$j);
		}
	    }
	}
    }
    my @CC = $G->connected_components();
    print $G, "\n" if($DEBUG);
    my $ccCount=0;
    for my $c (@CC) {
	print "Primary Connected Component: ", $G->subgraph($c), "\t" if($DEBUG);
	for my $v (@{$c}) {
	    print "w($v)=",scalar(@{$segments->[$v]})," " if($DEBUG);
	}
	print "\n" if($DEBUG);
	$ccCount++;
	my($knotsList,$warning) = getBestKnots($G,$c,$segments);
	for my $v (@{$knotsList}) {
	    #print "1. removing $v\n" if($DEBUG);
	    $G->delete_vertex($v);
	    $knot{$v}++;
	}
	$warnings .= $warning;
    }
    if(@{$segments} == keys(%knot)) {
	# everything is involved in a PK! 
	# assign the largest to not be a PK segment
	my $maxSize = 0;
	my $maxI;
	for(my $i=0;$i<@{$segments};$i++) {
	    if(@{$segments->[$i]} > $maxSize) {
		$maxI = $i;
		$maxSize = @{$segments->[$i]}
	    }
	}
	# remove largest segment from knot list:
	delete($knot{$maxI});
    }
    my @segments;
    my @knots;    
    for(my $i=0;$i<@{$segments};$i++) {
	if($knot{$i}) {
	    push(@knots,$segments->[$i]);
	} else {
	    push(@segments,$segments->[$i]);
	}
    }
    return(\@segments,\@knots,$warnings);
}

sub getMinVPair {
    my($c,$segments) = @_;
    my $v = $c->[0];
    my $w = $c->[1];
    # check the pk sequence
    if(@{$segments->[$v]} < @{$segments->[$w]}) {
	return ($v,"");
    } elsif(@{$segments->[$w]} < @{$segments->[$v]}) {
	return ($w,"");
    } else {
	my $Nv = @{$segments->[$v]};
	my $Nw = @{$segments->[$w]};
	print "checking $v and $w (CC of size 2), w($v)=$Nv and w($w)=$Nw\n" if($DEBUG);
	print "$v: " if($DEBUG);
	foreach my $pair (@{$segments->[$v]}) {
	    my($f,$t)=@{$pair};
	    print "$v: $f,$t\n" if($DEBUG);
	}
	print "$w: " if($DEBUG);
	foreach my $pair (@{$segments->[$w]}) {
	    my($f,$t)=@{$pair};
	    print "$w: $f,$t\n" if($DEBUG);
	}
	my $minV = min($v,$w);
	print "selected $minV\n" if($DEBUG);
	my $warning .= "#Warning: Structure contains linked PK-segments of same sizes $Nv and $Nw. Using PK brackets for the more 5' segment\n";	
	return($minV,$warning);
    }
}

sub getBestKnots {
    my($G,$c,$segments) = @_;
    my $warnings = "";
    # we have a "complex" connected component, and need heuristics
    my @knotsList;
    my $g = $G->subgraph($c);
    my $knotsRemain = 1;
    while($knotsRemain) {	    
	print "getBestKnots: Current graph: ", $g, "\n" if($DEBUG);
	$knotsRemain = 0;
	# recompute connected components, updated each iteration
	my @CC = $g->connected_components();
	for my $cc (@CC) {
	    $,=",";
	    print "Checking component: ", @{$cc}, "\n" if((@{$cc} > 1)&&($DEBUG));
	    my $componentSize = scalar(@{$cc});
	    if($componentSize == 2) {
		$knotsRemain = 1;
		print "2-checking: ", @{$cc}, "\n" if($DEBUG);
		my($minV,$warning) = getMinVPair($cc,$segments);
		$warnings .= $warning;
		$g->delete_vertex($minV);
		push(@knotsList,$minV);
		print "2-deleted $minV\n" if($DEBUG);
	    } elsif($componentSize > 2) {
		if(my $path = isPathGraph($g,$cc)) {
		    print "found a path: " if($DEBUG);
		    my @path = @{$path};
		    for(my $i=0;$i<@path;$i++) {
			my $v=$i+1;
			print "[ $v:$path[$i] " if($DEBUG);
			print "w($path[$i])=",scalar(@{$segments->[$path[$i]]})," ] " if($DEBUG);
		    }
		    print "\n" if($DEBUG);
		    my $w1 = scalar(@{$segments->[$path[0]]});
		    my $w2 = scalar(@{$segments->[$path[1]]});
		    my $w3 = scalar(@{$segments->[$path[2]]});		    
		    if((@path == 3)&&($w2 == $w1+$w3)) {
			my $v = $path[1];
			$g->delete_vertex($v);
			push(@knotsList,$v);
			print "3-path-deleted $v\n" if($DEBUG);
			$knotsRemain = 1;
		    } else {
			my @maxSet;
			$maxSet[0] = 0;
			# weight of first term in path:
			$maxSet[1] = scalar(@{$segments->[$path[0]]});
			for(my $i=2;$i<@path+1;$i++) {
			    # the @path index is 0-based, but maxSet is 1-based
			    # 0th term is 0.
			    $maxSet[$i] = max($maxSet[$i-1],$maxSet[$i-2] + scalar(@{$segments->[$path[$i-1]]}));
			    #print "maxSet $i $maxSet[$i]\n" if($DEBUG);
			}
			print @maxSet, "length=",scalar(@maxSet),"\n" if($DEBUG);
			my $i = @path;
			my %maxWeighted;
			# only apply path algorithm to all but last 2 members
			while($i >= 1) {		
			    my $weight1 = scalar(@{$segments->[$path[$i-2]]});
			    my $weight2 = scalar(@{$segments->[$path[$i-1]]});
			    print "checking $i: $maxSet[$i-2] + $weight2 vs $maxSet[$i-1]\n" if($DEBUG);
			    if($i == 2) {
				if($weight1 == $weight2) {
				    my $maxV = max($path[$i-2],$path[$i-1]);
				    $maxWeighted{$maxV} = 1;
				    $i--;
				    last;
				}
			    } 
			    if($i == @path) {
				if($weight1 == $weight2) { # if w(i-1) == w(i)
				    if($maxSet[$i] == $maxSet[$i-1]) {
					if($path[$i-1] > $path[$i-2]) {					    
					    $maxWeighted{$path[$i-1]} = 1;
					    $i -= 2;
					    next;
					}
				    }
				}
			    }
			    if($maxSet[$i] ==  $maxSet[$i-1]) {
				$i--;
			    } else {
				print "storing $i: $path[$i-1]\n" if($DEBUG);
				$maxWeighted{$path[$i-1]} = 1;
				$i -= 2;
			    }
			}
			for(my $i=0;$i<@path;$i++) {
			    my $v = $path[$i];
			    unless($maxWeighted{$v}) {
				$knotsRemain = 1;
				$g->delete_vertex($v);
				push(@knotsList,$v);
				print "path-deleted $v\n" if($DEBUG);
			    }
			}
		    }
		} else {
		    # complex graph with more than 2 nodes
		    $knotsRemain = 1;
		    my $minV = "";
		    my $minWeight = 10e10;
		    my $maxDegree = 0;
		    my $maxMaxDegreeScore = -1;
		    my $maxDegreeV;
		    my %nodeInfo;
		    # loop through all vertices of this connected component
		    for(my $i=0;$i<@{$cc};$i++) {
			my $v = $cc->[$i];   # node label
			my $d = $G->degree($v); # degree
			my $weight = @{$segments->[$v]};
			# keep track of the minimum weight node:
			if($weight < $minWeight) {
			    #print "min $v weight = $weight\n";
			    $minWeight = $weight;
			    $minV = $v;
			}
			# and also keep track of the highest degree
			if($d >= $maxDegree) {
			    # when examining a new max degree, reset minMaxDegreeWeight
			    if($d > $maxDegree) {
				$maxMaxDegreeScore = -1;
			    }
			    # compute the sum of the neighbor's weights
			    my $weightSum = 0;
			    for my $w ($g->neighbours($v)) {
				$weightSum += scalar(@{$segments->[$w]});
			    }
			    print "weight sum = $weightSum\n" if($DEBUG);
			    print "v=$v weight = ", scalar(@{$segments->[$v]}), "\n" if($DEBUG);
			    my $score = $weightSum - $weight;
			    $nodeInfo{$v} = [$d,$score,$weight];
			    # when tied for highest degree, take the largest score
			    if($score > $maxMaxDegreeScore) {
				$maxDegreeV = $v;
				$maxDegree = $d;
				$maxMaxDegreeScore = $score;
			    }
			}
		    }
		    my $count = 0;
		    for my $v (keys %nodeInfo) {
			my($d,$score,$weight) = @{$nodeInfo{$v}};
			if(($d == $maxDegree)&&($score == $maxMaxDegreeScore)) {
			    $count++;
			    #print "$v $d $score $weight\n";
			}
		    }
		    #$warnings .= "#Warning: Two nodes found with max degree and max score.\n" if($count > 1);
		    print "two nodes: $id\n" if($DEBUG);
		    if($maxMaxDegreeScore > 0) {
			$g->delete_vertex($maxDegreeV);
			push(@knotsList,$maxDegreeV);
			print "degree-deleted $maxDegreeV\n" if($DEBUG);
			foreach my $pair (@{$segments->[$maxDegreeV]}) {
			    my($f,$t)=@{$pair};
			    print "$maxDegreeV: $f,$t\n" if($DEBUG);
			}
		    } else {
			# otherwise, delete minimum weight node
			$g->delete_vertex($minV);
			push(@knotsList,$minV);
			print "minweight-deleted $minV\n" if($DEBUG);
		    }
		}
	    }
	}
    }
    print "All knots eliminted for this component: ", scalar(@knotsList), " found.\n" if($DEBUG);
    return(\@knotsList,$warnings);
}

sub isPathGraph {
    my($G,$c) = @_;
    # first find two nodes with 1 edge
    my @ends;
    if(scalar(@{$c}) == 1) {
	return 0;
    }
    for my $v (@{$c}) {
	my $d = $G->degree($v);
	if($d) {
	    if($d == 1) {
		push(@ends,$v)
	    }
	} else {
	    print "v=$v\n" if($DEBUG);
	    print $G if($DEBUG);
	    die "No degree for $v\n";
	}
    }
    if(@ends > 2) {
	# pretty sure more than 2 1-degree vertices means it is not a path graph.
	return 0;
    } elsif(@ends < 2) {
	# can't be a path graph
	return 0;
    } else {
	# has 2 ends.
	#@ends = sort {$b <=> $a} @ends;
	my $start = shift(@ends);
	my $end = shift(@ends);
	#print "Graph: $G\n";
	#print "start=$start, end=$end\n";
	my @path;
	push(@path,$start);	
	while(@path < @{$c}) {
	    $,=",";
	    #print "path so far: ", @path, "\n";
	    my @v = $G->neighbours($path[$#path]);
	    if(@v == 2) {
		my($a,$b)=@v;
		#print "found '$a','$b'\n";
		if($a == $path[$#path-1]) {
		    push(@path,$b);
		    #print "storing $b to path\n";
		} elsif($b == $path[$#path-1]) {
		    push(@path,$a);
		    #print "storing $a to path\n";
		} else {
		    die "something else! found that neither '$a' or '$b' was equal to '$path[$#path-1]'\n";
		}
	   } elsif(@v == 1) {
		# if on first one
		if($path[$#path] == $start) {
		    push(@path,$v[$#v]);
		# or if last one
		} elsif($v[$#v] == $end) {
		    push(@path,$v[$#v]);
		} else {
		    #print "neighbor = ", $v[$#v], "\n";
		    die "Unexpected break. found a node with ",scalar(@v)," nbors\n";
		}
	    } else {
		return 0;
	    }
	}
	return \@path;
    }
}

sub min {
    my($x,$y) = @_;
    if($x < $y) {
	return $x;
    } else {
	return $y;
    }
}

sub max {
    my($a,$b) = @_;
    if($b > $a) {
	return $b;
    } 
    return $a;
}

sub getExtremePositions {
    my($bp)= @_;
    if(my @keys = keys(%{$bp})) {
	my $length = @keys;
	my @sortedKeys = sort {$a<=>$b} @keys;
	my $firstKey = $sortedKeys[0];
	my $lastKey = $sortedKeys[$#sortedKeys];    
	return ($firstKey,$lastKey);
    } else {
	die "No basepairs found for $id\n";
    }
}

sub getNextPair {
    my($i,$bp,$lastPos,$knots) = @_;
    for(my $n=$i+1;$n<=$lastPos;$n++) {
	unless(inKnot($n,$knots)) {
	    if($bp->{$n}) {
		return $n;
	    }
	}
    }
    return 0;
}

sub getPrevPair {
    my($j,$bp,$firstPos,$knots) = @_;
    for(my $p=$j-1;$p>=$firstPos;$p--) {
	unless(inKnot($p,$knots)) {
	    if($bp->{$p}) {
		return $p;
	    }
	}
    }
    return 0;
}

sub nonCanonical {
    my($b1,$b2) = @_;
    my @b = sort ($b1,$b2);
    if(($b[0] eq "A")&&($b[1] eq "U")) {
	return 0;
    } 
    if(($b[0] eq "C")&&($b[1] eq "G")) {
	return 0;
    } 
    if(($b[0] eq "G")&&($b[1] eq "U")) {
	return 0;
    } 
    return 1;
}

sub inKnot {
    my($pos,$knots) = @_;
    if($knots) {
	for my $knot (@{$knots}) {
	    # treat $knot like a segment
	    #
	    # 5'start   3'stop
	    #       ACGUA
	    #       |||||
	    #       UGCAU
	    # 3'start   5'stop
	    #	    
	    my @knot = @{$knot};
	    my $first = shift(@knot);
	    my($k_5pStart,$k_3pStart) = @{$first};
	    my $last = @knot ? pop(@knot) : $first;
	    my($k_3pStop,$k_5pStop) = @{$last};
	    if(($k_5pStart <= $pos)&&($pos <= $k_3pStop)) {
		return 1;
	    }
	    if(($k_5pStop <= $pos)&&($pos <= $k_3pStart)) {
		return 1;
	    }
	}
    }
    return 0;
}

sub includesKnot {
    my($start,$stop,$knots) = @_;
    my @loopKnots;
    if($knots) {
	for(my $k=0;$k<@{$knots};$k++) {
	    # treat $knot like a segment
	    #
	    # 5'start   3'stop
	    #       ACGUA
	    #       |||||
	    #       UGCAU
	    # 3'start   5'stop
	    #	    
	    my @knot = @{$knots->[$k]};
	    my $first = shift(@knot);
	    my($k_5pStart,$k_3pStart) = @{$first};
	    my $last = @knot ? pop(@knot) : $first;
	    my($k_3pStop,$k_5pStop) = @{$last};
	    if(($start <= $k_5pStart)&&($k_3pStop <= $stop)) {
		# store the node label, which is $k+1, 1-based
		push(@loopKnots,$k+1);
	    }
	    if(($start <= $k_5pStop)&&($k_3pStart <= $stop)) {
		# store the node label, which is $k+1, 1-based
		push(@loopKnots,$k+1);
	    }
	}
    }
    return @loopKnots;
}

sub pairMap {
    my ($dotbracket) = @_;
    my @map;
    my %stack;
    for (my $i=0; $i<length($dotbracket); $i++) {
	my $c = substr($dotbracket,$i,1);
        if($c =~ /[\(\[\<\{A-Z]/){
	    # left symbols
            push(@{$stack{$c}},$i);
        } elsif($c =~ /[\)\]\>\}a-z]/){
	    # right symbols
            my $C = charMap($c);
            my $pos = pop(@{$stack{$C}});
	    # 0-based to 1-based. See print above.
            $map[$pos] = $i+1;
            $map[$i] = $pos+1;
	    #print "$pos ", $i+1, "\n";
        } elsif($c =~ /[-_,:.]/) {
            $map[$i] = "0";
        } else {
            die "Unknown character $c found in\n$dotbracket\n";
        }
    }
    return @map;
}

sub charMap {
    my($c) = @_;
    if($c eq ")") {
	return "(";
    } elsif($c eq "]") {
        return "[";
    } elsif($c eq ">") {
        return "<";
    } elsif($c eq "}") {
        return "{";
    } elsif($c =~ /[a-z]/) {
        return uc($c);
    } else {
        die "undefined character in charMap: $c\n";
    }
}

#####################
# INPUT SUBROUTINES #
#####################

sub isBPSeqFile {
    my($inputFile) = @_;
    open(BPS,$inputFile) or die "Could not open $inputFile\n";
    while(<BPS>) {
	unless(/^#/) {
	    chomp;
	    if(scalar(split()) != 3) {
		close(BPS);
		return 0
	    }
	}
    }
    close(BPS);
    return 1;
}

sub isDotBracketFile {
    my($inputFile) = @_;
    open(IN,$inputFile) or die "Could not open $inputFile for reading\n";
    my($defline,$sequence,$dotbracket);
    my @lines;
    while(<IN>) {
	unless(/^#/) {
	    chomp;
	    if($_) {
		push(@lines,$_)
	    }
	}
    }
    if(@lines == 2) {
	$sequence = $lines[0];
	$dotbracket = $lines[1];
    } elsif(@lines == 3) {
	$defline = $lines[0];
	$sequence = $lines[1];
	$dotbracket = $lines[2];
    } else {
	return 0;
    }
    if($defline) {
	unless(substr($defline,0,1) eq ">") {
	    return 0;
	}
    }
    if(length($sequence) == length($dotbracket)) {
	return 1;
    } else {
	return 0;
    }
}

sub readBPSeqFile {
    my($bpseqFile) = @_;
    my %bp;
    my %bpCheck;
    my $seq = "";
    open(BPS,$bpseqFile) or die "Could not open $bpseqFile\n";
    while(<BPS>) {
	unless(/^#/) {
	    chomp;
	    die "bad data (need only 3 columns) on line $. of $bpseqFile:\n$_\n" if(scalar(split()) != 3);
	    my($i,$b,$j) = $_ =~ /(\d+)\s([a-zA-Z~\._])\s(\d+)/;
	    $seq .= $b;
	    if(defined($bp{$i})) {
		die "Fatal error: Position $i is paired to both $bp{$i} and $j: Line $. of $bpseqFile\n";
	    }
	    if($bpCheck{$j}) {
		die "Fatal error: Position $j is paired to both $bpCheck{$j} and $i: Line $. of $bpseqFile.\n";
	    }
	    if($j) {
		die "Fatal error: Position $i is paired to itself in $bpseqFile: Line $. of $bpseqFile\n" if($i == $j);
		$bp{$i}=$j;
		$bpCheck{$j}=$i;
	    }
	    if($bpCheck{$i}) {
		die "Fatal error: bpseq file at positions $i paired to $bpCheck{$i} and $j. Caught on line $. of $bpseqFile\n" unless($j == $bpCheck{$i});
	    }	    
	}
    }
    return(\%bp,$seq);
}

sub readDotBracketFile {
    my($dotbracketFile)=@_;
    open(IN,$dotbracketFile) or die "Could not open $dotbracketFile for reading\n";
    my($defline,$sequence,$dotbracket);
    my @lines;
    while(<IN>) {
	unless(/^#/) {
	    chomp;
	    if($_) {
		push(@lines,$_)
	    }
	}
    }
    if(@lines == 2) {
	$sequence = $lines[0];
	$dotbracket = $lines[1];
    } elsif(@lines == 3) {
	$defline = $lines[0];
	$sequence = $lines[1];
	$dotbracket = $lines[2];
    }
    close(IN);
    my @map=pairMap($dotbracket);
    my %bp;
    for(my $i=0;$i<@map;$i++) {
	$bp{$i+1} = $map[$i];
    }
    return(\%bp,$sequence);
}

sub dotBracketToStructureArray {
    my($seq,$dotbracket) = @_;

    # if there are basepairs
    my %bp;
    my @map=pairMap($dotbracket);
    for(my $i=0;$i<@map;$i++) {
	$bp{$i+1} = $map[$i];
    }
    
    if(keys %bp) {
	my $allSegments = getSegments(\%bp);
	my($segments,$knots,$warnings) = separateSegments($allSegments);
	my $bp = filterBasePairs(\%bp,$knots);
	$segments = getSegments($bp);
	my($dotbracket,$s,$k,$structureTypes,$pageNumber) = buildStructureMap($segments,$knots,$bp,$seq);
	$s = join("",@{$s});
	return $s
    } else {
	# default to all external loops
	my $s = "E" x length($seq);
	return $s
    }
}
