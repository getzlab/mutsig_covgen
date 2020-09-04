#!/usr/bin/perl
use POSIX;

#2016-08-12//Julian Hess

#flush buffer to disk; reserve partially filled overhanging byte, if it exists
sub flush {
	my $buf = \@{$_[0]}; my $OUTFWB = $_[1];

	my $end = $#{$buf} - ($_[2] ? 0 : ($#{$buf} + 1) % 8);

	print $OUTFWB pack("B*", join("", @{$buf}[0 .. $end])); 
	print "Flushed " . ($end + 1) . " bits.\n";

	@{$buf} = @{$buf}[$end + 1 .. $#{$buf}];
	#splice(@{$buf}, 0, $end + 1);  #is this faster than slicing?
}

if($#ARGV != 1) {
	print "Usage: mutect_wig2fwb.pl <input wig> <output FWB stem>\n";
	exit 1;
}

$inwig = $ARGV[0];
$outfwb = $ARGV[1];

if($outfwb =~ /\.fwb$/i) {
	$outfwb =~ s/\.fwb//;
}

%chrH = ( X => 23, Y => 24);

open(INWIG, $inwig) || die "Cannot open input wig $inwig";
open($OUTFWBFH, ">", $outfwb . ".fwb") || die "Cannot open output FWB $outfwb.fwb";
open(OUTFWI, ">", $outfwb . ".fwi") || die "Cannot open output FWI $outfwb.fwi";

$chr = 0; $start = 0;
while(<INWIG>) {
	next if($_ =~ /^track/ || $_ eq "\n");

	if($_ =~ /^fixedStep chrom=([^ ]+) start=(\d+)/) {
		$contigM = $1; $startM = $2;
		next if($contigM !~ /^[\dXY]+$/);

		$pchr = $chr; $pstart = $start;
		$chr = $contigM; $start = $startM;

		$chr = $chrH{$chr} if($chr =~ /^[XY]$/);

		#handle the first interval
		if(!$fflag) {
			print OUTFWI $chr . "\t" . $start;
			$pchr = $chr; $pstart = $start; $ostart = $start;
			$fflag = 1;
		}

		#if we've skipped more than 100kb or switched chromosomes
		if($chr != $pchr || $start - ($pstart + $nchunk) > 1e5) {
			print OUTFWI "\t" . ($ostart + $nwritten - 1) . "\n";
			print OUTFWI $chr . "\t" . $start;
			print $chr . "\t" . $start . "\n";
			$ostart = $start;
			$nwritten = 0;
		} else { #pad out with zeros
			$npad = $start - ($pstart + $nchunk);
			push(@buf, (0) x $npad);
			$nwritten += $npad;
		}
		$nchunk = 0;

		flush(\@buf, $OUTFWBFH) if($#buf > 10e6);
	} else {
		chomp;
		push(@buf, $_); 
		$nwritten++; $nchunk++;
	}
}

#final buffer flush
print OUTFWI "\t" . ($ostart + $nwritten - 1) . "\n"; 
flush(\@buf, $OUTFWBFH, 1);
