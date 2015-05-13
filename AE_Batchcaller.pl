#! /usr/bin/perl -w
use strict;
use IPC::Open3;
use IO::Select;
use JSON qw( decode_json );

my ($out, $err) = ( '', '' );

my $loci_f = shift;
my %loci = ();
my @loci_arr = ();
my %gene = ();
my %transcript = ();

&read_locilist ($loci_f, \@loci_arr, \%loci, \%gene, \%transcript);

my $curr_loci = shift(@loci_arr); 

print $curr_loci." curr_loci ".$gene{"chr17 78064145"}."\n";

# my ($chr, $begin, $end) = ("chr1", 100, 200);
# $curr_loci = "chr1 201";
# print "($chr, $begin, $end $curr_loci) ".&isMatch($chr, $begin, $end, $curr_loci)."\n";

# ($chr, $begin, $end) = ("chr1", 100, 200);
# $curr_loci = "chr1 200";
# print "($chr, $begin, $end $curr_loci) ".&isMatch($chr, $begin, $end, $curr_loci)."\n";

# ($chr, $begin, $end) = ("chr1", 100, 200);
# $curr_loci = "chr1 150";
# print "($chr, $begin, $end $curr_loci) ".&isMatch($chr, $begin, $end, $curr_loci)."\n";


# ($chr, $begin, $end) = ("chr1", 100, 200);
# $curr_loci = "chr1 50";
# print "($chr, $begin, $end $curr_loci) ".&isMatch($chr, $begin, $end, $curr_loci)."\n";

# ($chr, $begin, $end) = ("chr1", 100, 200);
# $curr_loci = "chrX 200";
# print "($chr, $begin, $end $curr_loci) ".&isMatch($chr, $begin, $end, $curr_loci)."\n";

# ($chr, $begin, $end) = ("chr1", 100, 200);
# $curr_loci = "chr2 200";
# print "($chr, $begin, $end $curr_loci) ".&isMatch($chr, $begin, $end, $curr_loci)."\n";

# ($chr, $begin, $end) = ("chrX", 100, 200);
# $curr_loci = "chrX 205";
# print "($chr, $begin, $end $curr_loci) ".&isMatch($chr, $begin, $end, $curr_loci)."\n";


my $var_f = shift;
my $AE_path = shift;
print "perl $AE_path/AE_xs.pl\n";
my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, "perl $AE_path/AE_xs.pl")
#my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, 'cat')
    or die "open3() failed $!";

my $pid_var = open3(\*VARCHLD_IN, \*VARCHLD_OUT, \*VARCHLD_ERR, "bzcat $var_f")
#my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, 'cat')
    or die "open3() failed $!";

my $reader = IO::Select->new(\*CHLD_OUT, \*CHLD_ERR);
my $reader_var = IO::Select->new(\*VARCHLD_OUT);

my @readyvar = $reader_var->can_read();
my $fhvar = $readyvar[0];

my $match_f = "match.txt";
my $ann_f = "annotation.txt";
my $rcd_f = "rcd.txt";

my $cnt = 0;
open (M, ">$match_f") || die;
open (ANN, ">$ann_f") || die;
open (RCD, ">$rcd_f") || die;
my %annotated = ();


#4858615	2	1	chr5	77396835	77396838	ref	TTC	TTC	564	563	PASS
#4858615	2	2	chr5	77396835	77396838	del	TTC		564	563	PASS
my $prv_loci = "";
my $line1 = "";
while (<$fhvar>) {
	#print $_;
	if (/^\d+\t+(\d+)\t+(\S+)\t+(chr\S+)\t+(\d+)\t+(\d+)\t(\S*?)\t(\S*?)\t(.*?)\t\d*\t\d*\t(\S*)/) {
		$cnt ++;
		my ($ploidy, $allele, $chr, $begin, $end, $type, $ref, $var, $qual) = ($1, $2, $3, $4, $5, $6, $7, $8, $9);
		
		if ($cnt%10000 == 0) {
			print "$cnt finished ($chr, $begin, $end, $ref, $var)";
		}

		#if ($end != 35683240) {
		#	next;
		#}

		print "line $_\n";

		my $match_flag = &isMatch($chr, $begin, $end, $curr_loci);
		if ($match_flag > 0) {
			next;
		} elsif ($match_flag < 0) {
			if (@loci_arr == 0) {
				last;
			}

			while ($match_flag < 0) {
				print "($chr, $begin, $end, $curr_loci) $match_flag bf\n";
				if (@loci_arr > 0) {
					$curr_loci = shift(@loci_arr);
				} else {
					last;
				}
				$match_flag = &isMatch($chr, $begin, $end, $curr_loci);
				print "($chr, $begin, $end, $curr_loci) $match_flag af\n";
			}
		}

		if ($match_flag == 0) {
			print M $_;
			print "match ($chr, $begin, $end, $ref, $var)\n";

			print CHLD_IN &inputjson($chr, $begin, $end, $ref, $var);
			my $tmpline = <CHLD_OUT>;

			my $decoded_json = decode_json($tmpline);

			if ($allele eq 'all') {
				if (($ref eq '=') && ($var eq '=')) {
					printf ANN "$chr\t$begin\t$end\t$ref\t$var\tref\t$loci{$curr_loci}\t$gene{$loci{$curr_loci}}\tN/A\t$transcript{$loci{$curr_loci}}\n";
				} elsif ($var eq '?') {
					print "NO CALL $curr_loci\n";
					print ANN "$chr\t$begin\t$end\t$ref\t$var\t$type\t$loci{$curr_loci}\t$gene{$loci{$curr_loci}}\tN/A\t$transcript{$loci{$curr_loci}}\n";
				} else {
					print "Impossible allele all\n";
				}
			} elsif (($allele == 1) && ($ploidy == 1)) {
				my $output_line = &printAnnotationLine($chr, $begin, $end, $ref, $var, $type, $qual,
														\%loci, \%gene, \%transcript, 
														$curr_loci, 
														$decoded_json);
				print ANN $output_line."\n";
			} elsif (($allele == 1) && ($ploidy == 2)) {
				$prv_loci = "$chr $begin $end $ref $var $type";
				$line1 = '';

				if ($ref eq $var) {
					next;
				} else {
					$line1 = &printAnnotationLine($chr, $begin, $end, $ref, $var, $type, $qual,
														\%loci, \%gene, \%transcript, 
														$curr_loci, 
														$decoded_json);

					#print "line1 $line1\n";
					next;
				}
			} elsif (($allele == 2) && ($ploidy == 2)) {
				my ($prv_chr, $prv_begin, $prv_end, $prv_ref, $prv_var, $prv_type) = split (/\s/, $prv_loci);
				my $zygosity = '';

				#print "pre ($prv_chr, $prv_begin, $prv_end, $prv_ref, $prv_var, $prv_type) $prv_loci"."\n";

				if (($prv_chr eq $chr) && ($prv_begin == $begin) && ($prv_end == $end)) {
					if ($prv_var eq $var) {
						$zygosity = "Homozygous";
						print ANN $line1."\t".$zygosity."\n";
					} else {
						$zygosity = "Heterozygous";
						my $line2 = ''; 

						if ($ref ne $var) {
							$line2 = &printAnnotationLine($chr, $begin, $end, $ref, $var, $type, $qual,
														\%loci, \%gene, \%transcript, 
														$curr_loci, 
														$decoded_json);
						}
							
						if ($line1 ne '') {
							print ANN $line1."\t".$zygosity."\n";
						}
						
						if ($line2 ne '') {
							print ANN $line2."\t".$zygosity."\n";
						}
					}
				}	
			} else {
				print "Impossible\n";
			}
		}
	}
#	exit;
}
close (RCD);
close (ANN);
close (M);
exit;

sub printAnnotationLine() {
	my ($chr, $begin, $end, $ref, $var, $type, $qual, $loci_r, $gene_r, $transcript_r, $curr_loci, $decoded_json) = @_;
	my $rtn = '';

	if ($type eq 'no-call') {
		$rtn = sprintf("%s", "$chr\t$begin\t$end\t$ref\t$var\t$type\t$loci_r->{$curr_loci}\t$gene_r->{$loci{$curr_loci}}\t$transcript_r->{$loci{$curr_loci}}\tno-call\t$qual");
	} else {
		my $var_name = '';
		my $tx_orientation = '';
		my $rpt_varname = $var_name = $decoded_json->{"var"}{"VarName"};
		if (exists ($decoded_json->{"trInfo"}{$transcript_r->{$loci_r->{$curr_loci}}})) {
			$var_name = $decoded_json->{"trInfo"}{$transcript_r->{$loci_r->{$curr_loci}}}{"TranscriptVarName"};
			$tx_orientation = $decoded_json->{"trInfo"}{$transcript_r->{$loci_r->{$curr_loci}}}{"TranscriptOrientation"};
		} else {
			$var_name = $decoded_json->{"var"}{"VarName"};
			print "varname $curr_loci $var_name\n";
			my $txname = (split (/\s*\(/, $var_name))[0];
			$tx_orientation = $decoded_json->{"trInfo"}{$txname}{"TranscriptOrientation"};
			print "varname $var_name $txname $tx_orientation \n";
		}
		
		if ($var_name eq $rpt_varname) {
			$rpt_varname = 'same';
		}

		my $tenK_AF = "N/A";
		if (exists ($decoded_json->{"var"}{"1000G_AF"})) {
			$tenK_AF = $decoded_json->{"var"}{"1000G_AF"};
		}

		my $dbsnpIds = "N/A";
		if (exists ($decoded_json->{"var"}{"dbsnpIds"})) {
			$dbsnpIds = $decoded_json->{"var"}{"dbsnpIds"};
		}

		$rtn = sprintf("%s", "$chr\t$begin\t$end\t$ref\t$var\t$type\t$loci_r->{$curr_loci}\t$gene_r->{$loci{$curr_loci}}\t$qual\t$transcript_r->{$loci{$curr_loci}}\t$var_name\t$rpt_varname\t$tx_orientation\t$tenK_AF\t$dbsnpIds");
	}

	return $rtn;
}

sub isMatch () {
	my ($chr, $begin, $end, $curr_loci) = @_;
	my $rtn = 0;

	my ($loci_chr, $loci_end) = split (/\s+/, $curr_loci);

	$chr =~ s/chrX/chr23/g;
	$loci_chr =~ s/chrX/chr23/g;
	$chr =~ s/chr//g;
	$loci_chr =~ s/chr//g;

	if ($chr < $loci_chr) {
		$rtn = 1;
	} elsif ($chr > $loci_chr) {
		$rtn = -1;
	} else {
		if ($end < $loci_end) {
			$rtn = 1;
		} elsif ($end == $loci_end) {
			$rtn = 0;
		} else {
			if ($begin >= $loci_end) {
				$rtn = -1;
			} else {
				$rtn = 0;
			}
		}
	}

	return $rtn;
}

exit;


sub parse_varfileLine () {
	my ($line) = @_;
	my ($chr, $begin, $end, $ref, $var);

	if ($line =~ /^\d+\s+\d+\s+\S+\s+(chr\S+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
		($chr, $begin, $end, $ref, $var) = ($1, $2, $3, $4, $5);
	}

	return ($chr, $begin, $end, $ref, $var);
}


sub read_locilist () {
	my ($loci_f, $loci_arr_r, $loci_r, $gene_r, $transcript_r) = @_;

	open (LOCI, "<$loci_f") || die;
	while (<LOCI>) {
		chomp;

		if (/^(\S+)\t(HGNC\:\d+)\t(chr\S+)\t(\d+)\t(\S+)/) {
			my ($gene, $HGNC, $chr, $end, $tx) = ($1, $2, $3, $4, $5);
			$loci_r->{$chr." ".$end} = $gene;
			push (@$loci_arr_r, $chr." ".$end);
			$gene_r->{$gene} = $HGNC;
			$transcript_r->{$gene} = $tx;
		}
	}
	close (LOCI);
}

sub inputjson () {
	my ($chr, $begin, $end, $ref, $var) = @_;

	my $rtn = "\'\{\"chr\":\"$chr\",\"begin\":\"$begin\",\"end\":\"$end\",\"referenceSequence\":\"$ref\",\"variantSequence\":\"$var\"\}\'\n";

	return $rtn;
}

