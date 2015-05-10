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
my $output_line = "";
while (<$fhvar>) {
	#print $_;
	if (/^\d+\t+\d+\t+(\S+)\t+(chr\S+)\t+(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t(.*?)\t/) {
		$cnt ++;
		my ($allele, $chr, $begin, $end, $type, $ref, $var) = ($1, $2, $3, $4, $5, $6, $7);
		
		if ($cnt%10000 == 0) {
			print "$cnt finished ($chr, $begin, $end, $ref, $var)";
		}

		if ($cnt == 5) {
			#exit;
		}

		#print "line $_\n";

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

			my $zygosity = '';
			if ($prv_loci eq '') {
				$prv_loci = "$chr $begin $end $ref $var $type";
			} else {
				my ($prv_chr, $prv_begin, $prv_end, $prv_ref, $prv_var, $prv_type) = split (/\s+/, $prv_loci);

				print "($prv_chr, $prv_begin, $prv_end, $prv_ref, $prv_var, $prv_type)"."\n";

				if (($prv_chr eq $chr) && ($prv_begin == $begin) && ($prv_end == $end)) {
					if ($prv_var eq $var) {
						$zygosity = "Homozygous";
						print ANN $output_line."\t".$zygosity."\n";
					} else {
						$zygosity = "Heterozygous";
						if ($type ne 'ref') {
							my $tmp = sprintf("%s", "$chr\t$begin\t$end\t$ref\t$var\t$type\t$loci{$curr_loci}\t$gene{$loci{$curr_loci}}\t$transcript{$loci{$curr_loci}}\t".$decoded_json->{"var"}{"VarName"});
							print ANN $tmp."\t".$zygosity."\n";
						} else {
							print ANN $output_line."\t".$zygosity."\n";
						}
					}
					$prv_loci = '';
					$output_line = '';
				}
			}

			if ($ref eq $var) {
				if (($ref eq '=') && ($var eq '=')) {
					printf ANN "$chr\t$begin\t$end\t$ref\t$var\tref\t$loci{$curr_loci}\t$gene{$loci{$curr_loci}}\t$transcript{$loci{$curr_loci}}\n";
					$prv_loci = '';
					$output_line = '';
				} else {
					print "here a\n";
					$output_line = '';
				}
				next;
			} elsif ($var eq "?") {
				print "NO CALL $curr_loci\n";
				print ANN "$chr\t$begin\t$end\t$ref\t$var\t$type\t$loci{$curr_loci}\t$gene{$loci{$curr_loci}}\t$transcript{$loci{$curr_loci}}\n";
				next;
			} else {

			#if (exists $annotated {$decoded_json->{"var"}{"VarName"}}) {
			#	next;
			#} else {
			#	$annotated {$decoded_json->{"var"}{"VarName"}} = $decoded_json;
				print "snp\n";
				$output_line = sprintf("%s", "$chr\t$begin\t$end\t$ref\t$var\t$type\t$loci{$curr_loci}\t$gene{$loci{$curr_loci}}\t$transcript{$loci{$curr_loci}}\t".$decoded_json->{"var"}{"VarName"});

				if ($chr eq 'chrX') {
					print ANN $output_line."\n";
					$prv_loci = '';
					$output_line = '';
				}
			}
		}
	}	
#	exit;
}
close (RCD);
close (ANN);
close (M);
exit;

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
# 	print "here\n";
# 	foreach my $fh (@ready) {
# 		my $line = <$fh>;
# 		if (!defined $line) {
# 			$reader->remove($fh);
# 			$fh->close();
# 		} else {

			
# 			if (fileno($fh) == fileno(\*VARCHLD_OUT)) {
# 				if ($line =~ /^\d+\s+\d+\s+\S+\s+(chr\S+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
# 					my ($chr, $begin, $end, $ref, $var) = ($1, $2, $3, $4, $5);
# 					print "($chr, $begin, $end, $ref, $var)\n";
# 				}
# 			}
# 		}
# 	}
# }


# while ( my @ready = $reader_var->can_read() ) {
# 	print "here\n";
# 	foreach my $fh (@ready) {
# 		my $line = <$fh>;
# 		if (!defined $line) {
# 			$reader->remove($fh);
# 			$fh->close();
# 		} else {


# 			if (fileno($fh) == fileno(\*VARCHLD_OUT)) {
# 				if ($line =~ /^\d+\s+\d+\s+\S+\s+(chr\S+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
# 					my ($chr, $begin, $end, $ref, $var) = ($1, $2, $3, $4, $5);
# 					print "($chr, $begin, $end, $ref, $var)\n";
# 				}
# 			}
# 		}
# 	}
# }

exit;

for(my $i=1;$i<12;$i++) {
	print "$i\n";
	print CHLD_IN "\'\{\"chr\":\"chr8\",\"begin\":\"24811065\",\"end\":\"24811066\",\"referenceSequence\":\"G\",\"variantSequence\":\"A\"\}\'\n";
	my $cnt = 0;
	while ( my @ready = $reader->can_read() ) {
		$cnt++;
		print "$cnt\n";
		foreach my $fh (@ready) {
			my $line = <$fh>;
			if (!defined $line) {
				$reader->remove($fh);
				$fh->close();
			} else {
				print fileno(\*CHLD_OUT)."  id\n";
				print fileno($fh)."  fh id\n";
				if (fileno($fh) == fileno(\*CHLD_OUT)) {
					$out .= $line;

					my $decoded = decode_json($line);
					print $decoded->{'var'}{'chr'}."  chr \n";

					print "here\n";
					print STDOUT $line;
				}
				elsif (fileno($fh) == fileno(\*CHLD_ERR)) {
					$err .= $line;
					print STDERR $line;
				}
			}
		}
		last;
	}	
	print "$i last\n";
	#print CHLD_IN "$i\n";

    #print "Got $r from child\n";
}
waitpid($pid, 1);


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

