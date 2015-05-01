#! /usr/bin/perl -w
use strict;
use IPC::Open3;
use IO::Select;
use JSON qw( decode_json );

my ($out, $err) = ( '', '' );

my $AE_path = shift;
print "perl $AE_path/AE_xs.pl\n";
my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, 'perl AE_xs.pl')
#my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, 'cat')
    or die "open3() failed $!";

my $r;

my $reader = IO::Select->new(\*CHLD_OUT, \*CHLD_ERR);

	while ( my @ready = $reader->can_read()) {
		print CHLD_IN "\'\{\"chr\":\"chr8\",\"begin\":\"24811065\",\"end\":\"24811066\",\"referenceSequence\":\"G\",\"variantSequence\":\"A\"\}\'\n";
		
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

waitpid($pid, 1);

