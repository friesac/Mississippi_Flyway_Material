# Involved in the blastparse_miss -> blastparse_parsetime -> distance -> blastparse_constel -> blastparse_constel_meta -> elim_dupl_comparisons sequence of scripts
# associate with similarity comparisons of AIV in a region - NOTE: Several intermediate editing steps to files must be conducted prior to each script.
##############################################################################################################################################
## This script should take all observed geneflow events typically only done on a segment of the AIV virus; Then will make a new file that contains only the geneflow events
## that include viruses that were sampled on the same day or a later day (i.e. isolate 1 was sampled before isolate 2. This eliminates duplicates that are the opposite way
## The elim_dupl_comparisons.pl eliminates those comparisons that happen on the same day.  That script is run seperately later on in the analysis.
##############################################################################################################################################
## By the end of this file then you can take that lat long data into a new txt file and run it through the standard distance.pl

use strict;
use warnings;

my $line;
my @timestocompare;
my $time1;
my $time2;
my $time3;
my $OUTPUT3;


open (GENEFLOWPARSETIME, "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_IRDBlast_geneflowbytimeandspace_2nt_IRD.txt")||die "cannot open file";
open ($OUTPUT3, ">", "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_geneflowfinal_2nt.txt")||die "cannot open file";

while (<GENEFLOWPARSETIME>){ # Enter into the file that contains ALL similarity events regardless of timing
	$line = $_; # assign the first line
	@timestocompare = split(/\t/, $_); # split the line into an array by tab delimiter
	
	if ($timestocompare[18] =~ m/IRD/) {
		print $OUTPUT3 "$line";
	}
	
	else {
	$time1 = $timestocompare[6]; # Assign isolation date 1
	$time2 = $timestocompare[14]; # Assign isolation date 2
	$time3 = $time2-$time1; # Assign the difference between the two events
		
		if ($time2 >= $time1) { #Only keep those events that have a >= 0 comparison
			print $OUTPUT3 "$line"; # Put those comparisons into the geneflow final txt file.
		}	
	}
}
