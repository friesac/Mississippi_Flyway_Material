# Involved in the blastparse_miss -> blastparse_parsetime -> distance -> blastparse_constel -> blastparse_constel_meta -> elim_dupl_comparisons sequence of scripts
# associate with similarity comparisons of AIV in a region - NOTE: Several intermediate editing steps to files must be conducted prior to each script.
##############################################################################################################################################
## This block is for removing those instances from the PARSED constellation file where two isolates were matched on the same day and the duplicates have not been
## removed because Time2 is not greater than Time1 which is inherent in our "Gene Flow" criteria for this paper.
##############################################################################################################################################
## MAKE SURE YOUR PARSED FILE IS SORTED ALPHABETICALLY BY THE FIRST ISOLATE COLUMN OR ELSE THE HASH WILL HAVE AN EXISTING KEY THAT YOU STILL WANT!

use strict;
use warnings;

my @record1;
my $match1;
my $match2;
my $OUTPUT1;
my %duplicates;

open (PARSEOUTDUPLICATES, "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/ALL_1nt_IRD/ALL_1nt_IRD_constgeneflow_events_PARSED.txt")||die "cannot open file"; # The file that was created above
open ($OUTPUT1, ">", "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/ALL_1nt_IRD/ALL_1nt_IRD_constgeneflow_events_PARSED_dupl_removed.txt")||die "cannot open file";

while (<PARSEOUTDUPLICATES>){ # Enter into line by line to the PARSED file that contains all similarity events by segment, time, dist and number of constellation shares
	chomp; #eliminate newline
	@record1 = split(/\t/, $_); # parse the first line into an array by tab delimiter
	$match1 = $record1[2]; # Assign the first isolate
	$match2 = $record1[10]; # Assign the second isolate
	
	if ($record1[19] == 0) { # If the time column is 0 then we may have to seperate out duplicates 
		$duplicates{"$match2$match1"} = 1; # This places the OPPOSITE comparison that is likely to be seen later in the file into a key for the hash... this key is likely to be written over several times
	}
	my $uniquematch = "$match1$match2"; # Define the appropriate key which would represent the original comparison... the OPPOSITE is defined in the if statement above
	
	print $OUTPUT1 "$_\n" if not exists $duplicates {$uniquematch} ; # If the appropriate key does not exist in the hash then print that line... if it does exist then don't print.
}