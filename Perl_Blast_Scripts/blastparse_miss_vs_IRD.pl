## Involved in the blastparse_miss -> blastparse_parsetime -> distance -> blastparse_constel -> blastparse_constel_meta -> elim_dupl_comparisons sequence of scripts
## associate with similarity comparisons of AIV in a region - NOTE: Several intermediate editing steps to files must be conducted prior to each script.
##############################################################################################################################################
## This script is designed to take the output from a TAB DELIMITED blast output and parse  out and include two more columns that contain the states from each isolate
## (in the first block of code) and then go into our sequence database of all isolates and extract out the metadata that we are interested in... here that will be time, lat, long,
## latitudinal bin and isolate names.
##############################################################################################################################################
## This file is where you must make sure that the geneflow_events file with states included is the same length of the by times and state file... this is a problem because
## our regex code is not perfect. You will have to change some names in the blast result files (i.e. 9999) to account for this discrepancy. Then run this script again until perfect.

use strict;
use warnings FATAL => 'uninitialized';

my $line;
my @record;
my $match1;
my @match1parse;
my $size1;
my $state1;
my $match2;
my @match2parse;
my $size2;
my $state2;
my $identity;
my $length;
my $OUTPUT1;
my @timestoget;
my $isolate1;
my $isolate2;
my $timematch1;
my $timematch2;
my $isoline;
my @isolatedb;
my $OUTPUT2;
my $id1;
my @id2parse;
my $id2;
my @NCBIrecord;
my $NCBIdate;
my $NCBIstate;
my $NCBIcountry;
my $NCBIsubtype;
my $OUTPUT3;
my $ncbiline;
my @ncbilinerecord;
my $OUTPUT4;
my $bit;
my @bitparse;
my $NCBIConvert;
my $ncbimatch;

open (BLASTRESULTS, "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/mp_irdblast_results_2nt.txt")||die "cannot open file";
open ($OUTPUT1, ">" . "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_IRDBlast_geneflow_events_2nt.txt")||die "cannot open file";
open ($OUTPUT4, ">" . "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_Blastparse_same_isolate_comparison_2nt.txt");

while (<BLASTRESULTS>){ #getting the state distinction for each pairwise blast comparison and eliminating same comparisons
	chomp;
	$line = $_; #grab the first line of results
	@record = split(/\t/, $_); # parse the first line into an array by tab delimiter
	$match1 = $record[0]; # assign first isolate name to variable
	@match1parse = split(/_/,$match1); # split the isolate name by underscore delimiter
	$size1=scalar(@match1parse); # acquire the size of the new 1st isolate array
	$state1= $match1parse[$size1-4]; # in this array the 4th to last element should be the state
	$match2 = $record[1]; # assign second isolate name to variable
	@match2parse = split(/\|/,$match2); # split the isolate name from IRD based on "|" delimiter
	$size2=scalar(@match2parse); # acquire the size of the new 2nd isolate array
	$state2= $match2parse[$size2-5]; # in this array the 4th to last element should be the state also... careful, this is only true if using our naming system
		if ($state2 eq 'NA') { # this is if the isolate is from a different country
			$state2 = $match2parse[$size2-4];
		}
	$identity = $record[2]; # record the percent identity from the blast results
	$length = $record[3]; # record the length of the alignment
	$id1 = $match1parse[$size1-3]; # This should grab the OSU accession number... i.e. 10OS4750
	@id2parse = split(/\//, $match2parse[1]); #This should grab the second element from IRD which should be the strain name
	$id2 = $id2parse[3]; # This should grab the OSU accession number from the strain name.
	
		if ($id1 ne $id2) { # i.e. If the first isolate is not the same as the second WITHOUT THE STATE CONSARISON
		print $OUTPUT1 "$state1\t$state2\t$match1\t$match2\t$identity\t$length\n"; # then print state;state;isolate1;isolate2;percentident;length
		}
		elsif ($id1 eq $id2) {
			print $OUTPUT4 "$id1\t$id2\t$record[1]\n"; # the same isolate list... check this to make sure all of your MS Flyway samples are in your IRD Database... also has the IRD name to add to the geneflow file
		}
	
}	#End this loop

print "done with upper loop\n";

close $OUTPUT1;
close $OUTPUT4;
close BLASTRESULTS;

########################################################################################################################################################
##This block is the block that will go through the events folder and add genbank names to the MS flyway samples for further processing
########################################################################################################################################################
open (GENEFLOW, "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_IRDBlast_geneflow_events_2nt.txt")||die "cannot open file"; # open up the geneflow file you just made
open ($OUTPUT3, ">" . "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_IRDBlast_geneflow_events_genbanknames_2nt.txt")||die "cannot open file"; 

while (<GENEFLOW>){ # enter into the geneflow file that has all isolates that aren't identical 
		chomp;
		$line = $_; # assign the first line
		@timestoget = split(/\t/, $_); # split the line into an array by tab delimiter
		$isolate1 = $timestoget[2]; # the 2nd element should be the first isolate
		$isolate1 =~ /(\S{8}_[1:2][9][0-9][0-9]|\S{8}_[2][0][0-9][0-9])/; # match the id of the isolate
		$timematch1 = $1; # assign the id of the isolate
		$isolate1 =~ /(\A)/;
		
		open (SAME, "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_Blastparse_same_isolate_comparison_2nt.txt")||die "cannot open file"; # open up the duplicates file

		while (<SAME>) {	# enter into the identical isolates file to grab genbank names to add to the geneflow document
				chomp;
				$bit = $_; # assign line
				@bitparse = split(/\t/, $_);	# split on tab
				$NCBIConvert = $bitparse[2]; 	# grab genbank/IRD name
				$NCBIConvert =~ s/\//_/g;		# replace / with _
				$NCBIConvert =~ /(\S{8}_[1:2][9][0-9][0-9]|\S{8}_[2][0][0-9][0-9])/;
				$ncbimatch = $1; 
				
				if ($ncbimatch eq $timematch1) {
					print $OUTPUT3 "$timestoget[0]\t$timestoget[1]\t$bitparse[2]\t$timestoget[3]\t$timestoget[4]\t$timestoget[5]\n";
					last;
				}
		}
}
close GENEFLOW;
close SAME;
close $OUTPUT3;		

########################################################################################################################################################
##This block is the block referenced in the header of this file that is going to go into our database and extract metadata
########################################################################################################################################################
open (GENEFLOW, "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_IRDBlast_geneflow_events_genbanknames_2nt.txt")||die "cannot open file"; # open up the geneflow file you just made
open (DATABASE, "Tony/Influenza_Genomics/JCVI_OSU_Sequencing_12_12_2013.txt")||die "cannot open file"; # open up the database of isolates
open ($OUTPUT2, ">", "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_IRDBlast_geneflowbytimeandspace_2nt.txt")||die "cannot open file"; # open up a new output file
		
while (<GENEFLOW>){ # enter into the geneflow file that has all isolates that aren't identical 
		chomp;
		$line = $_; # assign the first line
		@timestoget = split(/\t/, $_); # split the line into an array by tab delimiter
		$isolate1 = $timestoget[2]; # the 2nd element should be the first isolate
		$isolate1 =~ s/\//_/g;
		$isolate1 =~ /(\S{8}_[1:2][9][0-9][0-9]|\S{8}_[2][0][0-9][0-9])/; # match the id of the isolate
		$timematch1 = $1; # assign the id of the isolate
		$isolate2 = $timestoget[3]; # the 3rd element should be the second isolate
		$isolate2 =~ s/\//_/g; #This is replacing the "/" that NCBI places into there strain names and dates... will ultimately end up capturing 2 variables but the OSU accession number should be first.
		$isolate2 =~ /(\S{8}_[1:2][9][0-9][0-9]|\S{8}_[2][0][0-9][0-9])/; # match the id of the isolate
		$timematch2 = $1; # assign the id of this isolate
		$isolate2 =~ /(\A)/; # This is to reset the capture variable.. because $1 will be maintained through the next while loop until a successful match is formed.
		my $matchtest = 0;
		

		open (DATABASE, "Tony/Influenza_Genomics/JCVI_OSU_Sequencing_12_12_2013.txt")||die "cannot open file"; # open up the database of isolates
		
		while (<DATABASE>){  # Now you enter into the mainframe of the database line by line... this is the most important file to keep up to date.
			$isoline = $_;  # Assigning the entire line to the $isoline varible
			@isolatedb = split(/\t/, $_); # Putting the specific line of the database (i.e. one isolate record) into an array
			$isoline = /(\S{8}_[1:2][9][0-9][0-9]|\S{8}_[2][0][0-9][0-9])/; # same capture idea as above except $1 is empty before entering the loop
			my $timematch3 = $1; # assign the id of the isolate from the line of the database
						
			if ($timematch1 eq $timematch3){ # if they match each other
				print $OUTPUT2 "$timestoget[4]\t$timestoget[5]\t$timestoget[2]\t$isolatedb[13]\t$isolatedb[6]\t$isolatedb[7]\t$isolatedb[8]\t$isolatedb[9]\t$isolatedb[10]\t$isolatedb[4]\t"; #then print to the new output file
			}
		}	
		
		open (DATABASE, "Tony/Influenza_Genomics/JCVI_OSU_Sequencing_12_12_2013.txt")||die "cannot open file"; # open up the database of isolates
		
		while (<DATABASE>){  # Now you enter into the mainframe of the database line by line... this is the most important file to keep up to date.
			$isoline = $_;  # Assigning the entire line to the $line varible
			@isolatedb = split(/\t/, $_); # Putting the specific line of the database (i.e. one isolate record) into an array
			$isoline = /(\S{8}_[1:2][9][0-9][0-9]|\S{8}_[2][0][0-9][0-9])/; # same capture idea as above except $1 is empty before entering the loop
			my $timematch3 = $1; # assign the id of the isolate from the line of the database

			if ($timematch2 eq $timematch3){ #if they match each other 
				$matchtest++;
				print $OUTPUT2 "$timestoget[3]\t$isolatedb[13]\t$isolatedb[6]\t$isolatedb[7]\t$isolatedb[8]\t$isolatedb[9]\t$isolatedb[10]\t$isolatedb[4]\tUSA\n"; # then print
			}
		}
		if ($matchtest == 0) {
			@NCBIrecord = split(/\|/, $timestoget[3]);
			$NCBIstate = $NCBIrecord[3];
			$NCBIsubtype = $NCBIrecord[6];
			$NCBIdate = $NCBIrecord[7];
			$NCBIcountry = $NCBIrecord[4];
			print $OUTPUT2 "$timestoget[3]\t$NCBIsubtype\tMigYear\t$NCBIdate\tDays\tLat\tLong\t$NCBIstate\t$NCBIcountry\n"; # then print
		}
}

