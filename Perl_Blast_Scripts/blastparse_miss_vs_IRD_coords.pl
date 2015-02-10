## Involved in the blastparse_miss_vs_IRD -> blastparse_miss_vs_IRD_coords -> blastparse_parsetime -> distance -> blastparse_constel -> blastparse_constel_meta -> elim_dupl_comparisons sequence of scripts
## associate with similarity comparisons of AIV in a region - NOTE: Several intermediate editing steps to files must be conducted prior to each script.
##############################################################################################################################################
## This script is designed to take the output from blastparse_miss_vs_IRD.pl in the form of geneflobytimeandspace_%.txt and replace the locations taken from each NCBI record and 
## input latitude and longitude coords for each record. .
##############################################################################################################################################

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
my $OUTPUT3;
my $ncbiline;
my @ncbilinerecord;
close DATABASE;
close GENEFLOW;

open (TIMEANDSPACE, "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_IRDBlast_geneflowbytimeandspace_2nt.txt")||die "cannot open file"; # open up the file you just wrote
open ($OUTPUT3, ">", "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_IRDBlast_geneflowbytimeandspace_2nt_IRD.txt")||die "cannot open file"; #new output file

## my %NCBILocations = (
		## Minnesota => '46.042736\t-94.573975',
		## California => '36.844461\t-119.49463',
		## New_Jersey => '40.061257\t-74.381104',
		## Arkansas => '34.759666\t-92.508545',
		## North_Dakota => '47.620975\t-100.440674',
		## Delaware => '38.976492\t-75.457764',
		## Louisiana => '31.503629\t-92.596436',
		## Texas => '31.484893\t-98.726807',
		## Alaska => '65.403445\t-151.538088',
		## New_Brunswick => '46.649436\t-66.163331',
		## Nova_Scotia => '45.02695\t-63.26294',
		## Alberta => '55.40407\t-114.792481',
		## Quebec => '50.736455\t-73.923341',
		## MB => '53.981935\t-98.137208',
		## Prince_Edward_Island => '46.346928\t-63.366852',
		## Nunavet => '65.622023\t-97.302247',
		## Interior_Alaska => '65.403445\t-151.538088',
		## Guatemala => '15.538376\t-90.069581',
	## );	
	
while (<TIMEANDSPACE>) {
	chomp;
	$ncbiline = $_;
	@ncbilinerecord = split(/\t/, $ncbiline);
	if ($ncbilinerecord[17] eq 'Minnesota') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t46.042736\t-94.573975\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'California') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t36.844461\t-119.49463\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'Indiana') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t39.622615\t-86.167603\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'Arizona') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t34.415973\t-111.66687\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'New_Mexico') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t34.452218\t-105.975952\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'Ohio') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t40.254377\t-82.760925\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'New_Jersey') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t40.061257\t-74.381104\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'Arkansas') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t34.759666\t-92.508545\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'North_Dakota') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t47.620975\t-100.440674\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'South_Dakota') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t44.559163\t-100.482788\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'Delaware') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t38.976492\t-75.457764\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'Louisiana') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t31.503629\t-92.596436\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'Texas') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t31.484893\t-98.726807\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'Alaska') {
		print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t65.403445\t-151.538088\t$ncbilinerecord[17]\tIRD\n";
	}
	elsif ($ncbilinerecord[17] eq 'NA') {
		my @parserecord = split(/\|/, $ncbilinerecord[10]);
		my @parsestrain = split(/\//, $parserecord[1]);
		
		if ($parsestrain[2] eq 'New_Brunswick') {
			print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t46.649436\t-66.163331\t$ncbilinerecord[17]\tIRD\n";
		}
		elsif ($parsestrain[2] eq 'Nova_Scotia') {
			print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t45.02695\t-63.26294\t$ncbilinerecord[17]\tIRD\n";
		}
		elsif ($parsestrain[2] eq 'Alberta') {
			print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t55.40407\t-114.792481\t$ncbilinerecord[17]\tIRD\n";
		}
		elsif ($parsestrain[2] eq 'Quebec') {
			print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t50.736455\t-73.923341\t$ncbilinerecord[17]\tIRD\n";
		}
		elsif ($parsestrain[2] eq 'MB') {
			print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t53.981935\t-98.137208\t$ncbilinerecord[17]\tIRD\n";
		}
		elsif ($parsestrain[2] eq 'Prince_Edward_Island') {
			print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t46.346928\t-63.366852\t$ncbilinerecord[17]\tIRD\n";
		}
		elsif ($parsestrain[2] eq 'Nunavet') {
			print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t65.622023\t-97.302247\t$ncbilinerecord[17]\tIRD\n";
		}
		elsif ($parsestrain[2] eq 'Interior_Alaska') {
			print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t46.042736\t-94.573975\t$ncbilinerecord[17]\tIRD\n";
		}
		elsif ($parsestrain[2] eq 'Guatemala') {
			print $OUTPUT3 "$ncbilinerecord[0]\t$ncbilinerecord[1]\t$ncbilinerecord[2]\t$ncbilinerecord[3]\t$ncbilinerecord[4]\t$ncbilinerecord[5]\t$ncbilinerecord[6]\t$ncbilinerecord[7]\t$ncbilinerecord[8]\t$ncbilinerecord[9]\t$ncbilinerecord[10]\t$ncbilinerecord[11]\t$ncbilinerecord[12]\t$ncbilinerecord[13]\t$ncbilinerecord[14]\t15.538376\t-90.069581\t$ncbilinerecord[17]\tIRD\n";
		}
	}
	else {
		print $OUTPUT3 "$_\n";
	}
}