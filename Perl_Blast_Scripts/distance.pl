###########################################################################################################################################
## This is probably one of my most used scripts. It is very generic. All you need is the lat and long of some starting object in the first two columns of an excel file. The 
## lat long of the second object is followed in columns three and four. Save as tab delimited and change the input and output files and you are ready to go.
###########################################################################################################################################

##REMEMBER TO ELIMINATE HEADERS FROM THE FILE AND MAKE SURE THAT THERE ARE NO UNPOPULATED CELLS!!

use strict;
use warnings;
use GIS::Distance;

open (COORD, "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_coords_2nt.txt") || die "cannot open file";

print "COORD opened successfully\n\n"; # Just to make sure you entered into the program successfully

my $gis = GIS::Distance->new(); #Here you are defining the default GIS:: Distance object just editing for the Haversine formula
$gis->formula('Haversine');

my $OUTPUT_FILE;
my @distances;
my $distance;

while (<COORD>) { #Enter into your tab file that has lat1 lon1 lat2 long2
	$_ =~ m/(^\S\d+\.\d+)\s+(\S\d+\.\d+)\s+(\S\d+\.\d+)\s+(\S\d+\.\d+)/; # Capture the variable for each column
	$distance = $gis->distance ($1, $2 => $3, $4); # This is part of the GIS::Distance object... with do the calculations within the object
	push @distances, $distance->kilometers(); # just tell it the type of distance you want to output
}

open ($OUTPUT_FILE, ">" . "Tony/Influenza_Genomics/OSU_Surveillance_Analysis/IRD_Blast/MP/MP_coords_distances_2nt.txt" || die "cannot open");

foreach my $dist (@distances) {
	print $OUTPUT_FILE "$dist\n";	# This loop just takes you through and allows you to put the distance on each line. It is important that the number of lines matches your original length of input
}
close COORD;