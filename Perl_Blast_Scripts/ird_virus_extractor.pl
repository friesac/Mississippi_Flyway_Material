#########################################################################################################################
#########################################################################################################################
# This takes two excel sheets of all sequence names and 8 columns and the same 8 columns transposed and creates different tables that 
# parse out all the information from the names. MUST BE IN IRD FORMAT WHEN INPUT. Used for the MS Flyway Supplementatal table!
#########################################################################################################################
#########################################################################################################################
use strict;
use warnings;
#########################################################################################################################
my $OUTPUT1; 
my $OUTPUT2;
my $OUTPUT3;
my $OUTPUT4;
my $line;
my $pb2;
my $pb1;
my $pa;
my $ha;
my $np;
my $na;
my $mp;
my $ns;
my $isolate;
my $gene;
my $num;
my $match;
my $virus;
my @record;
my @list;
my @ird_pb2;
my @ird_pb1;
my @ird_pa;
my @ird_ha;
my @ird_np;
my @ird_na;
my @ird_mp;
my @ird_ns;
my @tmp;
my @AoA;
my @metadata;
my %myhash;
my %table;
#########################################################################################################################
open (VIRUSESTITLES, "IRD_METADATA/ird_virus_titles.txt")||die "cannot open file"; 
open ($OUTPUT1, ">" . "IRD_METADATA/ird_virus_segment_counts.txt")||die "cannot open file"; # File with key and segment and number of occurences... should be around 8 of each isolate.

while (<VIRUSESTITLES>){

	chomp;
	$line = $_;
	@record = split(/\t/, $_); 			# parse the line into an array by tab delimiter	
	$pb2 = $record[0];					# assign the first element to a dummy scalar for the segment column
	@ird_pb2 = split(/\|/, $pb2);		# split that dummy scalar based on the "|" delimiter... have all the metadata in array now.
	$myhash {$ird_pb2[1]}{'pb2'}++;		# assign the second element of that array to the key of the hash with the number of segments.
	$pb1 = $record[1];					# repeat the rest of the way down.
	@ird_pb1 = split(/\|/, $pb1);
	$myhash {$ird_pb1[1]}{'pb1'}++;
	$pa = $record[2];
	@ird_pa = split(/\|/, $pa);
	$myhash {$ird_pa[1]}{'pa'}++;
	$ha = $record[3];
	@ird_ha = split(/\|/, $ha);
	$myhash {$ird_ha[1]}{'ha'}++;
	$np = $record[4];
	@ird_np = split(/\|/, $np);
	$myhash {$ird_np[1]}{'np'}++;
	$na = $record[5];
	@ird_na = split(/\|/, $na);
	$myhash {$ird_na[1]}{'na'}++;
	$mp = $record[6];
	@ird_mp = split(/\|/, $mp);
	$myhash {$ird_mp[1]}{'mp'}++;
	$ns = $record[7];
	@ird_ns = split(/\|/, $ns);
	$myhash {$ird_ns[1]}{'ns'}++;	
}
my $hashdimensioncheck = scalar keys %myhash;				# check for the number of keys in the hash by turning it scalar. number 
print "The number of keys is: $hashdimensioncheck";

while ( my ($key, $value) = each (%myhash)){ 				# now go through that hash and get each key and value
	while (my  ($segment, $count) = each (%$value)){		# now get the count associated witht the value... this is a basic subroutine i've	
		print $OUTPUT1 "$key\t$segment\t$count\n";			# stolen from a script i've seen to count occurences in a hash.
	}
}
close $OUTPUT1;
#########################################################################################################################################
open (TABLE, "IRD_METADATA/ird_virus_segment_counts.txt")||die "cannot open file";
open ($OUTPUT2, ">" . "IRD_METADATA/ird_viruses.txt")||die "cannot open file";

while (<TABLE>) {					# take that same file you just made and go through line by line	
	chomp;
	@list= split(/\t/, $_); 		# parse the first line into an array by tab delimiter
	$isolate = $list[0];			# grab the key
	$gene = $list[1];				# grab the value
	$num = $list[2];				# grab the count of the key->value pair	
	$table{$isolate} = [] unless exists $table{$isolate};	# put the isolate key from the previous hash into a new hash of arrays unless it exists
	push @{$table{$isolate}}, "$gene\t";	# push all of the gene values from the previous hash into the array value of the new hash of arrays.
}
	
foreach $isolate (sort keys %table) { 	# assign the keys to a new dummy variable
	print $OUTPUT2 "$isolate\t";		# print the keys	
	my @gene = @{$table{$isolate}};		# assign the dereferenced array value from the hash of arrays to an actual array for printing
	print $OUTPUT2 sort "@gene\t$num\n";		# send all of the segments to the line with the key and print
}
close TABLE;
close $OUTPUT2;
#########################################################################################################################################
open (VIRUSARRAY, "IRD_METADATA/ird_virus_titles_transposed.txt")||die "cannot open file";
open ($OUTPUT3, ">", "IRD_METADATA/ird_virus_table.txt")||die "cannot open file";				# creates the table of accession numbers
open ($OUTPUT4, ">", "IRD_METADATA/ird_virus_table_meta.txt")||die "cannot open file";			# creates the table of metadata to add to accession table

while (<VIRUSARRAY>) { 								# enter into a transposed file with sequence titles in a row for each segment (i.e. 8 by ~1600)
	@tmp = split(/\t/, $_);							# assign the first row to a temp array that contains \t delimited virus names for the segment	
	push @AoA,[@tmp];								# push that array to an array... i.e. an array of arrays	
}

foreach my $key (sort keys %table) {				# take all of the unique virus titles regardless of how many segments each has
	print $OUTPUT3 "$key\t";						# print that key to a final table file	
	foreach my $segment (@AoA) {					# enter into the array of arrays... $segment is a referenced array!!
		foreach my $virus (@$segment) {				# now take each element of the dereferenced array which will be a virus title from IRD
			if ($virus =~ m/($key)/) {				# for each of the IRD titles see if you match the current key
				$match = $1;						# if you find the key in the IRD title then capture that key... kinda of redundant but needed for test a few lines down
				@metadata = split(/\|/, $virus);	# take the IRD title and split on the "|"	
				print $OUTPUT3 "$metadata[0], ";	# print the accession number of the IRD title to the line with the key(i.e. isolate name) and comma in case there are more than one
			}
		}											# grab the next isolate and start the matching all over	
		if (defined $match) {						# if the $match variable had been defined then you will print a "\t" after the accession number 
			print $OUTPUT3 "\t";					# before moving onto the next virus segment (i.e. pb1)
		}												
		else {										# if you never found a match in that gene segment array		
			print $OUTPUT3 "X\t";					# print "X\t" before moving on.	
		}											
		$match =~ /(\A)/;							# clear the capture variable	
	}												# move onto the nec segment (ie 2 of 8 then 3 of 8... each representing a gene)	
	print $OUTPUT3 "\n";							# print "\n" before going onto the next unique key.
}													# go to the next key

foreach my $key (sort keys %table) {				# go back into the above loop but grap additional metadata.
	print $OUTPUT4 "$key\t";	
	foreach my $segment (@AoA) {					# enter into the array of arrays... $segment is a referenced array!!
		foreach my $virus (@$segment) {				# now take each element of the dereferenced array which will be a virus title from IRD
			if ($virus =~ m/($key)/) {				# for each of the IRD titles see if you match the current key
				$match = $1;						# if you find the key in the IRD title then capture that key... kinda of redundant but needed for test a few lines down
				@metadata = split(/\|/, $virus);	# take the IRD title and split on the "|"	
				print $OUTPUT4 "$metadata[2]\t$metadata[3]\t$metadata[4]\t$metadata[5]\t$metadata[6]\t$metadata[7]";	# print the accession number of the IRD title to the line with the key(i.e. isolate name) and comma in case there are more than one
			}
		}											# grab the next isolate and start the matching all over	
		if (defined $match) {						# if the $match variable had been defined then you will print a "\t" after the accession number 
			print $OUTPUT4 "\t\t\t\t\t";					# before moving onto the next virus segment (i.e. pb1)
		}												
		else {										# if you never found a match in that gene segment array		
			print $OUTPUT4 "X\t";					# print "X\t" before moving on.	
		}											
		$match =~ /(\A)/;							# clear the capture variable	
	}												# move onto the nec segment (ie 2 of 8 then 3 of 8... each representing a gene)	
	print $OUTPUT4 "\n";							# print "\n" before going onto the next unique key.
}	
	