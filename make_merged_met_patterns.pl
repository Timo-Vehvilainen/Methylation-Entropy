use warnings;
use lib '/home/timo/Desktop/code/';
open my $read_file,"../files/promoter_bam_reads.sam" or die "read_file:$!\n";
open my $OUT, ">../files/met_patterns_merged.txt" or die "OUT:$!\n";

my $promoter_reads = "";
my $promoter_name = "";
my %methylation_sites;
my $location;
my $met_data;
my $state;
my @strands= ("GA", "CT");

#First store all the reads of one promoter into a string
while (my $line = <$read_file>) {

	#When on a read line, add it to the long string	
	if($line=~/^SRR/){
		$promoter_reads = $promoter_reads.$line;
	}
	#when not on an empty or a read line, store the name of the promoter
	elsif ($line ne "\n") {
		$promoter_name = $line;
		#print $OUT ">$line\n";
	}
	#else (when on empty line), stop reading and do analysis on the obtained string
	else{

		#Analyze each promoter twice; once for each strand
		foreach my $strand_var (@strands) {
			
			print $OUT ">$promoter_name\n";
			#First go through the string one through to record all the methylation positions
			foreach (split("\n", $promoter_reads)){
				chomp($_);
	
				my @data = split("\t", $_);
				my $chr = $data[2];
				my $read_start = $data[3];
				my $strand = substr($data[14], 5);
				$met_data = substr($data[16], 5);
				if ($strand eq $strand_var){
					foreach (split('', $met_data)){
						if (($_ eq 'z') or ($_ eq 'Z')){
							if (not exists($methylation_sites{$read_start})){
								$methylation_sites{$read_start} = 1;
							}
						}
					$read_start += 1;
					}
				}
			}
			#Then go through them again, this time comparing each read to the 
			#map of all methylation sites and printing each methylation site result
			foreach (split("\n", $promoter_reads)){
				chomp($_);
	
				my @data = split("\t", $_);
				my $chr = $data[2];
				my $read_start = $data[3];
				my $strand = substr($data[14], 5);
				$met_data = substr($data[16], 5);
				if ($strand eq $strand_var){
					foreach $location (sort keys %methylation_sites) {
						if (($location - $read_start) < 90 and ($location - $read_start) >= 0){
							$state = substr($met_data, ($location - $read_start), 1);
							if ($state eq "z"){
								print $OUT "0";
							}			
							elsif ($state eq "Z"){
								print $OUT "1";
							}
							else{
								print $OUT "-";
							}		
						}
						else {
							print $OUT "-";
						}
					}
				print $OUT "\n";
				}	
			}
			%methylation_sites = ();
			print $OUT "\n";
		}
		$promoter_reads = "";
		$promoter_name = "";
	}


}
