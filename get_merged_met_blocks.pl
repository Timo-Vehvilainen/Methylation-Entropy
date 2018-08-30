#use warnings;

sub get{
	($str, $y, $x) = @_;
	my @lines1 = split /\n/, $str;
	return(substr($lines1[$y], $x, 1));
}

sub get_blocks {
	my ($corners, $str) = @_;
	my $box = "";
	my @boxes = split /\n/, $corners;
	my $ret = "";
	my $i;
	my $j;
	foreach $box (@boxes){
		my @coords = split /\t/, $box;
		for $i (($coords[0]) ... ($coords[2])){
			for $j ($coords[1] ... $coords[3]){
				$ret .= get($str, $i, $j);
			}
			$ret .= "\n";
		} 
		$ret .= "\n";
	}
	return $ret;
}

sub get_corners {
	my ($ar) = @_;
	my @lines = split /\n/, $ar;
	my $cols = length($lines[1]);
	my $rows = @lines;
	my $dead_end = 0;
	my $corners = "";
	my $width = 4;
	my $min_height = 16;
	my $i = 0;
	my $j = 0;
	my $m = 0;
	my $n = 0;
	my $m2 = 0;
	my $bottom = 0;
	
	foreach $i (0...($rows-1)) {
	#print("Found a total of $rows rows and $cols cols\n");
		foreach $j (0...($cols-1)){
			#print("checking for a box starting from ($i, $j)\n");
			if ((get($ar, $i, $j)) ne "-"){
				#print("\tposition ($i, $j) is NOT a dash, starting to count laterally\n");
				#first check if there is enough data laterally
				$m = $j;
				while (($m < ($cols)) and ((get($ar, $i, $m)) ne '-') and (($m - $j) <= $width)){
					$m = $m + 1;
				} 
				$m = $m - 1;
				#print("\t\tFound ", ($m - $j), " numbers laterally, from column $j to column $m\n");
				#if enough latitudal data, check longitudal for each column
				if (($m - $j) == ($width-1)){
					#print("\t\t\tStarting to count horizontally from row $i\n");
					$bottom = $rows;
					foreach $m2 ($j ... $m){ 
						#print("\t\t\tColumn $m2:\n");
						$n = $i;
						while (($n < $bottom) and (get($ar, $n, $m2)) ne '-'){	
							$n = $n + 1;
						}
						$n = $n - 1;
						#print("\t\t\t\tFound ", ($n-$i)," numbers horizontally, from row $i to row $n, m2 = $m2, j = $j\n");
						
						if (($n - $i) < $min_height) {
							#print("\t\t\t\t\tDead end, setting up flag to move to next box start position\n");
							$dead_end = 1;
							last;
						}
						elsif ($m2 eq $j) {
							#print("\t\t\t\t\tThis is the first checked column for this box, setting bottom from $bottom to ",($n+1),"\n");
							$bottom = $n+1;
						}
						elsif ($m2 eq $m) {
							#print("\t\t\t\t\tThis is the last checked column for this box, saving corners ($i, $j), ($n, $m)\n");
							#print("\t\t\t\t\tchecking if a box with a shared corner already exists....");
							if (not(($corners =~/($i, $j)/) or ($corners =~/($n, $m)/))){
								#print("\t\t\t\t\tNo duplicates found, saving box ($i, $j), ($n, $m)!\n");
								$corners.= "$i\t$j\t$n\t$m\n";
							} 
							else {
								#print("\t\t\t\t\tFound a duplicate already in the boxes. Not saving box ($i, $j), ($n, $m)\n");
							}
						}
					}
				}
				if ($dead_end eq 1){
					#print("\t\t\t\t\tDetected flag, going to next position\n");
					$dead_end = 0;
					next;
				}
			}
		}		
	}
	#print("Done!\n");
	return $corners;
}

open $IN,"../files/met_patterns_CT.txt" or die "IN:$!\n";
open $OUT,">../files/trimmed_met_blocks_CT_4.txt" or die "OUT:$!\n";
#open $OUT_lev,">../files/promoter_met_levels_merged_4.txt" or die "OUT:$!\n";

my $str = "";
my $promoter_name = "";
my $chomped_promoter = "";
my $met_level_average = 0;
my $ones = 0;
my $zeros = 0;
my $GA_blocks = "";
while (my $input = <$IN>){
	if ($input =~/^>/){
		if ($str ne ""){
			chomp($str);
			my $rev = reverse($str);
			chomp($rev);
			$str = reverse($rev);
			$ones += ($str =~ tr/1//);
			$zeros += ($str =~ tr/0//);
			#print $OUT_lev "$promoter_name\n\nmethylated = $ones\nunmethylated = $zeros\nmethylation level = $met_level_average\n\n";
			
			my $corners = get_corners($str);
			my $blocks = get_blocks($corners, $str);
			if ($corners ne ""){
				if ($input eq $promoter_name) {
					$GA_blocks = $blocks;;
				}
				else {
					if (($ones + $zeros) != 0) {
						$met_level_average = ($ones / ($ones + $zeros));
					}
					else {
						$met_level_average = 0;
					}
					
					print $OUT "$promoter_name\n$GA_blocks\n$blocks\n#Promoter methylation average = $met_level_average\n\n";
					$ones = 0;
					$zeros = 0;
				}
			}
			print("Done with promoter $promoter_name\n");
			$str = "";		
		}
		$promoter_name = $input;
	}
	elsif ($input ne "" and $input ne "\n"){
		$str.=$input;
	
	}
}

print $OUT "$blocks\n";

close $IN;
close $OUT;
#close $OUT_lev
