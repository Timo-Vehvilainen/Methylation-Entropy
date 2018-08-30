#package Methylation_entropy;
#use strict;
use warnings;

#use Exporter;

#our @ISA = qw(Exporter);

#our @EXPORT = qw( log10 log2 freg_table ME);

sub log10 {
	my $num = @_;
	print("num == $num\n");
	return log($num)/log(10);
}

sub log2{
	my $num = @_;
	return log($num)/log(2);
}

sub sum{
	my @ar = @_;
	my $sum = 0;
	for (@ar) {
		$sum += $_;
	}
	return $sum;
}

sub mean {
	my @ar = @_;
	my $len = @ar;
	my $mean = sum(@ar) / $len;
	return $mean;
}

sub var{
	my @ar = @_;
	my $len = @ar;
	my @vars;
	my $mean = mean(@ar);
	for (@ar) {
		push(@vars, (($_ - $mean)**2))
	}
	my $var_sum = sum(@vars);
	my $variance = $var_sum / $len;
	return $variance;
}

sub met_level {
	my @seqs = @_;
	my $sum = 0;
	my $count = 0;
	foreach my $val (@seqs) {
		my @chars = split //, $val;
		foreach my $char (@chars) {
				$count += 1;
				$sum += $char;		
		}
	}
	my $average = $sum / $count;
	return $average;
}

#Takes in same input as ME.
sub freq_table{
	my @seqs = @_;
	my %count;
	foreach my $val (@seqs) {
		$count{$val}++;
	}
	return values %count;
}

sub num_of_uniq{
	my @seqs = @_;
	my @unique_patterns = do{my %seen; grep {!$seen{$_}++} @seqs };
	my $len = @unique_patterns;
	return $len;
}

#Takes in an array of binary methylation strings, where 1 represents a methylated
#CpG site, and a 0 an unmethylated CpG site. For example:
# @seqs = ['0110', '1111', '0100', '0110', '0001', '1110', '1111'];
sub ME {
	my @seqs = @_;
	#my $e = log(10)/log(2);
	my $e = 1;	
	my $b = length($seqs[0]);
	my @n = freq_table(@seqs);
	my $N = @seqs;
	#print ("b = $b\ne = $e\nN = $N \nn = @n\n");
	
	my $sum = 0;
	foreach my $n_i (@n) { 
		$sum += ((((-1)* $n_i) / $N) * (log($n_i/$N)/log(2)));
	}
	return (($e/$b) * $sum);

}


open $IN,"../files/trimmed_met_blocks_CT_4.txt" or die "IN:$!\n";
open $OUT,">../files/methylation_entropies_CT_4.txt" or die "OUT:$!\n";


my @entropy_averages;
my @level_averages;
my @promoter_mets;
my @promoter_levs;

my $str = "";
my $promoter_name = "";
my $met_sum = 0;
my $met_count = 0;
my $met_level_sum = 0; 
my $promoter_level_average = "";
while (my $input = <$IN>){

	if ($input =~/^#/){
		$promoter_level_average = $input;
		next;
	}
	#if come across a promoter name
	if ($input =~/^>/){

		#Calculate and print the numbers from the previous block
		if ($promoter_name ne ""){
			my $met_average = $met_sum / $met_count;
			my $met_level_average = $met_level_sum / $met_count;
			my $variance = var(@promoter_mets);
			my $std = sqrt($variance);
			push(@entropy_averages, $met_average);
			push(@level_averages, $met_level_average);
			$met_sum = 0;
			$met_count = 0;
			$met_level_sum = 0;
			print $OUT "\nMethylation entropy mean = $met_average\n";
			print $OUT "Methylation entropy variance = $variance\n";
			print $OUT "Methylation entropy std = $std\n";
			print $OUT "Methylation level block average = $met_level_average\n";
			print $OUT "$promoter_level_average";
			@promoter_mets = ();
			@promoter_levs = ();
		}

		#Set the current line to be the next promoter
		$str = "";
		$promoter_name = $input;
		print $OUT "\n$promoter_name\n";
		
	}
	
	elsif ($input ne "" and $input ne "\n"){
		$str .= $input;
	}
	#parse data line
	if (($input eq "\n") and ($str ne "")){
		my @seq = split /\n/, $str;
		my $met_val = ME(@seq);
		push(@promoter_mets, $met_val);
		my $met_lev = met_level(@seq);
		push(@promoter_levs, $met_lev);
		#my $uniq_val = num_of_uniq(@seq);
		print $OUT "$met_val\t$met_lev\n"; 
		$met_sum += $met_val;
		$met_level_sum += $met_lev;
		$met_count += 1;
		$str = "";
	}
}

my $met_average = $met_sum / $met_count;
my $met_level_average = $met_level_sum / $met_count;
my $variance = var(@promoter_mets);
my $std = sqrt($variance);
push(@entropy_averages, $met_average);
push(@level_averages, $met_level_average);
$met_sum = 0;
$met_count = 0;
$met_level_sum = 0;
print $OUT "\nMethylation entropy mean = $met_average\n";
print $OUT "Methylation entropy variance = $variance\n";
print $OUT "Methylation entropy std = $std\n";
print $OUT "Methylation level block average = $met_level_average\n";
print $OUT "$promoter_level_average";

close $IN;
close $OUT;
