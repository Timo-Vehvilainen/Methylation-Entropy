use warnings;
use strict;
open my $read_file,"../files/bam/SRR847314-SRR847317_merged_sorted_withRG_mdups_karyotypic_recal.bam" or die "read_file:$!\n";
open my $snp_file, "../files/SNPs_not_hypermethylated.txt"  or die "snp_file:$!\n";
open my $OUT, ">../files/SNPs_not_hypermethylated_context.fa" or die "OUT:$!\n";

my $reads = "";
my $wanted_strand = 'CT';

while (my $line = <$snp_file>) {
	
	chomp($line);
	my @line_data=split /\t/,$line;
	my $chr = $line_data[0];
	my $position = $line_data[1];
	my $context_start = $position - 5;
	my $context_end = $position + 5;
	my $region = $chr.':'.$context_start.'-'.$context_end;
	my $out = `../samtools-1.3.1/samtools view ../files/bam/SRR847314-SRR847317_merged_sorted_withRG_mdups_karyotypic.bam $region`;
	my @read_lines = split /\n/, $out;
	foreach my $read (@read_lines) {
		chomp($read);
		my @read_data = split /\t/,$read;
		my $read_start = $read_data[3];
		my $ncl_strand = substr $read_data[14], -2;
		if (($read_start <= $context_end - 90) or ($read_start > $context_start) or
			($ncl_strand ne $wanted_strand)) {
			next;
		} 
		my $ncl = $read_data[9];
		my $context = substr $ncl, ($context_start - $read_start), 11;
		print $OUT "$context\n";
	}
}
close $read_file;
close $OUT;
close $snp_file;
