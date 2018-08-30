open IN,"../files/SNP_dataset.txt" or die "IN:$!\n";
open OUT,">../files/P_Value.txt" or die "OUT:$!\n";
open OUT_hyper,">../files/SNPs_hypermethylated.txt" or die "OUT:$!\n";
open OUT_not_hyper,">../files/SNPs_not_hypermethylated.txt" or die "OUT:$!\n";


print OUT "CHR\tPOS\tNAME\tREF\tALT\tP-VALUE\n";
while($line=<IN>)
{
	chomp $line;
	@read=split/\t/,$line;
	if($line!~/^>chr/ && $line!~/^SNP/)
	{
		if($read[0] eq $snp1)
		{
			$para[3]=$read[1];
			$para[1]=$read[3];
			$count++;
		}
		elsif($read[0] eq Total2)
		{
			$para[2]=$read[1];
			$para[0]=$read[3];
			$count++;
		}
	}
	
	elsif($line=~/^>chr/)
	{
		if($count==2)
		{
			$p_val = &p_value_calculation($para[0], $para[1], $para[2], $para[3]);
			if($p_val<0.05&&$p_val>0)
			{
				$null++;
			}
			else
			{
				$alternative++;
			}
				
			print OUT "$p_val\n";
			if($p_val<0.05&&$p_val>0)
			{
				if ($para[3] > ($para[2] - $para[3])) {
					print OUT_hyper "$snp_name\t$p_val\n";
				}
				elsif ($para[3] < ($para[2] - $para[3])) {
					print OUT_not_hyper "$snp_name\t$p_val\n";
				}
			}
#			print OUT "$para[0]\t$para[1]\t$para[2]\t$para[3]\n";
			$count=0;
		}
		@snp=split/\t/,$line;
		$snp1=@snp[3];
		$snp2=@snp[4];
		print OUT "$line\t";
		$snp_name = $line;
		
		$snpCount++;	
	}
}
print OUT "$p_val\n";
print OUT "\n\n$snpCount\t$null\t$alternative\n";
sub p_value_calculation{
	my($N, $n, $M, $x) = @_;
	$c1 = 0;
	$c2 = 0;
	$c3 = 0;
	$p_value = 0;
	if($n < $M){
		$min = $n;
	}
	else{
		$min = $M;
	}
	for($i = $x; $i <= $min; $i++){
		for($j=1; $j <= $i; $j++){
			$c1 = $c1 + log($M-$j+1) - log($j);
		}
		for($j=1; $j <= $n-$i; $j++){
			$c2 = $c2 + log($N-$M-$j+1) - log($j);
		}
		for($j=1; $j <= $n; $j++){
			$c3 = $c3 + log($N-$j+1) - log($j);
		}
		$p_value += exp($c1+$c2-$c3);
		$c1 = 0, $c2 = 0, $c3 = 0;
	}
	return $p_value;
}
