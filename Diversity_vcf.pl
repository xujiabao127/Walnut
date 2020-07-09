use strict;
use File::Basename;
use Getopt::Long;

############################## Caculate the theta pi tajiamd Fst

my ($vcf,$pop1,$pop2,$output,$window,$step,$snpnum);
GetOptions(
		"vcf=s" => \$vcf,
		"pop1=s"   => \$pop1,
		"pop2=s"   => \$pop2,
		"output=s"   => \$output,
		"window=i" =>\$window,
		"step:i" =>\$step,
		"snpnum:i" =>\$snpnum
		);
################  perl Diversity_vcf.pl -vcf Chr31.snp_result.vcf.gz -pop1 pop1.txt  -pop2 pop2.txt -window 10000 -step 1000 -snpnum 1 -output  Diversity_Fst.window.xls

$window||=10000;
$step||=1000;
$snpnum||=1;

print "$pop1 start\n";
open OUT,">$output" or die $!;
print OUT "Chr\tPos\tSamNum1\tSNP1\tTheta_w1\tTheta_p1\tTajima_D1\tSamNum2\tSNP2\tTheta_w2\tTheta_p2\tTajima_D2\tFst_Hundson\n";
my($samnum1,$samnum2);
$samnum1=$samnum2=0;
my %pop1;
my %pop2;
open POP1,"$pop1" or die $!;
while (<POP1>) {
	chomp;
	next if($_!~/\w/);
	my @a=split/\s+/,$_;
	my $sam=$a[0];
	$pop1{$sam}=0;
	$samnum1++;
}
close POP1;
print "$pop2  start\n";
open POP2,"$pop2" or die $!;
while (<POP2>) {
	chomp;
	next if($_!~/\w/);
	my @a=split/\s+/,$_;
	my $sam=$a[0];
	$pop2{$sam}=0;
	$samnum2++;
}
close POP2;
print "$samnum1\t$samnum2\n";

print "$vcf start\n";
my (%snpdb1,%snpdb2,%bin2loc);
my $line=0;
open VCF,$vcf=~/gz$/?"zcat $vcf|":"<$vcf" or die $!;
my @sams;
while (<VCF>) {
	chomp;
	if($_=~/^#/){
		if($_=~/\#CHROM/){
			my @tmp=split/\s+/,$_,10;
			@sams=split/\s+/,$tmp[9];
		}
		next;
	}
	$line++;
	my @a=split/\s+/,$_,10;
	my $chr=$a[0];
	my $pos=$a[1];
	my $ref=$a[3];
	my $alt=$a[4];

	my @nts=split/\s+/,$a[9];
	my %sam2nt;
	for (my $i=0 ;$i<@nts ;$i++) {
		my $sam=$sams[$i];
		my $nts=(split/\:/,$nts[$i])[0];
		$sam2nt{$sam}=$nts;
	}
	my (@nts1,@nts2);
	foreach my $sam (sort keys %pop1) {
		push @nts1,$sam2nt{$sam};
                if(!exists $sam2nt{$sam}){
                        print "1\t$chr\t$pos\t$sam\n";
                }
	}
	foreach my $sam (sort keys %pop2) {
		if(!exists $sam2nt{$sam}){
			print "2\t$chr\t$pos\t$sam\n";
		}
		push @nts2,$sam2nt{$sam};
	}
#	print "1\t$chr\t$pos\t@nts1\n";
#	print "2\t$chr\t$pos\t@nts2\n";
	my $nts1=join " ",@nts1;
	my $nts2=join " ",@nts2;
	$snpdb1{$chr}{$pos}=$nts1;
	$snpdb2{$chr}{$pos}=$nts2;
	my $site=int($pos/$window)*$window;
	# print "$chr\t$pos\t$nts1\tT\t$nts2\n";
	for (my $loc=($site-$window);$loc<($site+$window);$loc+=$step) {
		my $sts=$loc;
		my $ens=$loc+$window;
		if($sts<0){$sts=0};
		if($pos>=$sts && $pos<=$ens){
			$bin2loc{$chr}{$sts}{$pos}=0;
			#print "$chr\t$sts\t$ens\t$pos\n";
		}
	}
	if($line%10000==0){
		print "$line\t$chr\t$pos\n";
	};
	#if($pos>100000){last};
}
close VCF;
print "Diversity for sliding window start\n";
foreach my $chr (sort keys %bin2loc) {
	foreach my $bin (sort {$a<=>$b} keys %{$bin2loc{$chr}}) {
		my %divs;
		my $between_pi=0;
		my $within_pi=0;
		my $snpnum_real=scalar(keys %{$bin2loc{$chr}{$bin}});
		if($snpnum_real<$snpnum){next};
		$divs{pop1}{snp}=$divs{pop1}{theta_w}=$divs{pop1}{theta_p}=0;
		$divs{pop2}{snp}=$divs{pop2}{theta_w}=$divs{pop2}{theta_p}=0;
		foreach my $loc (sort {$a<=>$b} keys %{$bin2loc{$chr}{$bin}}) {
			#print "$chr\t$bin\t$loc\t$snpdb1{$chr}{$loc}\tT\t$snpdb2{$chr}{$loc}\t";
			my @nts1=split/\s+/,$snpdb1{$chr}{$loc};
			my @nts2=split/\s+/,$snpdb2{$chr}{$loc};
			my @nts=split/\s+/,$snpdb1{$chr}{$loc}."\t".$snpdb2{$chr}{$loc};
			###################### all 
			my (%baseall,%base1,%base2);
			foreach my $nt (@nts) {
				if($nt ne "./."){
					my @nt=split/\//,$nt;
					$baseall{$nt[0]}++;
					$baseall{$nt[1]}++;
				}
			}
			my @bases=sort {$baseall{$b}<=>$baseall{$a}} keys %baseall;
			next if(scalar(@bases)!=2);
			my $ma_nt=$bases[0];
			my $mi_nt=$bases[1];
			###################### pop1  pop2 
			foreach my $nt (@nts1) {
				if($nt ne "./."){
					my @nt=split/\//,$nt;
					$base1{$nt[0]}++;
					$base1{$nt[1]}++;
				}
			}
			foreach my $nt (@nts2) {
				if($nt ne "./."){
					my @nt=split/\//,$nt;
					$base2{$nt[0]}++;
					$base2{$nt[1]}++;
				}
			}
			################## 
			foreach my $base (@bases) {
				if(!exists $base1{$base}){$base1{$base}=0};
				if(!exists $base2{$base}){$base2{$base}=0};
			}
			if(($base1{$ma_nt}+$base1{$mi_nt})>0){
				my $seg1=&snp($base1{$ma_nt},$base1{$mi_nt});
				my $theta_w1=&theta_w($base1{$ma_nt},$base1{$mi_nt});
				my $theta_p1=&theta_p($base1{$ma_nt},$base1{$mi_nt});
				$divs{pop1}{snp}+=$seg1;
				$divs{pop1}{theta_w}+=$theta_w1;
				$divs{pop1}{theta_p}+=$theta_p1;
			}
			if(($base2{$ma_nt}+$base2{$mi_nt})>0){
				my $seg2=&snp($base2{$ma_nt},$base2{$mi_nt});
				my $theta_w2=&theta_w($base2{$ma_nt},$base2{$mi_nt});
				my $theta_p2=&theta_p($base2{$ma_nt},$base2{$mi_nt});
				$divs{pop2}{snp}+=$seg2;
				$divs{pop2}{theta_w}+=$theta_w2;
				$divs{pop2}{theta_p}+=$theta_p2;
			}
			if(($base1{$ma_nt}+$base1{$mi_nt})==0 or ($base2{$ma_nt}+$base2{$mi_nt})==0){
				#print "$chr\t$loc\t@nts1\t\t\t@nts2\n";
				next;
			}		
			my $bp=&between_pi($base1{$ma_nt},$base1{$mi_nt},$base2{$ma_nt},$base2{$mi_nt});
			$between_pi+=$bp;
			# print "$seg1\t$theta_w1\t$theta_p1\t$seg2\t$theta_w2\t$theta_p2\t$bp\t$base1{$ma_nt},$base1{$mi_nt},$base2{$ma_nt},$base2{$mi_nt}\n";
		}
		next if($divs{pop1}{snp} == 0 and $divs{pop2}{snp}==0);
		if($between_pi==0){next};
		my $tajimaD1="NA";
		my $tajimaD2="NA";
		if($divs{pop1}{snp} >1){
			$tajimaD1=&tajimaD($divs{pop1}{snp},$divs{pop1}{theta_w},$divs{pop1}{theta_p},$samnum1);
		}
		if($divs{pop2}{snp} >1){
			$tajimaD2=&tajimaD($divs{pop2}{snp},$divs{pop2}{theta_w},$divs{pop2}{theta_p},$samnum2);
		}
		###  if snp number more than 1 ,TajimaD will be caculated
		# between_pi and within_pi the average number of pairwise differences between two individuals sampled from different sub-populations or from the same sub-population
		my $fst_hudson=1-($divs{pop1}{theta_p}+$divs{pop2}{theta_p})/2/$between_pi;
		print OUT "$chr\t$bin\t";
		print OUT "$samnum1\t$divs{pop1}{snp}\t$divs{pop1}{theta_w}\t$divs{pop1}{theta_p}\t$tajimaD1\t";
		print OUT "$samnum2\t$divs{pop2}{snp}\t$divs{pop2}{theta_w}\t$divs{pop2}{theta_p}\t$tajimaD2\t";
		print OUT "\t$fst_hudson\n";
	}
}
close OUT;



sub snp{
	my ($freq1,$freq2)=@_;
	if($freq1==0 or $freq2==0){
		return 0;
	}else{
		return 1;
	}
}
sub theta_w{
	my ($freq1,$freq2)=@_;
	if($freq1==0 or $freq2==0){
		return 0;
	}else{
		my $num=$freq1+$freq2;
		my $a1=0;
		for (my $si=1;$si<$num ;$si++) {
			$a1+=1/$si;
		}
		my $val=1/$a1;
		return $val;
	}
}
sub theta_p{
	my ($freq1,$freq2)=@_;
	my $num=$freq1+$freq2;
	my $val=2*$freq1*$freq2/($num*($num-1));
	return $val;
}
sub tajimaD{
	my ($snp,$theta,$pi,$samnum)=@_;
	my $n=$samnum*2;
	my $a1=0;
	my $a2=0;
	for (my $i=1;$i<$n ;$i++) {
		$a1+=1/$i;
		$a2+=1/$i**2;
	}
	my $b1=($n+1)/(3*($n-1));
	my $b2=2*($n**2+$n+3)/(9*$n*($n-1));
	my $c1=$b1-1/$a1;
	my $c2=$b2-($n+2)/($a1*$n)+$a2/$a1**2;
	my $e1=$c1/$a1;
	my $e2=$c2/($a1**2+$a2);
	my $d=$pi-$theta;
	my $dstd=$e1*$snp+$e2*$snp*($snp-1);
	my $dval=$d/$dstd**0.5;
	return $dval;
}
sub between_pi{
	my ($freq1,$freq2,$freq3,$freq4)=@_;
	my $val1=$freq1*$freq4+$freq2*$freq3;
	my $num1=$freq1+$freq2;
	my $num2=$freq3+$freq4;
	my $val=$val1/($num1*$num2);
	return $val;
}























