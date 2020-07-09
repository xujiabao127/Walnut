use strict;
use File::Basename;
use Getopt::Long;

############################## Caculate the theta pi  
my $pop=shift;
my $chrlist=shift;
my $divdir=shift;
my $div_out=shift;

############################## Head
print "$pop start\n";
my %pop2sam;
open POP,"$pop" or die $!;
while (<POP>) {
        chomp;
        my @a=split/\s+/,$_;
        my $sam=$a[0];
        my $pop=$a[1];
        if($pop!~/JY/){
                $pop2sam{$pop}{$sam}=0;
        }
}
close POP;
my %chrlist;
open my $CL,"$chrlist" or die $!;
while (<$CL>) {
        chomp;
        my @a=split/\t+/,$_;
        my $chr=$a[0];
        my $len=$a[1];
        $chrlist{$chr}=$len;
}
close $CL;

my %divs;
my %thetas;
my %snps;
foreach my $chr (sort keys %chrlist) {
        my $div_file="$divdir/$chr.div.xls";
        my $num=0;
        print "$div_file\n";
        open DIV,"$div_file" or die $!;
        my $head=<DIV>;
        chomp $head;
        my @head=split/\s+/,$head,4;
        my @bins=split/\s+/,$head[3];
        print "$chr\t@bins\n";
        while (<DIV>) {
                chomp;
                $num++;
                my @a=split/\s+/,$_,4;
                my @val=split/\s+/,$a[3];
                for (my $i=0;$i<@val ;$i++) {
                        my $bin=$bins[$i];
                        my $val=$val[$i];
                        my ($snp,$pi)=(split/\|/,$val)[0,1];
                        $snps{$bin}+=$snp;
                        $divs{$bin}+=$pi;
                        if($snp==1){
                                my $freqs=(split/\|/,$val)[2];
                                my @freq=split/\:/,$freqs;
                                my $sum=$freq[0]+$freq[1];
                                my $a1=0;
                                for (my $j=1;$j<$sum ;$j++) {
                                        $a1+=(1/$j);
                                }
                                $thetas{$bin}+=(1/$a1);
                        }
                }
                if($num%100000==0){
                        print "$a[0]\t$a[1]\t$num\t".localtime()."\n";
                }
        }
        close DIV;
        print "$chr\tAllsnp\t$num\n";
}

open OUT,">$div_out" or die $!;
foreach my $bin (sort keys %snps) {
        my $snp=$snps{$bin};
        my $pi=$divs{$bin};
        my $theta=$thetas{$bin};
        print OUT "Sum_Div:\t$bin\t$snp\t$theta\t$pi\n";
}
close OUT;
