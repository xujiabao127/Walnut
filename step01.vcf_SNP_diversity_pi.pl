
use strict;
use File::Basename;
use Getopt::Long;

############################## Caculate the theta pi tajiamd Fst

my ($vcf,$pop,$output);
GetOptions(
                "vcf=s" => \$vcf,
                "pop=s"   => \$pop,
                "output=s"   => \$output,
                );
################  perl Diversity_vcf.pl -vcf Chr01.snp_result.vcf.gz -pop sample2pop.tx -output  Diversity_SNP_pi.xls


print "$pop start\n";
my %pop2sam;
open POP,"$pop" or die $!;
while (<POP>) {
        chomp;
        next if($_!~/\w/);
        my @a=split/\s+/,$_;
        my $sam=$a[0];
        my $pop=$a[1];
        if($pop!~/JY/){
                $pop2sam{$pop}{$sam}=0;
        }
}
close POP;


open OUT,">$output" or die $!;
my @pops=sort keys %pop2sam;
my $popstr=join "\t",@pops;
print OUT "Chr\tPos\tRef\t$popstr\n";
my $line=0;
print "$vcf start\n";
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
        my %alles;
        for (my $i=0 ;$i<@nts ;$i++) {
                my $sam=$sams[$i];
                my $nts=(split/\:/,$nts[$i])[0];
                $sam2nt{$sam}=$nts;
                next if($nts eq "./.");
                my @alle=split/\//,$nts;
                $alles{$alle[0]}++;
                $alles{$alle[1]}++;
        }
        my @alles=sort {$alles{$b}<=>$alles{$a}}keys %alles;
        next if(scalar(@alles)!=2);
        my $mant=$alles[0];
        my $mint=$alles[1];
        my @popdiv;
        foreach my $pop (@pops) {
                my @sam=sort keys %{$pop2sam{$pop}};
                my %freq;
                foreach my $sam (@sam) {
                        next if($sam2nt{$sam} eq "./.");
                        my @base=split/\//,$sam2nt{$sam};
                        #print "$chr\t$pos\t$sam\t$pop\t$sam2nt{$sam}\t$base[0]\t$base[1]\n";
                        $freq{$base[0]}++;
                        $freq{$base[1]}++;
                }
                my @basenum=sort keys %freq;
                my $nums=scalar(@basenum);
                if (scalar(@basenum)!=2){
                        push @popdiv,"0|0|0:0";
                }else{
                        my $freq1=$freq{$mant};
                        my $freq2=$freq{$mint};
                        my $num=$freq1+$freq2;
                        my $pi=2*$freq1*$freq2/($num*($num-1));
                        push @popdiv,"1|$pi|$freq1:$freq2";
#                       print "2\t$pop\t$chr\t$pos\t@basenum\t$nums\t$pi\t$num\n";
                }
        }
        my $popdiv=join "\t",@popdiv;
        if($line%10000==0){
                print "$chr\t$pos\t$line\n";
        }
        print OUT "$chr\t$pos\t$ref\t$popdiv\n";
}
close VCF;
close OUT;
