## mergy the TPMs of BGCs from multiple samples
#!usr/bin/perl -w
use strict;

my $out=shift;
my %TPM;
my @ls;
open IN,"sample_List.txt" or die; 
while(<IN>){chomp;push @ls,$_;}close IN;
foreach my $i(@ls){
#       print "read $i\n";
        open IN,"contig.bgc.result/$i.bgc.TPM.sum.txt" or die;
        while(<IN>){
          #     next if $_=~/^#/;
                chomp;
          #     s/^\s+//;   # delete the space in front
                my @arr=split;
                $TPM{$i}{$arr[0]}=$arr[1];
        }close IN;
}

open OUT,">$out";
foreach my $s(@ls){print OUT "\t$s";}print OUT "\n";
open IN,"all.uniq.BGCs" or die;
while(<IN>){
        chomp;
        print OUT "$_\t";
        foreach my $s(@ls){
                        if(exists $TPM{$s}{$_}){
                                print OUT "$TPM{$s}{$_}\t";
                                           }
                        else{print OUT "0\t";
                            }
                     }
        print OUT "\n";
}
close IN;
close OUT;
