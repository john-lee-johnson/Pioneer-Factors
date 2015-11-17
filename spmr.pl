#!/usr/bin/perl -w

open(IN,"$ARGV[1]");
open(OUT,">$ARGV[2]");

while(<IN>) {
    if (/^\d+/){
        chomp;
        @_=split;
        printf OUT "%d\t%.4f\n",$_[0],$_[1]/$ARGV[0];
    }
    else{
        print OUT $_;
    }
}
