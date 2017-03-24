# This file is part of ONCOCNV - a tool to detect copy number alterations in amplicon sequencing data
# Copyright (C) 2016 OncoDNA

# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Author: Valentina BOEVA

#!/usr/bin/perl -w

=pod

=head1 NAME

ONCOCNV_getCounts.pl - Script to create a file with GC content for amplicon regions or exons

=head1 SYNOPSIS

perl createTargetGC.pl -bed <bed file with coordinates> -fi <fasta file with the genome> -od <tmp dir> -of <output file> 
     
=head1 DESCRIPTION

This is a command-line interface to createTargetGC.pl.


=head1 AUTHORS

Valentina Boeva <lt>valentina.boeva@inserm.fr<gt>

=cut

# -------------------------------------------------------------------

use POSIX;
use strict;
use warnings;
use diagnostics;


my $usage = qq{    

perl createTargetGC.pl -bed <bed file with coordinates> -fi <fasta file with the genome> -od <tmp dir> -of <output file> 
			
};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

## mandatory arguments

my $bed = "";
my $fastaGenome = "";
my $outputFname = "";
my $tmpDir = "";

## parse command line arguments
while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }	
    elsif ( $this_arg eq '-bed') {$bed = shift @ARGV;}
    elsif ( $this_arg eq '-fi') {$fastaGenome = shift @ARGV;}
    elsif ( $this_arg eq '-of') {$outputFname = shift @ARGV;}
    elsif ( $this_arg eq '-od') {$tmpDir = shift @ARGV;}

    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

if ( $fastaGenome eq "" ){
    die "you should specify .fasta file with the reference genome\n";
}
if( $tmpDir eq ""){
    $tmpDir="./";
}

if( $bed eq ""){
    die "you should specify .bed file with coordinates of the target regions\n";
}

if ( $outputFname eq ""){
    die "you should specify output filename\n";
}

my $bedToolsCommand="bedtools getfasta -bed $bed -fi $fastaGenome -fo $tmpDir/target.fasta 2>$tmpDir/tmpFile.txt; rm $tmpDir/tmpFile.txt";
system("$bedToolsCommand");

#check whether the file is empty because of "chr" not present in the genome:

open FILE, "< $tmpDir/target.fasta " || die "$tmpDir/target.fasta : $!\n";
my $string = <FILE>;
close FILE;

if (!defined($string) || length($string) ==0) {
    	print "..Oops.. File $tmpDir/target.fasta is empty!\n";
	print "..It seems that there is not 'chr' prefixes in your reference genome fasta file..\n..But no worries! OncoCNV will adjust for it\n";
	open FILE, "< $bed " || die "$bed : $!\n";
	open OUT, "> $tmpDir/noChr.bed " || die "$tmpDir/noChr.bed : $!\n";
	while(<FILE>){
	    print OUT substr $_, 3;
	}
	close FILE;
	close OUT;
	
	my $bedToolsCommand="bedtools getfasta -bed $tmpDir/noChr.bed -fi $fastaGenome -fo $tmpDir/target.fasta; rm $tmpDir/noChr.bed";
	system("$bedToolsCommand");
}

open OUT, "> $outputFname " || die "$outputFname : $!\n";
open FILE, "< $tmpDir/target.fasta " || die "$tmpDir/target.fasta : $!\n";
print OUT "region\tGC\n";
						
my ($name,$GC);
while(<FILE>){
    chomp;
    if (/>(.+)/) {
	$name=$1;
	$name="chr".$1 unless ($1 =~/chr/);
    } else {
	my $len = length($_);
	$GC=(tr/GCgc//) / ($len-(tr/nN//));
	print OUT "$name\t$GC\n";
    }
}
close FILE;
close OUT;

my $command="rm $tmpDir/target.fasta";
system("$command");

