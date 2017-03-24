# This file is part of ONCOCNV - a tool to detect copy number alterations in amplicon sequencing data
# Copyright (C) 2014 OncoDNA

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

ONCOCNV_getCounts.pl - Script to count reads in targeted regions from mapped sequencing data (.BAM)

=head1 SYNOPSIS

ONCOCNV_getCounts.pl getControlStats -m <mode "Ampli" or "Exon"> -c <control_1.bam,control_2.bam,...> -b <target.bed> -o <output file> 

or

ONCOCNV_getCounts.pl getSampleStats -m <mode "Ampli" or "Exon"> -c <control.txt> -s <sample_1.bam,sample_2.bam,...> -o <output file> 

      
    
=head1 DESCRIPTION

This is a command-line interface to ONCOCNV_getCounts.pl.


=head1 AUTHORS

Valentina Boeva <lt>valentina.boeva@inserm.fr<gt>

=cut

# -------------------------------------------------------------------

use POSIX;
use strict;
use warnings;
use diagnostics;


my $usage = qq{    

ONCOCNV_getCounts.pl getControlStats -m <mode "Ampli" or "Exon"> -c <control_1.bam,control_2.bam,...> -b <target.bed> -o <output file> 

or

ONCOCNV_getCounts.pl getSampleStats -m <mode "Ampli" or "Exon"> -c <control.txt> -s <sample_1.bam,sample_2.bam,...> -o <output file>

			
};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

## mandatory arguments

my $sample = "";
my $control = "";
my $output_fname = "";
my $bed = "";
my $base = "Exon";

## parse command line arguments

my $command = shift @ARGV; 

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-c') {$control = shift @ARGV;}
    elsif ( $this_arg eq '-s') {$sample = shift @ARGV;}
    elsif ( $this_arg eq '-o') {$output_fname = shift @ARGV;}
    elsif ( $this_arg eq '-b') {$bed = shift @ARGV;}
    elsif ( $this_arg eq '-m') {$base = shift @ARGV;}
    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

if ( $control eq "" ){
    die "you should specify at least one control .BAM file (getControlStats) or a processed control file (getSampleStats)\n";
}
if( $sample eq "" && $command eq "getSampleStats"){
    die "you should specify at least one sample filename\n";
}

if( $bed eq "" && $command eq "getControlStats"){
    die "you should specify .bed file with coordinates of the target regions\n";
}

if ( $output_fname eq ""){
    die "you should specify output filename\n";
}


my %bedcoord;
####################
#read .bed provided:
####################

my $numberOfTargetedRegions = 0;

if( $bed ne "" && $command eq "getControlStats"){

	open FILE, "< $bed " || die "$bed : $!\n";
	while(<FILE>){
    		chomp;
		if (/(chr\S+)\s(\d+)\s(\d+)\s([\w\-\.]+)\s\S+\s([\w\-\.]+)/) {
			$bedcoord{$1}->{$4}{start}=$2;
			$bedcoord{$1}->{$4}{end}=$3;
			$bedcoord{$1}->{$4}{gene}=$5;
			$numberOfTargetedRegions++;
		} elsif (/^(\S+)\s(\d+)\s(\d+)\s([a-zA-Z0-9_\-\.]+)\s\S+\s([a-zA-Z0-9_\-\.]+)/) {
			my $chr=$1;
			$chr="chr".$1 unless ($1 =~/chr/);			
			$bedcoord{$chr}->{$4}{start}=$2;
			$bedcoord{$chr}->{$4}{end}=$3;
			$bedcoord{$chr}->{$4}{gene}=$5;
			$numberOfTargetedRegions++;
		}
	}
	close FILE;
	print ("number of targeted regions: ",$numberOfTargetedRegions, "\n");
	#check region boundaries and merge overlapping regions if exon-based:
	if ($base eq "Exon") {
	   if ($numberOfTargetedRegions<50000) {
	 	for my $chr (keys (%bedcoord)) {
STARTCHR:
		  my @regions = sort keys (%{$bedcoord{$chr}});
		  for my $myRegion1 (@regions) {
			for my $myRegion2 (@regions) {
				if ($myRegion1 ne $myRegion2 && $bedcoord{$chr}->{$myRegion2}{end}>=$bedcoord{$chr}->{$myRegion1}{start} && $bedcoord{$chr}->{$myRegion2}{start}<=$bedcoord{$chr}->{$myRegion1}{end}) 				{		
					print "These two regions overlap: $chr $bedcoord{$chr}->{$myRegion1}{start} - $bedcoord{$chr}->{$myRegion1}{end} and $bedcoord{$chr}->{$myRegion2}{start} - $bedcoord{$chr}->{$myRegion2}{end}\n"; 
					my $regionName=$myRegion1.".".$myRegion2;
					my $gene = $bedcoord{$chr}->{$myRegion1}{gene};
					if ($bedcoord{$chr}->{$myRegion1}{gene} ne  $bedcoord{$chr}->{$myRegion2}{gene}) {
						$gene = $bedcoord{$chr}->{$myRegion1}{gene} .".". $bedcoord{$chr}->{$myRegion2}{gene}; 
					}
					$bedcoord{$chr}->{$regionName}{start}=&min($bedcoord{$chr}->{$myRegion1}{start},$bedcoord{$chr}->{$myRegion2}{start});
					$bedcoord{$chr}->{$regionName}{end}=&max($bedcoord{$chr}->{$myRegion1}{end},$bedcoord{$chr}->{$myRegion2}{end});				
					print "\tmerged into $regionName: $chr:$bedcoord{$chr}->{$regionName}{start}-$bedcoord{$chr}->{$regionName}{end}\n";
					delete($bedcoord{$chr}->{$myRegion1});
					delete($bedcoord{$chr}->{$myRegion2});
					$bedcoord{$chr}->{$regionName}{gene}=$gene;
					goto STARTCHR;

				}
			}	
	    	  }	

	 	}
	     }
	} elsif ($numberOfTargetedRegions<50000) { #merge amplicons overlapping for more than 75% ##### starting from v6.2. dont merge regions if there are more than 50000 of them!!!!
		for my $chr (keys (%bedcoord)) {
STARTCH2:
		  my @regions = sort keys (%{$bedcoord{$chr}});
		  for my $myRegion1 (@regions) {
			for my $myRegion2 (@regions) {
				if ($myRegion1 ne $myRegion2 && &overlap($bedcoord{$chr}->{$myRegion1}{start},$bedcoord{$chr}->{$myRegion1}{end},$bedcoord{$chr}->{$myRegion2}{start}, $bedcoord{$chr}->{$myRegion2}{end}) >= 0.75* &min(($bedcoord{$chr}->{$myRegion1}{end}-$bedcoord{$chr}->{$myRegion1}{start}),($bedcoord{$chr}->{$myRegion2}{end}-$bedcoord{$chr}->{$myRegion2}{start}))) {		
					print "These two regions overlap by more than 75%: $chr $bedcoord{$chr}->{$myRegion1}{start} - $bedcoord{$chr}->{$myRegion1}{end} and $bedcoord{$chr}->{$myRegion2}{start} - $bedcoord{$chr}->{$myRegion2}{end}\n"; 
					my $regionName=$myRegion1.".".$myRegion2;
					my $gene = $bedcoord{$chr}->{$myRegion1}{gene};
					if ($bedcoord{$chr}->{$myRegion1}{gene} ne  $bedcoord{$chr}->{$myRegion2}{gene}) {
						$gene = $bedcoord{$chr}->{$myRegion1}{gene} .".". $bedcoord{$chr}->{$myRegion2}{gene}; 
					}
					$bedcoord{$chr}->{$regionName}{start}=&min($bedcoord{$chr}->{$myRegion1}{start},$bedcoord{$chr}->{$myRegion2}{start});
					$bedcoord{$chr}->{$regionName}{end}=&max($bedcoord{$chr}->{$myRegion1}{end},$bedcoord{$chr}->{$myRegion2}{end});				
					print "\tmerged into $regionName: $chr:$bedcoord{$chr}->{$regionName}{start}-$bedcoord{$chr}->{$regionName}{end}\n";
					delete($bedcoord{$chr}->{$myRegion1});
					delete($bedcoord{$chr}->{$myRegion2});
					$bedcoord{$chr}->{$regionName}{gene}=$gene;
					goto STARTCH2;

				}
			}	
	    	  }	

	 	}


	}
}

###################################################
#read .txt for Controls and get region boundaries :
###################################################
if( $command eq "getSampleStats"){
	open FILE, "< $control " || die "$control : $!\n";
	while(<FILE>){
    		chomp;
		if (/(chr\S+)\s(\d+)\s(\d+)\s([\w\-\.]+)\s([\w\-\.]+)/) {
			$bedcoord{$1}->{$4}{start}=$2;
			$bedcoord{$1}->{$4}{end}=$3;
			$bedcoord{$1}->{$4}{gene}=$5;
		} elsif (/^(\S+)\s(\d+)\s(\d+)\s([\w\-\.]+)\s([\w\-\.]+)/) {
			$bedcoord{$1}->{$4}{start}=$2;
			$bedcoord{$1}->{$4}{end}=$3;
			$bedcoord{$1}->{$4}{gene}=$5;
		}
	}
	close FILE;

}

###################################################

print "\n------------------------\n\n";
print "\n--Coordinates are read--\n\n";
print "\n------------------------\n\n";


my %sortedRegionsByStartPos;
my %bedcoordArray;
for my $chr (keys (%bedcoord)) {
	my @regions = sort {$bedcoord{$chr}->{$a}{start} <=> $bedcoord{$chr}->{$b}{start}} keys (%{$bedcoord{$chr}});
	$sortedRegionsByStartPos{$chr} = \@regions;
	my $i=0;
	for my $myRegion (@regions) {
		$bedcoordArray{$chr}->[$i]{start}=$bedcoord{$chr}->{$myRegion}{start};
		$bedcoordArray{$chr}->[$i]{end}=$bedcoord{$chr}->{$myRegion}{end};
		$bedcoordArray{$chr}->[$i]{gene}=$bedcoord{$chr}->{$myRegion}{gene};
		$bedcoordArray{$chr}->[$i]{ampliName}=$myRegion;

		$i+=1;
	}
}

#################################
#read .bam for control sample(s):
#################################

if ($control=~/bam$/ && $command eq "getControlStats") { #start reading control .BAM

my @BAMControls = split(/\,/,$control);
my $Ncontrol=scalar(@BAMControls);

print "\tDetected $Ncontrol control sample(s)\n";
my %readNumbers;
my @controlFilenames;

for my $i (0..($Ncontrol-1)) {

 my $filename = $BAMControls[$i];
 open(FILE, "samtools view $filename|") or die "$0: can't open ".$filename.":$!\n";
 print "\treading $filename\n";
 $filename=~s/.bam//;
 $filename=~s/.*[\/]//;
 print "\tsample name: $filename\n";
 $controlFilenames[$i]=$filename;
 my $readCount = 0;$readNumbers{$i} = 0;
 my $currentRegion = "";
 my $currentRegionInd = 0;

 while(<FILE>){
    my @t=split;    
    next if ($t[0]=~/^@/);   
    next if scalar(@t)<5; 

    my ($chr,$spos,$readlen,$CIGAR)=($t[2],$t[3],length($t[9]),$t[5]);
    next if ($chr eq "*");
    $chr="chr".$chr unless ($chr =~/chr/);
    #my $mpos = $spos + floor($readlen/2);   
    my $epos =$spos+$readlen;
    $epos -=$1 if ($CIGAR=~/(\d+)S/);$readlen=$epos-$spos;
    $readCount+=1;
    print "read $readCount reads\n" if ($readCount%100000==0);    

    if ($base eq "Exon") {
	    if (exists($bedcoord{$chr}->{$currentRegion}) && $epos>=$bedcoord{$chr}->{$currentRegion}{start} && $spos<=$bedcoord{$chr}->{$currentRegion}{end}) {		
				$bedcoord{$chr}->{$currentRegion}{controlCount}{$i}+=1;   #print "Bingo\n";
				if ($bedcoord{$chr}->{$currentRegion}{start}>$spos)	{	
					#if ($bedcoord{$chr}->{$currentRegion}{start}-$spos>50) {
					#	print $bedcoord{$chr}->{$currentRegion}{start}-$spos,"\n";
					#	print "$currentRegion\n";print "$_\n";
					#}


					$bedcoord{$chr}->{$currentRegion}{start} = $spos; 
				}
				if ($bedcoord{$chr}->{$currentRegion}{end}<$epos) {		
					#if ($epos-$bedcoord{$chr}->{$currentRegion}{end}>50) {
					#	print $epos-$bedcoord{$chr}->{$currentRegion}{end},"\n";
					#	print "$currentRegion\n";print "$_\n";
					#}

					$bedcoord{$chr}->{$currentRegion}{end}=$epos;
				}
				$readNumbers{$i}+=1;
				next;
			

	    } else {
		my @regions = keys (%{$bedcoord{$chr}});
	    	if (scalar(@regions )<1) {next;}
	    	for my $myRegion (@regions) {
			if ($epos>=$bedcoord{$chr}->{$myRegion}{start} && $spos<=$bedcoord{$chr}->{$myRegion}{end}) {		
				$bedcoord{$chr}->{$myRegion}{controlCount}{$i}+=1;    #print "found in ", scalar(@regions)," regions\n";
				if ($bedcoord{$chr}->{$myRegion}{start}>$spos)	{	
					$bedcoord{$chr}->{$myRegion}{start} = $spos;
				}
				if ($bedcoord{$chr}->{$myRegion}{end}<$epos) {		
					$bedcoord{$chr}->{$myRegion}{end}=$epos;
				}
				$currentRegion = $myRegion;
				$readNumbers{$i}+=1;
				next;
			}	
	    	}	
	    }
     } else { #amplicon-based:
	 #need to find a region that the read overlaps the best
		my $maxOverlap=0; #print "Missed: $chr $spos",$currentRegion,"\n";
		my $bestIndex = -1;

		next unless (exists($sortedRegionsByStartPos{$chr}));

	 	for my $indexR (&max($currentRegionInd-5,0)..&min($currentRegionInd+5,scalar(@{$sortedRegionsByStartPos{$chr}})-1)) {
			my $o = overlap($spos,$epos,$bedcoordArray{$chr}->[$indexR]{start},$bedcoordArray{$chr}->[$indexR]{end});
			if ($o>$maxOverlap) {
				$maxOverlap = $o;
				$bestIndex=$indexR;								
			} elsif ($o==$maxOverlap && $o>0) {
				my $sPrev=abs($spos-$bedcoordArray{$chr}->[$bestIndex]{start});
				my $ePrev=abs($epos-$bedcoordArray{$chr}->[$bestIndex]{end});
				my $sCur=abs($spos-$bedcoordArray{$chr}->[$indexR]{start});
				my $eCur=abs($epos-$bedcoordArray{$chr}->[$indexR]{end});
				if (&min($sPrev,$ePrev)>&min($sCur,$eCur)) {
					$bestIndex=$indexR;
				}
			}
		}
		if ($maxOverlap==0) {
			#check all options;
			my @indexRtoCheck = &max(0,$currentRegionInd-20)..(&min($currentRegionInd+20,scalar(@{$sortedRegionsByStartPos{$chr}})-1)); 
			push(@indexRtoCheck,0..&min(5,scalar(@{$sortedRegionsByStartPos{$chr}})-1));

			for my $indexR (@indexRtoCheck) {
				next if ($epos<$bedcoordArray{$chr}->[$indexR]{start}); #since reads are ordered
				my $o = overlap($spos,$epos,$bedcoordArray{$chr}->[$indexR]{start},$bedcoordArray{$chr}->[$indexR]{end});
				if ($o>$maxOverlap) {
					$maxOverlap = $o;
					$bestIndex=$indexR;
					#print "found\n";				
				} elsif ($o==$maxOverlap && $o>0) {
					my $sPrev=abs($spos-$bedcoordArray{$chr}->[$bestIndex]{start});
					my $ePrev=abs($epos-$bedcoordArray{$chr}->[$bestIndex]{end});
					my $sCur=abs($spos-$bedcoordArray{$chr}->[$indexR]{start});
					my $eCur=abs($epos-$bedcoordArray{$chr}->[$indexR]{end});
					if (&min($sPrev,$ePrev)>&min($sCur,$eCur)) {
						$bestIndex=$indexR;
					}
				}
			}
		}

		next if ($maxOverlap==0); 
		my $bestRegion=$bedcoordArray{$chr}->[$bestIndex]{ampliName};
		$bedcoord{$chr}->{$bestRegion}{controlCount}{$i}+=1;    #print "found in ", scalar(@regions)," regions\n";				
		$currentRegionInd = $bestIndex;
		$readNumbers{$i}+=1;
		#print "$maxOverlap $readlen $t[0] $spos\n";
		next;
	

     }  
 }
 close FILE;

}

#check region boundaries and merge overlapping ones:
if ($base eq "Exon") {
for my $chr (keys (%bedcoord)) {
STARTCH3:
	my @regions = sort keys (%{$bedcoord{$chr}});
	for my $myRegion1 (@regions) {
		for my $myRegion2 (@regions) {
			if ($myRegion1 ne $myRegion2 && $bedcoord{$chr}->{$myRegion2}{end}>=$bedcoord{$chr}->{$myRegion1}{start} && $bedcoord{$chr}->{$myRegion2}{start}<=$bedcoord{$chr}->{$myRegion1}{end}) 				{		
				print "These two regions overlap: $chr $bedcoord{$chr}->{$myRegion1}{start} - $bedcoord{$chr}->{$myRegion1}{end} and $bedcoord{$chr}->{$myRegion2}{start} - $bedcoord{$chr}->{$myRegion2}{end}\n"; 
				my $regionName=$myRegion1.".".$myRegion2;

				my $gene = $bedcoord{$chr}->{$myRegion1}{gene};
				if ($bedcoord{$chr}->{$myRegion1}{gene} ne  $bedcoord{$chr}->{$myRegion2}{gene}) {
					if ($bedcoord{$chr}->{$myRegion1}{gene}=~m/$bedcoord{$chr}->{$myRegion2}{gene}/) {
						$gene = $bedcoord{$chr}->{$myRegion1}{gene};
					}elsif ($bedcoord{$chr}->{$myRegion2}{gene}=~m/$bedcoord{$chr}->{$myRegion1}{gene}/) {
						$gene = $bedcoord{$chr}->{$myRegion2}{gene};
					} else {
						$gene = $bedcoord{$chr}->{$myRegion1}{gene} .".". $bedcoord{$chr}->{$myRegion2}{gene};
					}
				}
				$bedcoord{$chr}->{$regionName}{gene}=$gene;
				$bedcoord{$chr}->{$regionName}{start}=&min($bedcoord{$chr}->{$myRegion1}{start},$bedcoord{$chr}->{$myRegion2}{start});
				$bedcoord{$chr}->{$regionName}{end}=&max($bedcoord{$chr}->{$myRegion1}{end},$bedcoord{$chr}->{$myRegion2}{end});				
				print "\tmerged into $regionName: $chr:$bedcoord{$chr}->{$regionName}{start}-$bedcoord{$chr}->{$regionName}{end}\n";
				for my $i (0..($Ncontrol-1)) {
					unless(exists($bedcoord{$chr}->{$myRegion1}{controlCount}{$i})) {
						$bedcoord{$chr}->{$myRegion1}{controlCount}{$i}=0;
					}
					unless(exists($bedcoord{$chr}->{$myRegion2}{controlCount}{$i})) {
						$bedcoord{$chr}->{$myRegion2}{controlCount}{$i}=0;
					}
					$bedcoord{$chr}->{$regionName}{controlCount}{$i}=$bedcoord{$chr}->{$myRegion1}{controlCount}{$i}+$bedcoord{$chr}->{$myRegion2}{controlCount}{$i};
				}
				delete($bedcoord{$chr}->{$myRegion1});
				delete($bedcoord{$chr}->{$myRegion2});
				goto STARTCH3;

			}
		}	
    	}	

}
}
#getTotal target length:
my $targetLength= 0;
for my $chr (keys (%bedcoord)) {
	my @regions = sort keys (%{$bedcoord{$chr}});
	for my $myRegion (@regions) {
		$targetLength+=$bedcoord{$chr}->{$myRegion}{end}-$bedcoord{$chr}->{$myRegion}{start}+1;
	}
}
print "\tTotal target length: $targetLength\n";

#print out summary for the control files:

open(FILE, ">$output_fname") or die "$0: can't open ".$output_fname.":$!\n";
print FILE "chr\tstart\tend\tID\tgene\tmean";

print "processed $Ncontrol controls, @controlFilenames\n";

for my $i (0..($Ncontrol-1)) {print FILE "\t$controlFilenames[$i]";}

print FILE "\n";	

for my $chr (sort bychr keys (%bedcoord)) {
	my @regions = @{$sortedRegionsByStartPos{$chr}};
	for my $myRegion (@regions) {
		my @values;
		if ($bedcoord{$chr}->{$myRegion}{end} <= $bedcoord{$chr}->{$myRegion}{start}) {
			print ERR "Warning: Incorrect coordinates of a region: $chr:",$bedcoord{$chr}->{$myRegion}{start},"-",$bedcoord{$chr}->{$myRegion}{end},"\n";
			print ERR "Warning: This region is going to be excluded\n";
			next;
		}
		for my $i (0..($Ncontrol-1)) {
			unless(exists($bedcoord{$chr}->{$myRegion}{controlCount}{$i})) {
				$bedcoord{$chr}->{$myRegion}{controlCount}{$i}=0;
			}			
			push(@values,$bedcoord{$chr}->{$myRegion}{controlCount}{$i}/( $readNumbers{$i}/$targetLength*($bedcoord{$chr}->{$myRegion}{end}-$bedcoord{$chr}->{$myRegion}{start}) ));		
		}

		my $mean = &mean(@values);
		print FILE "$chr\t$bedcoord{$chr}->{$myRegion}{start}\t$bedcoord{$chr}->{$myRegion}{end}\t$myRegion\t$bedcoord{$chr}->{$myRegion}{gene}\t$mean\t";
		print FILE join ("\t",@values);
		print FILE "\n";	
    	}	

}
close FILE;

} #end reading control .BAM

########################################################

#getTotal target length:
my $targetLength= 0;
for my $chr (keys (%bedcoord)) {
	my @regions = @{$sortedRegionsByStartPos{$chr}};
	for my $myRegion (@regions) {
		$targetLength+=$bedcoord{$chr}->{$myRegion}{end}-$bedcoord{$chr}->{$myRegion}{start}+1;
	}
}
print "\tTotal target length: $targetLength\n";

#################################
#read .bam for tumor sample(s):
#################################

my @tumorFilenames;

if ($sample=~/bam$/ && $command eq "getSampleStats") { #start reading sample .BAM

	my @BAMSamples = split(/\,/,$sample);
	my $Nsample=scalar(@BAMSamples);
	print "\tDetected $Nsample tumor sample(s)\n";
	my %readNumbers;

	for my $i (0..($Nsample-1)) {
		my $currentRegion = "";
		my $currentRegionInd = 0;
 		my $filename = $BAMSamples[$i];
 		open(FILE, "samtools view $filename|") or die "$0: can't open ".$filename.":$!\n";
 		print "\treading $filename\n";
		$filename=~s/.bam//;$tumorFilenames[$i]=$filename;
 		my $readCount = 0;$readNumbers{$i} = 0;
 		while(<FILE>){
    			my @t=split;    
    			next if ($t[0]=~/^@/);   
    			next if scalar(@t)<5; 
    			my ($chr,$spos,$readlen,$CIGAR)=($t[2],$t[3],length($t[9]),$t[5]);
    			next if ($chr eq "*");
    			$chr="chr".$chr unless ($chr =~/chr/);
			$readlen -=$1 if ($CIGAR=~/(\d+)S/);
    			my $epos =$spos+$readlen;
    			#my $mpos = $spos + floor($readlen/2);   
			$readCount+=1;
    			print "read $readCount reads\n" if ($readCount%1000000==0);

			if ($base eq "Exon") {
				if (exists($bedcoord{$chr}->{$currentRegion}) && $epos>=$bedcoord{$chr}->{$currentRegion}{start} && $spos<=$bedcoord{$chr}->{$currentRegion}{end}) {		
					$bedcoord{$chr}->{$currentRegion}{sampleCount}{$i}+=1;  
					$readNumbers{$i}+=1;
					next;
				} else {
					my @regions = keys (%{$bedcoord{$chr}});
	    				if (scalar(@regions )<1) {next;}
					my $found = 0;							
	    				for my $myRegion (@regions) {
						if ($epos>=$bedcoord{$chr}->{$myRegion}{start} && $spos<=$bedcoord{$chr}->{$myRegion}{end}) {		
							$bedcoord{$chr}->{$myRegion}{sampleCount}{$i}+=1;    
							$currentRegion = $myRegion;
							$readNumbers{$i}+=1;
							next;
						}	
	    				}					
	    			}    

			} else {
				#need to find a region that the read overlaps the best
				my $maxOverlap=0; #print "Missed: $chr $spos",$currentRegion,"\n";
				my $bestIndex = -1;

				next unless (exists($sortedRegionsByStartPos{$chr}));

			 	for my $indexR (&max($currentRegionInd-5,0)..&min($currentRegionInd+5,scalar(@{$sortedRegionsByStartPos{$chr}})-1)) {
					my $o = overlap($spos,$epos,$bedcoordArray{$chr}->[$indexR]{start},$bedcoordArray{$chr}->[$indexR]{end});
					if ($o>$maxOverlap) {
						$maxOverlap = $o;
						$bestIndex=$indexR;								
					} elsif ($o==$maxOverlap && $o>0) {
						my $sPrev=abs($spos-$bedcoordArray{$chr}->[$bestIndex]{start});
						my $ePrev=abs($epos-$bedcoordArray{$chr}->[$bestIndex]{end});
						my $sCur=abs($spos-$bedcoordArray{$chr}->[$indexR]{start});
						my $eCur=abs($epos-$bedcoordArray{$chr}->[$indexR]{end});
						if (&min($sPrev,$ePrev)>&min($sCur,$eCur)) {
							$bestIndex=$indexR;
						}
					}
				}
				if ($maxOverlap==0) {
					#check all options;
					my @indexRtoCheck = &max(0,$currentRegionInd-20)..(&min($currentRegionInd+20,scalar(@{$sortedRegionsByStartPos{$chr}})-1)); 
					push(@indexRtoCheck,0..&min(5,scalar(@{$sortedRegionsByStartPos{$chr}})-1));

					for my $indexR (@indexRtoCheck) {
						next if ($epos<$bedcoordArray{$chr}->[$indexR]{start}); #since reads are ordered
						my $o = overlap($spos,$epos,$bedcoordArray{$chr}->[$indexR]{start},$bedcoordArray{$chr}->[$indexR]{end});
						if ($o>$maxOverlap) {
							$maxOverlap = $o;
							$bestIndex=$indexR;
							#print "found\n";				
						} elsif ($o==$maxOverlap && $o>0) {
							my $sPrev=abs($spos-$bedcoordArray{$chr}->[$bestIndex]{start});
							my $ePrev=abs($epos-$bedcoordArray{$chr}->[$bestIndex]{end});
							my $sCur=abs($spos-$bedcoordArray{$chr}->[$indexR]{start});
							my $eCur=abs($epos-$bedcoordArray{$chr}->[$indexR]{end});
							if (&min($sPrev,$ePrev)>&min($sCur,$eCur)) {
								$bestIndex=$indexR;
							}
						}
					}
				}

				next if ($maxOverlap==0); 
				my $bestRegion=$bedcoordArray{$chr}->[$bestIndex]{ampliName};
				$bedcoord{$chr}->{$bestRegion}{sampleCount}{$i}+=1;    				
				$currentRegionInd = $bestIndex;
				$readNumbers{$i}+=1;
				#print "$maxOverlap $readlen $t[0] $spos\n";
				next;
			}
		}
 		close FILE;

	}
	#print out summary for the sample files:

	open(FILE, ">$output_fname") or die "$0: can't open ".$output_fname.":$!\n";
	print FILE "chr\tstart\tend\tID";
	for my $i (0..($Nsample-1)) {print FILE "\t$tumorFilenames[$i]";}
	print FILE "\n";	

	for my $chr (sort bychr keys (%bedcoord)) {
		next unless(exists($sortedRegionsByStartPos{$chr}));
		my @regions = @{$sortedRegionsByStartPos{$chr}};
		for my $myRegion (@regions) {
			my @values;
			for my $i (0..($Nsample-1)) {
				unless(exists($bedcoord{$chr}->{$myRegion}{sampleCount}{$i})) {
					$bedcoord{$chr}->{$myRegion}{sampleCount}{$i}=0;
				}
				push(@values,$bedcoord{$chr}->{$myRegion}{sampleCount}{$i}/( $readNumbers{$i}/$targetLength*($bedcoord{$chr}->{$myRegion}{end}-$bedcoord{$chr}->{$myRegion}{start}) ));			
			}
			print FILE "$chr\t$bedcoord{$chr}->{$myRegion}{start}\t$bedcoord{$chr}->{$myRegion}{end}\t$myRegion\t";
			print FILE join ("\t",@values);
			print FILE "\n";	
    		}	

	}
	close FILE;

} #end reading Sample BAMs
########################################################
sub overlap {
	my($a,$b,$x,$y) = @_;
	return 0 if ($b<$x || $a>$y);
	return min($b,$y) - max($a,$x);
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub mean {
    my $result;
    foreach (@_) { $result += $_ }
    return $result / @_;
}


sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub round {
    my($number) = shift;
    return int($number + .5);
}

sub bychr {
 my $chr1 = $a;
 my $chr2 = $b;
 $chr1=~s/chr//;
 $chr2=~s/chr//;
 $chr1=~s/X/23/;
 $chr1=~s/Y/24/;
 $chr1=~s/M/25/;
 $chr2=~s/X/23/;
 $chr2=~s/Y/24/;
 $chr2=~s/M/25/;
 $chr1 <=>$chr2 ;
}

sub max ($$) { $_[$_[0] < $_[1]] }

sub min ($$) { $_[$_[0] > $_[1]] }


