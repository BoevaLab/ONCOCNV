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
use strict;
use POSIX;

my $usage = qq{    

$0 getControlStats -f <fastaFile>
	
};

my $filename = "";
while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV; 
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }
    elsif ( $this_arg eq '-f') {$filename = shift @ARGV;} 
    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}
print STDERR $usage if ($filename eq "");

open FILE, "< $filename " || die "$filename : $!\n";

print "region\tGC\n";

my ($name,$GC);
while(<FILE>){
    chomp;
    if (/>(.+)/) {
	$name=$1;
    } else {
	my $len = length($_);
	$GC=(tr/GCgc//) / ($len-(tr/nN//));
	print "$name\t$GC\n";
    }
}
close FILE;


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub round {
  my $number = shift || 0;
  my $dec = 10 ** (shift || 0);
  return int( $dec * $number + .5 * ($number <=> 0)) / $dec;
}
	
