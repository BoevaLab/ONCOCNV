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

VERSION=6.9

###############################################################################################################

# intput: .bam files 

#---------------------------------------------------------------------------------------------------------------
#  you may want to create .bai indexes to visualize .bam files in the IGV:
#     files="Control.1.bam,Control.2.bam,Control.3.bam,Sample.1.bam,Sample.2.bam,Sample.3.bam,Sample.4.bam"
#     for file in $files
#     do
#       samtools index $file
#     done
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#                              Set paths and arguments
#---------------------------------------------------------------------------------------------------------------

#alias bedtools=YOURPATH/BEDTools/bin/bedtools          #add alias to bedtools if needed
#alias samtools=YOURPATH/samtools/bin/samtools          #add alias to bedtools if needed


TOOLDIR=../scripts/               #set a path to the directory with CNV detection scripts. In your case, it can be, for example, ONCOCNV.v5.5/
DATADIR=../BAMfiles/              #set a path to the directory with .BAM files
OUTDIR=../outputDEEPCNA/          #set a path to the *existing* directory to write the output
GENOME=../hg19/hg19.fa            #set a path to the .fasta file with the reference genome

targetBed="4477685_CCP_designed.bed"      #.bed file with start and end position of amplicons (amplicons overlaping for more than 75% will be automatically merged)


#set control and test samples (located in $DATADIR):

controls="Control.1.bam,Control.2.bam,Control.3.bam"
tests="Sample.1.bam,Sample.2.bam,Sample.3.bam,Sample.4.bam"

#---------------------------------------------------------------------------------------------------------------
#                                   Run commands
#---------------------------------------------------------------------------------------------------------------

echo "ONCOCNV $VERSION - a package to detect copy number changes in Deep Sequencing data"

cd $DATADIR

echo "$DATADIR contains the following files:"
ls -l $DATADIR 

#get normalized read counts for the control samples
echo "running ONCOCNV_getCounts.pl on control samples"
perl $TOOLDIR/ONCOCNV_getCounts.pl getControlStats -m Ampli -b $targetBed -c $controls -o $OUTDIR/Control.stats.txt 
echo "$OUTDIR/Control.stats.txt was created"
ls -l $OUTDIR/Control.stats.txt

#get normalized read counts for the tumor samples
echo "running ONCOCNV_getCounts.pl on tumor samples"
perl $TOOLDIR/ONCOCNV_getCounts.pl getSampleStats -m Ampli -c $OUTDIR/Control.stats.txt -s $tests -o $OUTDIR/Test.stats.txt 
echo "$OUTDIR/Test.stats.txt was created"
ls -l $OUTDIR/Test.stats.txt

#create .bed file with targeted regions
echo "creating target.bed"
cat $OUTDIR/Control.stats.txt | grep -v start | awk '{print $1,$2,$3}' | sed "s/ /\t/g" >$OUTDIR/target.bed
ls -l $OUTDIR/target.bed
 
#get GC-content per targeted region 
echo "creating target.GC.txt"
perl $TOOLDIR/createTargetGC.pl -bed $OUTDIR/target.bed -fi $GENOME -od $OUTDIR -of $OUTDIR/target.GC.txt
ls -l $OUTDIR/target.GC.txt

#process control samples
echo "running processControl.R"
cat $TOOLDIR/processControl.R | R --slave --args $OUTDIR/Control.stats.txt $OUTDIR/Control.stats.Processed.txt $OUTDIR/target.GC.txt

#uncomment in case you want to manually limit the number of principal components
#PCtoKeep=1
#cat $TOOLDIR/processControl.v$VERSION.R | R --slave --args $OUTDIR/Control.stats.txt $OUTDIR/Control.stats.Processed.txt $OUTDIR/target.GC.txt $PCtoKeep
ls -l $OUTDIR/Control.stats.Processed.txt

#process test samples and predict CNA and CNVs:
echo "running processSamples.R"
#to use cghseg segmentation instead of circular binary segmentation (defaut) please add "chgseg" at the end of the command line:
cat $TOOLDIR/processSamples.R | R --slave --args $OUTDIR/Test.stats.txt $OUTDIR/Control.stats.Processed.txt $OUTDIR/Test.output.txt cghseg

######### segment data with binary segmentation (it is a default option but often cghseg performs better than the default circular binary segmentation): ########################################
#using DNAcopy circular binary segmentation (defaut) :
#cat $TOOLDIR/processSamples.R | R --slave --args $OUTDIR/Test.stats.txt $OUTDIR/Control.stats.Processed.txt $OUTDIR/Test.output.txt
#################################################################################################################################################################################################

