#Here is a *small* test set to test the functionality of ONCOCNV:

######################################################
#  to run the whole pipeline refer to "ONCOCNV.sh"   #
######################################################

#change paths if needed!!!

DATADIR=.
OUTDIR=.
TOOLDIR=..

cd $DATADIR


cat $OUTDIR/Control.stats.txt | grep -v start | awk '{print $1,$2,$3}' | sed "s/ /\t/g" >$OUTDIR/target.bed
 
cat $TOOLDIR/processControl.R | R --slave --args $OUTDIR/Control.stats.txt $OUTDIR/Control.stats.Processed.txt $OUTDIR/target.GC.txt

cat $TOOLDIR/processSamples.R | R --slave --args $OUTDIR/Test.stats.txt $OUTDIR/Control.stats.Processed.txt $OUTDIR/Test.output.txt cghseg
