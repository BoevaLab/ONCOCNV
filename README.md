# ONCOCNV


**ONCOCNV - a package to detect copy number changes in Deep Sequencing data**


REQUIREMENTS

0. Perl and R installed and added to the PATH
	E.g., export PATH=$PATH:YOURPATH/R/bin	
1. SAMtools (http://samtools.sourceforge.net/) installed and added to the PATH
	To add to PATH, type in the command line or add to "ONCOCNV.sh": 
	export PATH=$PATH:YOURPATH/samtools/bin
		or
	alias samtools=YOURPATH/samtools/bin/samtools	
2. BEDTools (http://bedtools.readthedocs.org/en/latest/) installed and added to the PATH
	To add to PATH, type in the command line or add to "ONCOCNV.sh": 
		export PATH=$PATH:YOURPATH/BEDTools/bin/
			or
	alias bedtools=YOURPATH/BEDTools/bin/bedtools
3. The following R libraries should be installed: MASS, mclust, PSCBS, DNAcopy, R.cache, scales, cwhmisc, fastICA, cghseg, digest
4. The fasta sequence (one file, unzipped; e.g. "hg19.fa") of the targeted genome should be downloaded from http://hgdownload.soe.ucsc.edu/downloads.html
5. You need to have your data aligned (.bam files)
6. You need to have at least **three** control files to construct a **reliable** baseline. However, ONCOCNV will run with only 2 controls starting from version 5.4 and with JUST **one** control starting from version 5.7. Yet, we recommend to have at least 3 control for good performance of the algorithm.

INSTALLATION

0. Download ONCOCNV.zip (or ONCOCNV.vX.X.zip)
1. Unzip files into detectory "scripts"
2. Check requirements (R + the necessary R packages must be installed)
   To install the necessary P packages (when R is installed), type in the command line:
   
	   R
	   install.packages("MASS")
	   install.packages("mclust")
	   install.packages("R.cache")
	   install.packages("scales")
	   install.packages("cwhmisc")
	   install.packages("fastICA")
	   install.packages("cghseg")
	   install.packages("digest")
	   source("http://bioconductor.org/biocLite.R")
	   biocLite("DNAcopy")
	   install.packages("PSCBS")
	   quit()


RUN ONCOCNV

0. Open "ONCOCNV.sh" with a text editor (gedit, textpad, etc.)
1. Set correct paths and filenames in the top part of the "ONCOCNV.sh"
2. Check properties of "ONCOCNV.sh"
	chmod +rwx PathToONCOCNV/scripts/ONCOCNV.sh
3. Check formats: 
	o	reads should be given in .BAM format
	o	amplicon coordinates should be given in .bed format (with or without the headline) and have amplicon ID in column 4 and gene symbol in column 6, e.g.:
		chr1	2488068	2488201	AMPL223847	0	TNFRSF14
		It is mandatory to provide gene names in the 6th column.
		
---------------------------------------------------------------------------------------------------------------------------

****************VERY IMPORTANT****************

		Please make sure that:
	-	There is no duplicates in the coordinates
	-	Coordinates are sorted
	-	Gene names are gene names in the sense that corresponding amplicons fall in the same genomic locus and not on different chromosomes
	-	Gene names cannot be the same as amplicon names or IDs because ONCOCNV assumes to have several amplicons per gene

---------------------------------------------------------------------------------------------------------------------------

4. Run "ONCOCNV.sh" from the command line:
	cd PathToONCOCNV/scripts
	./ONCOCNV.sh	
		or		
	. PathToONCOCNV/scripts/ONCOCNV.sh
		
HOW TO READ OUTPUT FILES

There are three output files per sample:
1. *.profile.png
   - Visual representation of normalized and annotated copy number profile
	Each dot corresponds to an amplicon; the X-axis is not up to scale.
	Color code:
		o	GREEN				one-point-outlier
		o	DARK GREY SURROUNDINGS		frequent one-point-outlier
		o	BROWN				>1 level gain
		o	BROWN SURROUNDINGS		1-level gain
		o	BLUE				>1 level loss
		o	BLUE SURROUNDINGS		1-level loss
	
2. *.summary.txt
   - predictions per gene

gene		gene name
chr		chromosome name
start		first amplicon start
end		last amplicon start
copy.number	predicted copy number (no normal contamination nor subclones is taken into accout)
p.value		p-value for the copy number status of the genomic region encompassing the gene 
q.values	q-value ("fdr"-corrected p-value) for the copy number status of the genomic region encompassing the gene 
comments	p-value for the hypothesis that the copy number of the gene does not match the copy number of the encompassing segment 
		(in the case of a break within a gene - it is the p-value for the break) 

3. *.profile.txt
   - predictions per amplicon (detailed information)

chr			chromosome name
start			first amplicon start
end			last amplicon start
gene			gene name
ID			amplicon ID
ratio			logarithm of the normalized read count (zero values correspond to the neutral copy number)
predLargeSeg		copy number predicted by segmentation of normalized read counts
predLargeCorrected	*final prediction for the copy number*
pvalRatioCorrected	p-value of the t-test to test the difference between the normalized read counts and the value expected from the segmentation or from the gene-based copy number assessment
perGeneEvaluation	copy number predicted per gene (unaware of the segementation)
pvalRatioGene		gene-based p-value of the t-test for the difference of the mean of the normalized read counts from zero
predPoint		predicted one-point-outlier 
predPointSusp		predicted (frequent) one-point-outlier 
comments		additional information: 
				SegRatio	mean value of the logarithms of the normalized read counts per segment
				AbsMeanSigma	normalized difference of the mean value (~z-score/sqrt(#amplicons in the segment))
				pvalue		p-value for AbsMeanSigma
				pvalueTTest	p-value of the t-test (per segment)
