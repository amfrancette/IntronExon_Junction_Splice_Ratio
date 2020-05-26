# this script separately counts total and unspliced reads of a bam file over the 5' splice site of annotated introns 
# individual count files are made for each sample alongside summary files for the total and unsplcied
# a sample name file is generated to keep track of which column in the output files are which samples
# the output should be directly usable in the accompanying R script to generate Intron Exon Junction Enrichment Scores 
# visualization will be up to the user but some examples are provided. 

# pulls the names of the parent features holding the introns 
# I should fix this up to be more reliable somehow
awk -F '[=;]' '{print $2}' intronBorder.gff | awk -F '_' '{print $1}' > res/featurenames.txt
# grep -f res/featurenames.txt sacharomyces_cerevisiae.gff

# adds labels to the features associated with the introns 
sed 's/^Y/gene	Y/' res/featurenames.txt > res/temp.tab
sed 's/^t/tRNA	t/' res/temp.tab > res/temp2.tab
sed 's/^s/snoRNA	s/'  res/temp2.tab > res/featurenames.txt

# removes any old sample_names file from previous run and makes a new, empty file
rm res/sample_names.tab
touch res/sample_names.tab

# builds skeleton of a summary file with just the feature names to later be appended with total and unspliced count data 
cat res/featurenames.txt > res/totalSpliceJunctionCount.tab
cat res/featurenames.txt > res/unsplicedSpliceJunctionCount.tab

# loops through folder of bamfiles (the variable bam_folder is specified in the makefile: line 4) and counts reads that fall on the +1 bp of all 5' splice junctions
for bamfile in `ls $bam_folder | grep -v "bai" | awk -F '[./]' '{print $1}'`
do
	# creates a txt file keeping track of order in which files are fed through this loop
	echo $bamfile >> res/sample_names.txt
	# counts all reads (aligned sense with respect to the transcript [-s]) that touch the annotated splice site
	bedtools coverage -s -counts -a intronBorder.gff -b $bam_folder/"$bamfile"*.bam | awk '{print $NF}' > res/"$bamfile"_total.tab
	# counts only unspliced [--split] reads (aligned sense with respect to the transcript [-s]) that touch the annotated splice site
	bedtools coverage -s -counts -split -a intronBorder.gff -b $bam_folder/"$bamfile"*.bam | awk '{print $NF}' > res/"$bamfile"_unspliced.tab
	# appends the total counts to a summary matrix
	paste -d '\t' res/totalSpliceJunctionCount.tab res/"$bamfile"*total.tab > res/temp.tab
	mv res/temp.tab res/totalSpliceJunctionCount.tab
	# appends the unspliced counts to a summary matrix
	paste -d '\t' res/unsplicedSpliceJunctionCount.tab res/"$bamfile"*unspliced.tab > res/tempfile.txt
	mv res/tempfile.txt res/unsplicedSpliceJunctionCount.tab
done
exit

# cleans up some temporary files
rm temp*.tab