# assists in finding formatting errors in make file: 
# cat -e -t -v  Makefile
annotFile:=saccharomyces_cerevisiae.gff
# please place bam files for analysis in directory called "bam" or else rename the bam_folder variable to the file path to a folder with bam files for analysis
bam_folder:=bam
export bam_folder

#help:  @List available tasks on this project
help: 
	@echo ">--------------------------------------------------------------------------------------------------------------------------------------------------<"
	@echo "This makefile will calculate an Intron/Exon Junction Splice Ratio (IEJR) from bam files." \
	"\nThis script uses bedtools to count reads spanning intron/exon boundaries."
	@echo "\nThe rules:"
	@grep -E '[a-zA-Z\.\-]+:.*?@ .*$$' $(MAKEFILE_LIST)| sort | tr -d '#'  | awk 'BEGIN {FS = ":.*?@ "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
	@echo "\nYou can specify a bam folder for analysis by changing the bam_folder variable in the makefile"
	@echo "bam_folder is currently: "${bam_folder}
	@echo ">--------------------------------------------------------------------------------------------------------------------------------------------------<"

#annotationFiles: @ Uses gff inputfile of introns and makes a new GFFs of bases directly after the 5' splice site 
annotationFiles:
	# please specify the annotation gff using the annotFile variable in the makefile default is "saccharomyces_cerevisiae.gff"
	# pulls intron features out of the SGD gff file w/o mito reads
	grep "	intron" ${annotFile} | grep -v "mt	" > intron.gff
	# these commands split reference gff files into + and - strands
	awk '$$7 == "-" {print}' intron.gff > intron_minus.gff
	awk '$$7 == "+" {print}' intron.gff > intron_plus.gff
	# this defines a gff file with the first intronic basepair in the feature
	awk '{$$4=$$5}1' OFS='\t' intron_minus.gff > intronBorder_minus.gff
	awk '{$$5=$$4}1' OFS='\t' intron_plus.gff > intronBorder_plus.gff
	cat intronBorder_minus.gff intronBorder_plus.gff | bedtools sort > intronBorder.gff
	-mkdir intermediateAnnotFiles
	mv intron.gff intermediateAnnotFiles
	mv *minus* intermediateAnnotFiles
	mv *plus* intermediateAnnotFiles

#count: @ Counts the number of total and unspliced reads at the +1 bp of the 5' splice site 
count:
	bash src/countIntersect.sh
	
	
#IEJR: @ Uses R to calculate IEJR. Finish visualization in R script. 
IEJR:
	Rscript src/calculateIEJR.R
	@echo "Please see src/calculateIEJR.R to find examples of IEJR visualization"
	