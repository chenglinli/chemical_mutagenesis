Authors: Cheng-Lin Frank Li (chenglil@bcm.edu, chenglinfli@gmail.com)  
Date: October 20, 2015

## Introduction
These scripts written in R are genetic variant annotation tools for _Dictyostelium_. It is designed for filtering and annotating a SNP calling output (.vcf file) by the GATK UnifiedGenotyper (version 3.4-46). The README file contains guidelines for 1) how to prepare the Variant Call Format (VCF) file and 2) examples of applying the scripts to VCF files. Be aware of that the tools are designed for the VCF file format version 4.1 and the format of the VCF file are subjected to changes. If you encounter any error, feel free to send me an email.

## Part I: Prepare the Variant Call Format (VCF) file
This section will walk you through the steps and tools to generate a VCF file (.vcf) that is compatible with these scripts. If you create a VCF file by using different tools, you may encounter unexpected errors. For more information about the VCF file, please read the [document](http://gatkforums.broadinstitute.org/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it). 

**1. Align the sequence reads**  
First, you will align the whole-genome sequence data (.bam files) to the reference genome (mask the duplication on chromosome 2 from 3016083 to 3768654) by the BWA (0.7.5a). The chromosome 2 duplication region will be addressed later. If you use Bowtie to align your sequence reads, you may have to run `java -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fa -I input.bam -o output.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60` to covert the quality score before variant calling.

**2. Mark duplicated reads**  
You will use the `MarkDuplicates` in the picard-tools (version 1.119) to mark duplicated reads. These duplicated reads will later be excluded from the variant calling.

**3. Add read group**
The GATK requires to have read group information in BAM files. So you have to run the `AddOrReplaceReadGroups` in the picard-tools (version 1.119) to add read groups. Here is an example of the read groups.

RGID | RGLB | RGPL | RGPU | RGSM
--- | --- | --- | --- | --- 
strain01 | aggless | illumina | 2x100bp | AX4
strain02 | aggless | illumina | 2x100bp | mutant01
strain03 | aggless | illumina | 2x100bp | mutant02

As you might have noticed, the RGID and RGSM names must be unique for different samples. The strain01 is the parental strain, AX4. The strain02 and strain03 are chemically mutagenized strains that carry unknown mutations. The RGSM will be the column names of individual strains in the VCF file. Please name the RGSM of your mutant strains with a common prefix (e.g. “mutant” in this example). 

**4. Variant calling**  
Lastly, running the UnifiedGenotyper on multiple samples of the same experiment to generate a single VCF file (**Note1**). You will generate two VCF files for **SNVs** and **indels**, respectively, across the whole genome excluding the chromosome 2 duplication. Here are examples of the commands for running GATK in MAC OS X (v10.10.5). If you use other system or cloud server, make sure you have the same values for the two parameters (-stand_call_conf and -ploidy).  
**Note1**: Do not run separate variant calling on individual samples and combine the individual vcf files later.

**The VCF file contain the single nucleotide variants (SNVs)**  
Calling SNVs for a haploid genome. You will later use the `snv.R` to filter and annotate the VCF file.
```
java -jar gatk_directory/GenomeAnalysisTK.jar -T UnifiedGenotyper -R reference_fasta_file -I input1.bam [-I input2.bam ...] -o snv.vcf -stand_call_conf 30 -ploidy 1
```

**The VCF file contain the indels**  
Calling indels for a haploid genome. You will later use the `indel.R` to filter and annotate the VCF file.
```
java -jar gatk_directory/GenomeAnalysisTK.jar -T UnifiedGenotyper -glm INDEL -R reference_fasta_file -I input1.bam [-I input2.bam ...] -o indel.vcf -stand_call_conf 30 -ploidy 1
```

## Part II: Examples of applying the scripts to VCF files.

**1. The reference files**  
There are six reference files providing essential information for annotating variants. The reference files have to be placed in a given folder.

* gene_09-26-2015.txt  
This file provides the gene model. It looks like this:
```
chromosome start  end strand        ddb_g         genename
chr1 1890 3287      + DDB_G0267178 DDB_G0267178_RTE
```

* dicty_primary_cds.txt  
This file provides the dictybase ID and the coding sequence of each gene. It looks like this:
```
DDB0216524|DDB_G0267364
"a" "t" "g" "g" "t" "c" "g" "t" "a" "c" "c" "a" "g" "a" "c" "t" "g" "g" "t" "t" "a" "c" "a" "a" "a" "c" "t" "g" "a" "t" "g"
 ....
```

* cds_09-26-2015.gff  
This file is a General Feature Format (GFF) file of coding DNA sequence (CDS). It looks like this:
```
seqname            source feature start  end score strand frame         attribute chromosome
DDB0232428 Sequencing Center     CDS  1890 3287     .      +     . Parent=DDB0216437       chr1
```

* splice_acceptor.gff  
This file is a General Feature Format (GFF) file of splicing acceptor sites of the introns. It looks like this:
```
seqname	source	feature	start	end	score	strand	frame	attribute	chromosome	n
DDB0232428	Sequencing Center	CDS	11262	11263	.	+	.	DDB0216442	chr1	3
```

* splice_donor.gff  
This file is a General Feature Format (GFF) file of splicing donor sites of the introns. It looks like this:
```
seqname	source	feature	start	end	score	strand	frame	attribute	chromosome	n
DDB0232428	Sequencing Center	CDS	11200	11201	.	+	.	DDB0216442	chr1	3
```

* dd_codon_bias.txt  
This file provides the frequency of codon usage.
```
Codon	AA	Number	freq	freq_abs
GCA	Ala	109859	1.67	1.668609801
```

**2. Applying the scripts to your VCF files**  
Here is an exmaple of applying the `snv.R` to filter and annotate the `example_snv.vcf`. The VCF file contains 39 SNVs as a result of whole-genome SNVs calling. First, let's take a look of the `example_snv.vcf`. If you are not familiar with the format, please read the [document](http://gatkforums.broadinstitute.org/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it)

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	AX4	mutant01	mutant02	mutant03	mutant04	mutant05	mutant06
chr6	773798	.	A	G	5181.7	.	AC=1...	GT:AD:DP:GQ:PL	0:68,0:68:99:0,2891	...
```

The parental strain is the "AX4". Mutant strains is named "mutant01" to "mutant06". To filter and annotate the SNVs. You will run the following commands in R.
```
source("/path_to_the_folder/example_snv.vcf")
snv(input_file = "/path_to_the_folder/example_snv.vcf", ref_file = "/path_to_the_folder_of_reference_files/", parental_strain = "par", mutant_strain = "mutant")
```

The script creates a new folder `example_snv.vcf_5`. The suffix `5` indicates that the SNVs must be covered by >= 5 reads. You could change the threshold by adding the argument `read_depth = 5` in the function `snv()`. The script generates 8 files.

**File name** | **Description** 
--- | --- 
summary.txt | Summarize the input parameters and results.
variant_filtered.txt | A data frame of SNVs that passed the filters.
variant_mult_alt.txt | The SNVs that contain more than two alternative variants. These SNVs are likely false positive but saved in the file just in case.
gene_list_all.txt | Genes that are mutated.
gene_list_top.txt | Genes that are mutated at least twice.
mutations_by_chr.txt |  List mutations in individual chromosomes.
mutations_by_strain.txt |  List mutations in individual strains.
strain_by_gene.txt | List mutants that carrying mutations in individual genes.

* summary.txt 
The file should look like this:
```
input_file = ...
parental_strain = AX4
mutant_strain = mutant
...
88: total number of variant sites before filtering
87: remaining variants after removing sites where the parental strain has no read (filter 1).
86: remaining variants after filter 2 (genotype filter). Allowed genotypes : ., 0, 1.
66: remaining variants after filter 3.1 (frequency filter)
42: remaining variants after filter 3.2 (freq_diff & read_depth filters).
1: variants with multiple alternative alleles that pass filter 3.2 (freq_diff & read_depth filters).
...
```
There are 88 SNVs before filtering. One SNVs did not pass the filters because it contains multiple alternative alleles. 42 SNVs pass all filters (36 SNVs are found in one mutant and 3 SNVs are found in two mutants, 36 + 3*2 = 42). The 42 SNVs are listed in the `variant_filtered.txt`. The results can be found in the folder `example_snv.vcf_5`.

Similarly, you can apply the `indel.R` to filter and annotate the `example_indel.vcf`. The `example_indel.vcf` contains 90 indels but only 2 of them pass the filters. No variant with multiple alternative alleles are found. The results can be found in the folder `example_indel.vcf_5`.

## Part III: Chromosome 2 duplication
In the Part I & II, we treat the Dictyostelium genome exluding the chromosome 2 duplication as a haploidy. In Part III, we will do variant calling, filtering and annotation on chromosome 2 duplication region by considering it as a diploidy. All the steps are basically the same as the Part I & II. The differences are adding the arguments `-ploidy 2 -L chr2:2263132-3015703` during vairant calling and use the R scripts (`snv_chr2.R` and `indel_chr2.R`) for filtering and annotation.

**1. Variant calling** 
**The VCF file contains single nucleotide variants (SNVs) on the chromosome 2 duplication**  
Calling SNVs for a diploid genome (There are two copies of each gene in the chromosome 2 duplication). You will later use the  `snv_chr2.R` to filter and annotate the VCF file.
```
java -jar gatk_directory/GenomeAnalysisTK.jar -T UnifiedGenotyper -R reference_fasta_file -I input1.bam [-I input2.bam ...] -o snv_chr2.vcf -stand_call_conf 30 -ploidy 2 -L chr2:2263132-3015703
```

**The VCF file contains indels on the chromosome 2 duplication**  
Calling indels for a diploid genome (There are two copies of each gene in the chromosome 2 duplication). You will later use the  `indel_chr2.R` to filter and annotate the VCF file.
```
java -jar gatk_directory/GenomeAnalysisTK.jar -T UnifiedGenotyper  -glm INDEL -R reference_fasta_file -I input1.bam [-I input2.bam ...] -o indel_chr2.vcf -stand_call_conf 30 -ploidy 2 -L chr2:2263132-3015703
```

No example is provided for this part.
