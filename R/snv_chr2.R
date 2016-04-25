#' Annotate SNVs
#'
#' @param input_file file path of the vcf file (e.g. input_file = "./test.vcf")
#' @param parental_strain the name of the parental strain. The name should not contain -. (e.g. parental_strain = "par")
#' @param mutant_strain please use a prefix to name mutant strains. Specify the common prefix here (e.g. mutant_strain = "mut" given mutant strains named as mut1, mut2, etc.)
#' @param read_depth at a given variant, the read depth (coverage) must be >= read_depth. Default value = 5.
#' @param var_thld if the reference allele was detected in <= var_thld (default value = 0.2) of the reads, it is
#'                  considered as having a variant allele.
#' @param ref_thld if the reference allele was detected in >= ref_thld (default value = 0.8) of the reads, it is 
#'                 considered as having a reference allele.
#' @param freq_diff at a given variant, the porportion of reads supports the reference allele in the parental
#'                     and the mutant strain should differ >= freq_diff (default value = 0.8).
#'
#' @return None
#'
#' @author Cheng-Lin Frank Li, \email{chenglinfli@gmail.com}
#'
#' @examples
#' indel(input_file = "/path_to_the_folder/snv.vcf", ref_file = "/path_to_the_folder_of_reference_files/",
#'      parental_strain = "AX4", mutant_strain = "mutant")
#'
#' @export
#' @import Rsamtools
#' @import IRanges
#' @import GenomicRanges
#' @importFrom seqinr read.fasta comp translate
#' @importFrom dplyr select mutate filter
#' @importFrom tidyr separate
#' @importFrom reshape2 colsplit
#' @importFrom stringr str_count
#'
snv_chr2 <- function(input_file = NULL, ref_file = system.file("extdata", package="chemut"), parental_strain = NULL, mutant_strain = NULL,
                     read_depth = 5, var_thld = 0.2, ref_thld = 0.8, freq_diff= 0.4){
    for(rd in read_depth){
        read_depth = rd
        
        StartTime = Sys.time()
        
        # Check input parameters
        if(is.null(input_file)) {print("Please specify the path of input_file")}
        if(is.null(parental_strain)) {"Please specify the name of the parental strain"}
        if(is.null(mutant_strain)) {"Please specify the name of the mutant strains"}
        
        # The first column in the vcf file is named "#CHROM". The hashtag(#) should be removed.
        ## the new vcf is saved as xxx_edited.tab
        input_file_edited = paste(input_file,"_edited.tab",sep="")
        if(file.exists(input_file_edited)) {
            vcf_ <- read.table(input_file_edited, header = T, sep = "\t", stringsAsFactors = F)
            print("Load edited vcf file")
        } else {
            vcf_ <- readLines(input_file) # Read the vcf file
            pos_ <- grep("#CHROM",vcf_) # Find the column contains "#CHROM"
            if (length(pos_) != 0) {
                vcf_[pos_] <- sub(pattern = "#CHROM", replacement = "CHROM", vcf_[pos_])
            } else {print("Header had been edited.")}
            
            vcf_edit <- paste(input_file,"_edited.tab",sep="") # Name of edited vcf file
            write.table(vcf_, vcf_edit, sep="\t", quote = F, row.names=F, col.names = F) 
            vcf_ <- read.table(vcf_edit, header = T, sep = "\t", stringsAsFactors = F)
        }
        
        # Create directory that contains the output files
        out_dir = paste(input_file,as.character(read_depth), sep="_")
        dir.create(out_dir)
        print(paste("Results will be written to :", out_dir, sep=""))
        
        # Names of strains
        parental_strain = as.character(parental_strain) # Parental strain
        mutant_strain = as.character(mutant_strain)  # Mutant strain
        # Vectors contain the name of the parental strain or all mutant strains
        names_mut <- grep(mutant_strain, colnames(vcf_), value = T) # Names of the mutant strains
        if(length(names_mut)==0){stop("check the name of the mutant strains")}
        names_par <- grep(parental_strain, colnames(vcf_), value = T) # Names of parental strain
        if(length(names_par)==0){stop("check the name of the parental strains")}
        names_all <- c(names_par,names_mut) # A vector contains names of all strains
        
        # Load the gene position data, which contains the details of each gene (chromosome, gene_start, gene_stop and strand).
        dd_gene_pos <- read.table(file.path(ref_file, "gene_09-26-2015.txt.gz"), sep="\t", header=T, stringsAsFactors = FALSE)
        dd_gene_pos <- dd_gene_pos[grep("DDB_G", dd_gene_pos$ddb_g),]
        dd_gr <- GenomicRanges::GRanges(dd_gene_pos$chromosome, IRanges(start = dd_gene_pos$start, end = dd_gene_pos$end))
        
        # write information to the summary.txt
        write(paste("input_file =", input_file),paste(out_dir,"/","summary.txt",sep=""), append=T)
        write(paste("parental_strain =",parental_strain),paste(out_dir,"/","summary.txt",sep=""), append=T)
        write(paste("mutant_strain =",mutant_strain),paste(out_dir,"/","summary.txt",sep=""), append=T)
        write(paste("freq_diff =",freq_diff),paste(out_dir,"/","summary.txt",sep=""), append=T)
        write(paste("read_depth =",read_depth),paste(out_dir,"/","summary.txt",sep=""), append=T)
        write(paste("var_thld = ",var_thld),paste(out_dir,"/","summary.txt",sep=""), append=T)
        write(paste("ref_thld = ",ref_thld),paste(out_dir,"/","summary.txt",sep=""), append=T)
        write("", paste(out_dir, "/", "summary.txt",sep=""), append=T)
        
        
        print(paste(dim(vcf_)[1], ": total number of variant sites before filtering", sep=""))
        write(paste(dim(vcf_)[1], ": total number of variant sites before filtering", sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
        
        # variants with >= 2 alternative alleles were saved to variant_mult_alt.txt
        alt <- grepl(",",vcf_$ALT)
        
        vcf_one_alt <- vcf_[(!alt),]
        vcf_multi_alt <- vcf_[alt,]
        
        # variants with one alternative allele
        vcf_split <- dplyr::select(vcf_one_alt, CHROM:INFO) # create a data frame with columns from CHROM to INFO
        # Extract the MQ and QD from the INFO column
        vcf_split$INFO <- gsub("[a-zA-Z0-9;,\\.\\=\\-]+MQ=(\\d+\\.\\d+)[a-zA-Z0-9;,\\.\\=\\-]+QD=(\\d+\\.\\d+)[a-zA-Z0-9;,\\.\\=\\-]+",
                               "\\1;\\2", vcf_split$INFO)
        vcf_split <- tidyr::separate(vcf_split,col = INFO, into=c("MQ","QD"),sep = ";",remove = T, fill = "right")
        vcf_split$MQ <- as.numeric(vcf_split$MQ); vcf_split$MQ[is.na(vcf_split$MQ)]=0
        vcf_split$QD <- as.numeric(vcf_split$QD); vcf_split$QD[is.na(vcf_split$QD)]=0
        
        #### variants with >= two alternative alleles
        vcf_multi_alt_split <- dplyr::select(vcf_multi_alt, CHROM:INFO) # create a data frame with columns from CHROM to INFO
        if(dim(vcf_multi_alt_split)[1]>0){
            if(max(str_count(vcf_multi_alt_split$ALT,",")) == 1){ 
                n_alt <- 2
            } else if(max(nchar(vcf_multi_alt_split$ALT)) == 2){
                n_alt <- 3
            } else { n_alt <- 3}
            col_ <- c()
            for(i in 1:n_alt){
                col_ <- c(col_, paste("ALT",i,sep=""))
            }
            
            vcf_multi_alt_split <- tidyr::separate(vcf_multi_alt_split, col = ALT, into = col_, sep = ",",remove = T, fill = "right")
            # Extract the MQ and QD from the INFO column
            vcf_multi_alt_split$INFO <- gsub("[a-zA-Z0-9;,\\.\\=\\-]+MQ=(\\d+\\.\\d+)[a-zA-Z0-9;,\\.=-]+QD=(\\d+\\.\\d+)[a-zA-Z0-9;,\\.=-]+",
                                             "\\1;\\2", vcf_multi_alt_split$INFO)
            vcf_multi_alt_split <- tidyr::separate(vcf_multi_alt_split,col = INFO, into=c("MQ","QD"),sep = ";",remove = T, fill = "right")
            vcf_multi_alt_split$MQ <- as.numeric(vcf_multi_alt_split$MQ); vcf_multi_alt_split$MQ[is.na(vcf_multi_alt_split$MQ)]=0
            vcf_multi_alt_split$QD <- as.numeric(vcf_multi_alt_split$QD); vcf_multi_alt_split$QD[is.na(vcf_multi_alt_split$QD)]=0
        }
        ####
        
        # Splitting the sample-specific information (GT:AD:DP:GQ:PL) into 8 columns ("GT","ref","alt","DP","GQ","PL_AA","PL_AB","PL_BB").
        # Also calculate the proportion of reads that support the reference allele.
        for(n_ in names_all){
            # Split one column into eight columns
            vcf_tmp <- reshape2::colsplit(vcf_one_alt[,n_], ",|:",paste(n_, c("GT","ref","alt","DP","GQ","PL_AA","PL_AB","PL_BB"),sep="_"))
            for(i in 2:8) {vcf_tmp[,i] <- as.numeric(vcf_tmp[,i])} # Turn columns into numeric
            # Add a column calculating the frequency of the reference allele
            ref_freq_n <- paste(n_,"_ref_freq",sep="") 
            vcf_tmp[,ref_freq_n] <- vcf_tmp[,2]/vcf_tmp[,4] 
            vcf_split <- cbind(vcf_split,vcf_tmp)
        }
        vcf_split[is.na(vcf_split)] = 0 # turn NA into 0
        
        #### variants with >= two alternative alleles
        if(dim(vcf_multi_alt_split)[1]>0){
            for(n_ in names_all){
                # Split one column into eight columns
                vcf_tmp <- reshape2::colsplit(vcf_multi_alt[,n_], ",|:",paste(n_, c("GT","ref",col_,"DP","GQ","PL_AA","PL_AB","PL_BB"),sep="_"))
                for(i in 2:8) {vcf_tmp[,i] <- as.numeric(vcf_tmp[,i])} # Turn columns into numeric
                # Add a column calculating the frequency of the reference allele
                alt_freq_n <- paste(n_,"_",col_,"_freq",sep="")
                alleles <- grep(paste(c("ref",col_), collapse = "|"), colnames(vcf_tmp), value = T)
                alleles_freq <- vcf_tmp[,alleles]
                alleles_freq$sum <- rowSums(alleles_freq)
                vcf_tmp[,alt_freq_n] <- dplyr::select(alleles_freq/alleles_freq$sum, -contains("ref"), -sum)
                vcf_multi_alt_split <- cbind(vcf_multi_alt_split,vcf_tmp)
            }
            vcf_multi_alt_split[is.na(vcf_multi_alt_split)] = 0 # turn NA into 0
        }
        ####
        
        # FILTER 1: if the genetype of a given variant in the parental strain has missing values ("."), remove the variant
        print(paste(dim(vcf_split)[1], ": number of variant sites that have one alternative allele", sep=""))
        write(paste(dim(vcf_split)[1], ": number of variant sites that have one alternative allele", sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
        
        gt_missing <- "\\."
        vcf_filter1 <- vcf_split
        gt_missing_l <- grepl(pattern = gt_missing, vcf_filter1[,paste(names_par,"_GT",sep="")]) # Logical vector of where the parental genotype is missing
        vcf_filter1 <- vcf_filter1[!gt_missing_l,] # Remove variants in which the genotype of the parental strain are missing
        
        print(paste(dim(vcf_filter1)[1], ": remaining variants after removing sites where the parental strain has zero read coverage (filter 1).", sep=""))
        write(paste(dim(vcf_filter1)[1], ": remaining variants after removing sites where the parental strain has zero read coverage (filter 1).", sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
        
        # FILTER 2 : genotype filter
        # Allowed genotypes : "0/0", "0/1", "1/1". 
        gt_ <- "\\./\\.|0/0|0/1|1/1"
        gt_l <- logical(length = dim(vcf_filter1)[1])
        for(n_ in paste(names_all,"_GT",sep="")){
            gt_tmp <- grepl(pattern = gt_, vcf_filter1[, n_]) 
            gt_l <- gt_l + gt_tmp
        }
        gt_l <- gt_l == length(names_all)
        vcf_filter2 <- vcf_filter1[gt_l,]
        
        print(paste(dim(vcf_filter2)[1], ": remaining variants after filter 2 (genotype filter).", sep=""))
        write(paste(dim(vcf_filter2)[1], ": remaining variants after filter 2 (genotype filter). Allowed genotypes : ., 0/0, 0/1, 1/1", sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
        
        
        # FILTER 3.1 : the proportion of reads support the reference allele at a given variant should be >= ref_thld or <= var_thld.
        ref_freq <- dplyr::select(vcf_filter2, contains("ref_freq")) # columns of ref_freq
        ref_freq_n <- colnames(ref_freq)
        # All strains have freqency >= ref_thld (default=0.8), <= var_thld (default=0.2) or between 0.4~0.6.
        vcf_filter3_1 <- vcf_filter2[rowSums((ref_freq <= var_thld | ref_freq >= ref_thld | (ref_freq>=0.4 & ref_freq<=0.6))) == dim(ref_freq)[2],]
        
        print(paste(dim(vcf_filter3_1)[1], ": remaining variants after filter 3.1 (frequency filter).", sep=""))
        write(paste(dim(vcf_filter3_1)[1], ": remaining variants after filter 3.1 (frequency filter)", sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
        
        #### variants with >= two alternative alleles
        if(dim(vcf_multi_alt)[1]>0){
            alt_freq <- dplyr::select(vcf_multi_alt_split, contains("_freq")) # columns of ref_freq
            alt_freq_n <- colnames(alt_freq)
            # All strains have freqency >= ref_thld (default=0.8) or frequency <= var_thld (default=0.2).
            vcf_multi_alt_split <- vcf_multi_alt_split[rowSums((alt_freq <= var_thld | alt_freq >= ref_thld)) == dim(alt_freq)[2],]
        }
        ####
        
        # FILTER 3.2 : a given variant is considered as positive if the following two conditions are both true
        # 1) the porportion of reads supports the reference allele in the parental and the mutant strain should differ >= freq_diff (default value = 0.8).
        # 2) Read depth of both mutant and parental should be equal or larger than read_depth (default value = 5)
        freq_diff <- as.numeric(freq_diff) ; read_depth = as.numeric(read_depth)
        vcf_filter3_2 <- c()
        ref_freq <- dplyr::select(vcf_filter3_1, contains("ref_freq"))
        ref_freq_par <- dplyr::select(ref_freq, contains(names_par)) # frequency of the reference allele in the parental strain
        depth_ <- dplyr::select(vcf_filter3_1,contains("_DP"))
        depth_par <- dplyr::select(depth_,contains(names_par)) # read depth in the parental strain
        
        for(n_ in names_mut){
            ref_freq_tmp <- dplyr::select(ref_freq, contains(n_)) # frequency of the reference allele in a given mutant strain
            depth_tmp <- dplyr::select(depth_,contains(n_)) # read depth in a given mutant strain
            freq_diff_l <- abs(ref_freq_par - ref_freq_tmp) >= freq_diff & abs(ref_freq_par - ref_freq_tmp)<= 0.6 # the difference in allele frequency between the mutant and parental strains
            read_depth_l <- depth_par >= read_depth & depth_tmp >= read_depth # the read depth threshold
            pass_ <- freq_diff_l & read_depth_l # variants that pass both the frequency difference and read depth threshold
            vcf_filter_tmp <- vcf_filter3_1[pass_,]
            
            if(dim(vcf_filter_tmp)[1] > 0) {
                # Add a column with strain names
                vcf_filter_tmp$strain_id <- n_
                
                # Annotate variants (find if variants are within an ORF)
                variant_gr <- GRanges(vcf_filter_tmp$CHROM, IRanges(start = (vcf_filter_tmp$POS), end = (vcf_filter_tmp$POS)))
                grange_ <- findOverlaps(variant_gr, dd_gr)
                gr_variant_gr <- as.matrix(grange_)[,1]
                gr_dd_gr <- as.matrix(grange_)[,2]
                variant_orf <- cbind(vcf_filter_tmp[gr_variant_gr,], dplyr::select(dd_gene_pos[gr_dd_gr,], ddb_g:description, start:strand)) # Add gene information to the variants in a ORF
                
                # The rest of variants are in intergenic regions
                if(length(grange_)>0) {
                    variant_notorf <- vcf_filter_tmp[-gr_variant_gr,]
                } else {variant_notorf <- vcf_filter_tmp}
                # add NAs to gene information columns of variants in intergenic regions
                if (dim(variant_notorf)[1] > 0) {
                    variant_notorf$description = variant_notorf$genename = variant_notorf$ddb_g = variant_notorf$strand = variant_notorf$end = variant_notorf$start = NA
                    vcf_filter_tmp <- rbind(variant_orf, variant_notorf)
                } else { vcf_filter_tmp <- variant_orf }
                
                vcf_filter3_2 <- rbind(vcf_filter3_2, vcf_filter_tmp)
                
            } else { 
                print(paste("No variant satisfies the filter 3.2 (freq_diff & read_depth filters) in: ", n_, sep="")) 
                write(paste("No variant satisfies the filter 3.2 (freq_diff & read_depth filters) in: ", n_, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)  	
            }
        }
        print(paste(dim(vcf_filter3_2)[1], ": remaining variants after filter 3.2 (freq_diff & read_depth filters).", sep=""))
        write(paste(dim(vcf_filter3_2)[1], ": remaining variants after filter 3.2 (freq_diff & read_depth filters).", sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
        
        #### variants with >= two alternative alleles
        vcf_multiple_alt <- c()
        if(dim(vcf_multi_alt_split)[1]>0){
            
            alt_freq <- dplyr::select(vcf_multi_alt_split, contains("_freq"))
            alt_freq_par <- dplyr::select(alt_freq, contains(names_par)) # frequency of the reference allele in the parental strain
            depth_ <- dplyr::select(vcf_multi_alt_split,contains("_DP"))
            depth_par <- dplyr::select(depth_,contains(names_par)) # read depth in the parental strain
            
            for(n_ in names_mut){
                alt_freq_tmp <- dplyr::select(alt_freq, contains(n_)) # frequency of the reference allele in a given mutant strain
                depth_tmp <- dplyr::select(depth_,contains(n_)) # read depth in a given mutant strain
                freq_diff_l <- abs(as.matrix(alt_freq_par) - as.matrix(alt_freq_tmp)) >= freq_diff # the difference in allele frequency between the mutant and parental strains
                freq_diff_l <- rowSums(freq_diff_l) == length(col_)
                read_depth_l <- depth_par >= read_depth & depth_tmp >= read_depth # the read depth threshold
                pass_ <- freq_diff_l & read_depth_l # variants that pass both the frequency difference and read depth threshold
                vcf_filter_tmp <- vcf_multi_alt_split[pass_,]
                
                if(dim(vcf_filter_tmp)[1] > 0) {
                    # Add a column with strain names
                    vcf_filter_tmp$strain_id <- n_
                    
                    # Annotate variants (find if variants are within an ORF)
                    variant_gr <- GRanges(vcf_filter_tmp$CHROM, IRanges(start = (vcf_filter_tmp$POS), end = (vcf_filter_tmp$POS)))
                    grange_ <- findOverlaps(variant_gr, dd_gr)
                    gr_variant_gr <- as.matrix(grange_)[,1]
                    gr_dd_gr <- as.matrix(grange_)[,2]
                    variant_orf <- cbind(vcf_filter_tmp[gr_variant_gr,], dplyr::select(dd_gene_pos[gr_dd_gr,], ddb_g:description, start:strand)) # Add gene information to the variants in a ORF
                    
                    # The rest of variants are in intergenic regions
                    if(length(grange_)>0) {
                        variant_notorf <- vcf_filter_tmp[-gr_variant_gr,]
                    } else {variant_notorf <- vcf_filter_tmp}
                    # add NAs to gene information columns of variants in intergenic regions
                    if (dim(variant_notorf)[1] > 0) {
                        variant_notorf$description = variant_notorf$genename = variant_notorf$ddb_g = variant_notorf$strand = variant_notorf$end = variant_notorf$start = NA
                        vcf_filter_tmp <- rbind(variant_orf, variant_notorf)
                    } else { vcf_filter_tmp <- variant_orf }
                }
                vcf_multiple_alt <- rbind(vcf_multiple_alt, vcf_filter_tmp)
            }    
            
            if(dim(vcf_multiple_alt)[1]>0){
                vcf_multiple_alt <- dplyr::select(vcf_multiple_alt, CHROM:QD,strain_id:description)
                write.table(vcf_multiple_alt, paste(out_dir,"/","variant_mult_alt.txt",sep=""), sep="\t", quote=F, row.names=F) 
            } else {
                vcf_multiple_alt <- dplyr::select(vcf_multiple_alt, CHROM:QD)
                write.table(vcf_multiple_alt, paste(out_dir,"/","variant_mult_alt.txt",sep=""), sep="\t", quote=F, row.names=F) 
            }
        }
        print(paste(dim(vcf_multiple_alt)[1], ": variants with multiple alternative alleles that pass filter 3.2 (freq_diff & read_depth filters).", sep=""))
        write(paste(dim(vcf_multiple_alt)[1], ": variants with multiple alternative alleles that pass filter 3.2 (freq_diff & read_depth filters).", sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
        
        ####
        
        if(is.null(vcf_filter3_2)){
            write("No variant passes the filter 3.2 (freq_diff & read_depth filters)",paste(out_dir,"/","summary.txt",sep=""), append=T)
            
        } else {
            # ADDING ADDITIONAL INFORMATION ABOUT variants   
            # Add the dna_change column (nucleotide change)
            for( j in seq(dim(vcf_filter3_2)[1])){
                vcf_filter3_2$dna_change[j] <- paste(vcf_filter3_2$REF[j], ">", vcf_filter3_2$ALT[j], sep="")
            }
            
            # Count the number of mutation events on each gene (count one for each of 
            ## the variants on the same gene among all strains)
            gene_list <- data.frame(genename = names(table(as.character(vcf_filter3_2$genename))))
            gene_list$snp_occurrence <- as.integer(table(as.character(vcf_filter3_2$genename)))
            if(dim(gene_list)[1]>0){
                # Count the number of independent mutation sites of each gene  
                var_count_tmp <- sapply(tapply(X = vcf_filter3_2$POS, INDEX = vcf_filter3_2$genename, FUN = unique), length)
                var_count <- data.frame(genename = names(var_count_tmp))
                var_count$var_count <- var_count_tmp
                # The number of independent strains had mutations on each gene.
                strain_count_tmp <- sapply(tapply(X = vcf_filter3_2$strain_id, INDEX = vcf_filter3_2$genename, FUN = unique), length)
                strain_count <- data.frame(genename = names(strain_count_tmp))
                strain_count$strain_count <- strain_count_tmp
                # Count incidence of GA or CT mutations
                GAorTC_tmp = tapply(X = vcf_filter3_2$dna_change, INDEX = vcf_filter3_2$genename, FUN = list)
                GAorTC = data.frame(genename = names(GAorTC_tmp))
                GAorTC$GAorTC = sapply(GAorTC_tmp, function(x) ( length(grep("G>A",x)) + length(grep("A>G",x)) + length(grep("C>T",x)) + length(grep("T>C",x))))
                # Merge to the variant file 
                gene_list_all <- merge(merge(merge(merge(dd_gene_pos, gene_list, by= "genename"), var_count, by = "genename"), strain_count, by = "genename"), GAorTC, by = "genename")
            }    
            # SAVE THE DATA FRAME INTO FILES
            remov_ <- paste(parental_strain, "|", mutant_strain, sep = "")
            remov_index <- grep(remov_, colnames(vcf_filter3_2))
            var_filtered_concise <- vcf_filter3_2[,-remov_index]
            
            
            # ANNOTATING variants (MISENSE, NONSENSE, INTRONIC or INTERGENIC mutations)
            dd_pri_cds <- seqinr::read.fasta(file.path(ref_file, "dicty_primary_cds.txt.gz"), set.attributes = F)
            dd_gff <- read.table(file.path(ref_file, "cds_09-26-2015.gff.gz"), sep="\t", header=T)
            dd_splice_d <- read.table(file.path(ref_file, "splice_donor.gff.gz"), sep="\t", header=T)
            dd_sd <- GRanges(dd_splice_d$chromosome, IRanges(start = dd_splice_d$start, end = dd_splice_d$end))
            dd_splice_a <- read.table(file.path(ref_file, "splice_acceptor.gff.gz"), sep="\t", header=T)
            dd_sa <- GRanges(dd_splice_a$chromosome, IRanges(start = dd_splice_a$start, end = dd_splice_a$end))
            
            
            # Annotating variants (changes in DNA codon and amino acid)
            codon_dna = codon_dna_mut = codon_pro = codon_pro_mut = aa_pos = splice_variant = other_splice_variant =  c()
            
            for(i in seq(dim(var_filtered_concise)[1])){
                
                s_ <- var_filtered_concise[i,]
                s_gr <- GRanges(s_$CHROM, IRanges(start = (s_$POS), end = (s_$POS)))
                
                ov_orf <- findOverlaps(dd_gr, s_gr) # Find if variant lies on a ORF
                if (length(ov_orf) == 0) { # If the variant is not in any ORF
                    codon_dna <- c(codon_dna, "intergenic")
                    codon_dna_mut <- c(codon_dna_mut, "intergenic")
                    codon_pro <- c(codon_pro, "intergenic")
                    codon_pro_mut <- c(codon_pro_mut, "intergenic")
                    aa_pos <- c(aa_pos, NA)
                    splice_variant = c(splice_variant, NA)
                    other_splice_variant = c(other_splice_variant, NA)
                }
                else { # If the variant is in a ORF
                    g_ <- as.character(dd_gene_pos[as.matrix(ov_orf)[1,1],]$ddb_g) # DDB_G gene name of the ORF
                    pri_cds <- dd_pri_cds[grep(g_, names(dd_pri_cds))[1]] # use the first gene model if there are more than one alternative transcripts
                    pri_cds_all <- dd_pri_cds[grep(g_, names(dd_pri_cds))]
                    
                    # a vector contains dictybase ID of other splice variants if there is any 
                    if(length(pri_cds_all) > 1) {
                        tmp_ = c()
                        for(i in 2:length(pri_cds_all)) { 
                            tmp_ <- c(tmp_, strsplit(names(pri_cds_all)[i], split = "\\|")[[1]][1])
                        } 
                        tmp_ <- paste(tmp_, collapse=";")
                        other_splice_variant <- c(other_splice_variant, tmp_)
                    } else if (length(pri_cds_all) == 1) {
                        other_splice_variant <- c(other_splice_variant,  NA)
                    } else { other_splice_variant <- c(other_splice_variant, NA) }
                    
                    # non-coding gene or pseudogene
                    if(pri_cds == "NULL"){
                        codon_dna <- c(codon_dna, "NCG")
                        codon_dna_mut <- c(codon_dna_mut, "NCG")
                        codon_pro <- c(codon_pro, "NCG")
                        codon_pro_mut <- c(codon_pro_mut, "NCG")
                        aa_pos <- c(aa_pos, NA)
                        splice_variant = c(splice_variant, NA)
                    } else {
                        splice_id <- strsplit(names(pri_cds),split = "\\|")[[1]][1] # dictybase ID of the first splice variant
                        cds_ <- dd_gff[grepl(splice_id, dd_gff$attribute) & dd_gff$feature == "CDS",] # find CDS
                        cds_gr <- GRanges(cds_$chr, IRanges(start = cds_$start, end = cds_$end))
                        
                        ov_cds <- findOverlaps(cds_gr, s_gr) # Find if variant lies on a CDS
                        if (length(ov_cds) == 0) { # If the variant is not in a CDS
                            
                            ov_sd <- findOverlaps(dd_sd, s_gr) # splice donor
                            ov_sa <- findOverlaps(dd_sa, s_gr) # splice acceptor
                            
                            if(length(ov_sd) > 0) {
                                codon_dna <- c(codon_dna, "splice_donor")
                                codon_dna_mut <- c(codon_dna_mut, "splice_donor")
                                codon_pro <- c(codon_pro, "splice_donor")
                                codon_pro_mut <- c(codon_pro_mut, "splice_donor")
                                aa_pos <- c(aa_pos, NA)
                                splice_variant = c(splice_variant, splice_id)  
                            } else if (length(ov_sa) > 0){
                                codon_dna <- c(codon_dna, "splice_acceptor")
                                codon_dna_mut <- c(codon_dna_mut, "splice_acceptor")
                                codon_pro <- c(codon_pro, "splice_acceptor")
                                codon_pro_mut <- c(codon_pro_mut, "splice_acceptor")
                                aa_pos <- c(aa_pos, NA)
                                splice_variant = c(splice_variant, splice_id) 
                            } else {
                                codon_dna <- c(codon_dna, "intronic")
                                codon_dna_mut <- c(codon_dna_mut, "intronic")
                                codon_pro <- c(codon_pro, "intronic")
                                codon_pro_mut <- c(codon_pro_mut, "intronic")
                                aa_pos <- c(aa_pos, NA)
                                splice_variant = c(splice_variant, splice_id)
                            }
                        } else { # If the variant is in a CDS
                            splice_variant = c(splice_variant, splice_id)
                            if (unique(cds_$strand) == "-") { # variant is in the negative strand of a CDS
                                cds_ <- cds_[dim(cds_)[1]:1,]
                                cds_pos = c()
                                # get the nucleotide position of CDS
                                for(j in seq(dim(cds_)[1])){ cds_pos = c(cds_pos, seq(cds_[j,]$end, cds_[j,]$start)) }
                            }
                            else { # variant is in the positive strand of a CDS
                                cds_pos = c()
                                # get the nucleotide position of CDS
                                for(j in seq(dim(cds_)[1])){ cds_pos = c(cds_pos, seq(cds_[j,]$start, cds_[j,]$end)) }
                            }
                            
                            # the nuleotide position of tri-nucleotides (codon)
                            codon_idx = seq(1, length(cds_pos),3)
                            codon_pos = paste(cds_pos[codon_idx], cds_pos[(codon_idx+1)], cds_pos[(codon_idx+2)], sep=":")
                            
                            pri_cds_ <- names(dd_pri_cds)[grep(g_, names(dd_pri_cds))] # Find the matched primary CDS sequence
                            if (length(pri_cds_) != 1){
                                pri_cds_ <- pri_cds_[1] # If there are alternative transcripts, use the first one
                            }
                            
                            cds_seq <- dd_pri_cds[[pri_cds_]]
                            stopifnot(length(cds_pos) == length(cds_seq)) # Check if the length of the primary CDS sequence (cds_seq) is the same as the length inferred from start and stop positions of CDSs (cds_pos)
                            
                            g_coding <- c() # the trinucletide codons
                            for(k in codon_idx) { g_coding <- c(g_coding, paste(cds_seq[k:(k+2)], collapse="")) }
                            
                            c_id <- grep(s_$POS, cds_pos) %% 3 ; if(c_id == 0){c_id = 3}
                            #c_id = grep(as.character(s_$POS), strsplit(codon_pos[grep(as.character(s_$POS), codon_pos)],":")[[1]])
                            
                            codon_variant <- g_coding[grep(as.character(s_$POS), codon_pos)]
                            codon_variant <- strsplit(codon_variant,"")[[1]]
                            
                            codon_dna <- c(codon_dna, paste(codon_variant,collapse="")) # reference codon
                            codon_variant_ntg <- codon_variant
                            if (unique(cds_$strand) == "+") {
                                codon_variant_ntg[c_id] <- tolower(as.character(s_$ALT))
                            } else {codon_variant_ntg[c_id] <- tolower(seqinr::comp(as.character(s_$ALT)))}
                            
                            codon_dna_mut = c(codon_dna_mut, paste(codon_variant_ntg, collapse="")) # mutant codon
                            aa_ = grep(as.character(s_$POS), codon_pos) ; aa_pos = c(aa_pos, aa_) # amino acid position
                            
                            codon_variant_translate <- as.character(seqinr::translate(codon_variant))
                            codon_variant_ntg_translate <- as.character(seqinr::translate(codon_variant_ntg))
                            
                            codon_pro <- c(codon_pro, codon_variant_translate)
                            codon_pro_mut <- c(codon_pro_mut, codon_variant_ntg_translate)
                        }
                    }
                }
            }
            
            var_filtered_concise$codon_dna <- codon_dna ; var_filtered_concise$codon_dna_mut = codon_dna_mut
            
            # Add codon usage frequency
            codon_bias <- read.table(file = file.path(ref_file, "dd_codon_bias.txt.gz"), header = T, sep="\t", stringsAsFactors = F)
            codon_freq <- c(); codon_freq_mut <- c()
            for(i in codon_dna){
                c_ <- grep(toupper(i), codon_bias$Codon)
                if(length(c_) > 0){
                    codon_freq <- c(codon_freq, codon_bias$freq[c_])
                } else {codon_freq <- c(codon_freq, "NA")}
                
            } 
            var_filtered_concise$codon_freq <- codon_freq
            for(i in codon_dna_mut){
                c_ <- grep(toupper(i), codon_bias$Codon)
                if(length(c_) > 0){
                    codon_freq_mut <- c(codon_freq_mut, codon_bias$freq[c_])
                } else {codon_freq_mut <- c(codon_freq_mut, "NA")}
            } 
            var_filtered_concise$codon_freq_mut <- codon_freq_mut
            
            
            var_filtered_concise$codon_pro <- codon_pro ; var_filtered_concise$codon_pro_mut = codon_pro_mut
            var_filtered_concise$aa_pos <- aa_pos ; var_filtered_concise$splice_variant <- splice_variant
            var_filtered_concise$other_splice_variants <- other_splice_variant
            
            # Mark synonymous (syn), missense(mis), nonsense(non), intergenic (inter), intronic (intron), splice_donor, splice_acceptor
            for( i in seq(dim(var_filtered_concise)[1])){
                if(var_filtered_concise$codon_pro[i] ==  "intronic" & var_filtered_concise$codon_pro_mut[i] == "intronic"){
                    var_filtered_concise$mutation[i] <- "intronic"
                } 
                else if(var_filtered_concise$codon_pro[i] == "intergenic" & var_filtered_concise$codon_pro_mut[i] == "intergenic") {
                    var_filtered_concise$mutation[i] <- "intergenic"
                }
                else if(var_filtered_concise$codon_pro[i] == "splice_donor" & var_filtered_concise$codon_pro_mut[i] == "splice_donor"){
                    var_filtered_concise$mutation[i] <- "splice_donor"
                }
                else if(var_filtered_concise$codon_pro[i] == "splice_acceptor" & var_filtered_concise$codon_pro_mut[i] == "splice_acceptor"){
                    var_filtered_concise$mutation[i] <- "splice_acceptor"
                }
                else{
                    if(var_filtered_concise$codon_pro[i] == "NCG" & var_filtered_concise$codon_pro_mut[i] == "NCG"){
                        var_filtered_concise$mutation[i] <- "NCG"
                    } 
                    else {
                        if(grepl("[A-Z*]", var_filtered_concise$codon_pro_mut[i])){
                            if(var_filtered_concise$codon_pro_mut[i] == var_filtered_concise$codon_pro[i]){
                                var_filtered_concise$mutation[i] <- "synonymous"
                            } else if(var_filtered_concise$codon_pro_mut[i] != var_filtered_concise$codon_pro[i]){
                                if(var_filtered_concise$codon_pro_mut[i] == "*"){
                                    var_filtered_concise$mutation[i] <- "nonsense"
                                } else { var_filtered_concise$mutation[i] <- "missense"}
                            }
                        }
                        else { print(paste("Cannot recognize the mutation type in row: ", i, sep="")) }
                    }
                }
            }
            
            # Count GAorTC mutation by type (genic or intergenic)
            GAorTC_ <- tapply(X = var_filtered_concise$dna_change, INDEX = var_filtered_concise$mutation, FUN = list)
            
            GAorTC_genic <- sum(grepl("G>A|C>T|A>G|T>C", GAorTC_$missense)) + sum(grepl("G>A|C>T|A>G|T>C", GAorTC_$nonsense)) + 
                sum(grepl("G>A|C>T|A>G|T>C", GAorTC_$synonymous)) + sum(grepl("G>A|C>T|A>G|T>C", GAorTC_$intronic))
            other_genic <- sum(!grepl("G>A|C>T|A>G|T>C", GAorTC_$missense)) + sum(!grepl("G>A|C>T|A>G|T>C", GAorTC_$nonsense)) + 
                sum(!grepl("G>A|C>T|A>G|T>C", GAorTC_$synonymous)) + sum(!grepl("G>A|C>T|A>G|T>C", GAorTC_$intronic))
            GAorTC_intergenic <- sum(grepl("G>A|C>T|A>G|T>C", GAorTC_$intergenic))
            other_intergenic <- sum(!grepl("G>A|C>T|A>G|T>C", GAorTC_$intergenic))
            GAorTC_NCG <- sum(grepl("G>A|C>T|A>G|T>C", GAorTC_$NCG))
            other_NCG <- sum(!grepl("G>A|C>T|A>G|T>C", GAorTC_$NCG))
            GAorTC_total <- sum(grepl("G>A|C>T|A>G|T>C", var_filtered_concise$dna_change))
            other_total <- sum(!grepl("G>A|C>T|A>G|T>C", var_filtered_concise$dna_change))
            
            write("", paste(out_dir, "/", "summary.txt",sep=""), append=T)
            write(paste("The total number of GAorTC mutations in genic region = ", GAorTC_genic, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
            write(paste("The total number of other types of mutations in genic region = ", other_genic, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
            write(paste("The total number of GAorTC mutations in intergenic region = ", GAorTC_intergenic, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
            write(paste("The total number of other types of mutations in intergenic region = ", other_intergenic, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
            write(paste("The total number of GAorTC mutations in non-coding genes = ", GAorTC_NCG, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
            write(paste("The total number of other types of mutations in non-coding genes = ", other_NCG, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
            write(paste("The total number of GAorTC mutations = ", GAorTC_total, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
            write(paste("The total number of other types of mutations = ", other_total, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
            
            # Count the number of each mutation type of all strains
            mis_m <- length(grep("missense", var_filtered_concise$mutation))
            non_m <- length(grep("nonsense", var_filtered_concise$mutation))
            syn_m <- length(grep("synonymous", var_filtered_concise$mutation))
            intron_m <- length(grep("intronic", var_filtered_concise$mutation))
            inter_m <- length(grep("intergenic", var_filtered_concise$mutation))
            NCG_m <- length(grep("NCG", var_filtered_concise$mutation))
            sd_m <- length(grep("splice_donor", var_filtered_concise$mutation))
            sa_m <- length(grep("splice_acceptor", var_filtered_concise$mutation))
            
            write("", paste(out_dir, "/", "summary.txt",sep=""), append=T)
            write(paste("Annotating ",paste("variant_filtered.txt",sep=""), sep=""), paste(out_dir, "/", "summary.txt",sep=""), append=T)
            write(paste("The number of synonymous mutations: ",syn_m, sep=""), paste(out_dir, "/", "summary.txt",sep=""), append=T)
            write(paste("The number of missense mutations: ",mis_m, sep=""), paste(out_dir, "/", "summary.txt",sep=""), append=T)
            write(paste("The number of nonsense mutations: ",non_m, sep=""), paste(out_dir, "/", "summary.txt",sep=""), append=T)
            write(paste("The number of intronic mutations (not including splite donor or acceptor sites): ",intron_m, sep=""), paste(out_dir, "/", "summary.txt",sep=""), append=T)
            write(paste("The number of splice donor site mutations: ",sd_m, sep=""), paste(out_dir, "/", "summary.txt",sep=""), append=T)
            write(paste("The number of splice acceptor site mutations: ",sa_m, sep=""), paste(out_dir, "/", "summary.txt",sep=""), append=T)
            write(paste("The number of intergenic mutations: ",inter_m, sep=""), paste(out_dir, "/", "summary.txt",sep=""), append=T)
            write(paste("The number of mutations in non-coding genes: ",NCG_m, sep=""), paste(out_dir, "/", "summary.txt",sep=""), append=T)
            
            
            
            # Count the number of each mutation type by strains
            mut_list <- tapply(X = var_filtered_concise$strain_id, INDEX = var_filtered_concise$mutation, FUN = list)
            mut_count <- lapply(mut_list, FUN = table)
            mut_by_strain <- data.frame(strain_id = names_mut, stringsAsFactors = F)
            for(i in names(mut_count)){
                df_ <- data.frame(mut_count[i], stringsAsFactors = F); names(df_) <- c("strain_id",i)
                mut_by_strain <- merge(mut_by_strain, df_, by = "strain_id", all = T)
            }
            mut_by_strain[is.na(mut_by_strain)] = 0
            # Check if any mutation type is absesnt, fill it with 0
            mis_ <- setdiff(c('NCG','intergenic','intronic','splice_donor','splice_acceptor','synonymous','missense','nonsense'),names(mut_by_strain))
            if(length(mis_)>=1){
                for(i in mis_){
                    mut_by_strain[,i] <- 0   
                }
            }
            mut_by_strain <- dplyr::select(mut_by_strain, strain_id,intergenic,intronic,splice_donor,splice_acceptor,synonymous,missense,nonsense,NCG)
            mut_by_strain <- dplyr::mutate(mut_by_strain, total = (intergenic+intronic+splice_donor+splice_acceptor+synonymous+missense+nonsense+NCG))
            
            write.table(mut_by_strain, paste(out_dir,"/","mutations_by_strain.txt",sep=""), sep="\t", quote=F, row.names=F) 
            
            # Count the number of variants by chromosome
            chr_list <- tapply(X = var_filtered_concise$CHROM, INDEX = var_filtered_concise$mutation, FUN = list)
            chr_count <- lapply(chr_list, table)
            chr_mut <- data.frame(chromosome = c("chr1","chr2","chr3","chr4","chr5","chr6","chrM","chrR","chrBF","chr2F","chr3F"), stringsAsFactors = F)
            for(i in names(chr_count)){
                df_ <- data.frame(chr_count[i], stringsAsFactors = F); names(df_) <- c("chromosome",i)
                chr_mut <- merge(chr_mut, df_, by = "chromosome", all = T)
                # Check if any mutation type is absesnt, fill it with 0
                mis_ <- setdiff(c('NCG','intergenic','intronic','splice_donor','splice_acceptor','synonymous','missense','nonsense'),names(chr_count))
                if(length(mis_)>=1){
                    for(i in mis_){
                        chr_mut[,i] <- 0   
                    }
                }
                
            }
            chr_mut[is.na(chr_mut)] = 0
            chr_mut <- dplyr::select(chr_mut, chromosome,intergenic,intronic,splice_donor,splice_acceptor,synonymous,missense,nonsense,NCG)
            chr_mut <- dplyr::mutate(chr_mut, exonic = (synonymous+missense+nonsense))
            chr_mut <- dplyr::mutate(chr_mut, total = (NCG+intergenic+intronic+splice_donor+splice_acceptor+synonymous+missense+nonsense))
            
            write.table(chr_mut, paste(out_dir,"/","mutations_by_chr.txt",sep=""), sep="\t", quote=F, row.names=F) 
            
            if(dim(gene_list)[1]>0){
                # List the strains that carry individual mutations
                variant_by_gene = table(var_filtered_concise$genename,var_filtered_concise$strain_id)
                variant_tb <- data.frame(genename=rownames(variant_by_gene), strain_id="NA", stringsAsFactors = F)
                for(i in seq(dim(variant_by_gene)[1]) ){
                    variant_tb$strain_id[i] <- paste(colnames(variant_by_gene)[variant_by_gene[i,]>0], collapse=",")
                }
                write.table(variant_tb, paste(out_dir,"/strain_by_gene.txt",sep=""), row.names = F, quote = F, sep = "\t")
                
                # Count number of nonsense, missense, synonymous mutations of each gene
                mut_list1 <- tapply(var_filtered_concise$genename,var_filtered_concise$mutation,list)
                mut_list2 <- lapply(X = mut_list1, FUN = table)
                mut_by_gene <- data.frame(genename = gene_list_all$genename, stringsAsFactors = F)
                for(i in c("synonymous","missense","nonsense","intronic","splice_donor","splice_acceptor","NCG")){
                    df_ <- data.frame(mut_list2[i], stringsAsFactors = F)
                    if(dim(df_)[1]>0){
                        names(df_) <- c("genename",i)
                        mut_by_gene <- merge(mut_by_gene, df_, by = "genename", all = T)
                    }
                }
                # Check if any mutation type is absesnt, fill it with 0
                mis_ <- setdiff(c("synonymous","missense","nonsense","intronic","splice_donor","splice_acceptor","NCG"),colnames(mut_by_gene))
                if(length(mis_)>=1){
                    for(i in mis_){
                        mut_by_gene[,i] <- 0   
                    }
                }
                mut_by_gene[is.na(mut_by_gene)] = 0
                gene_list_all = merge(gene_list_all,mut_by_gene, by="genename", all = T)
                
                # SAVE THE DATA FRAMES 
                
                write.table(gene_list_all, paste(out_dir,"/","gene_list_all.txt",sep=""), sep="\t", quote=F, row.names=F)
                
                gene_list_top <- dplyr::filter(gene_list_all, snp_occurrence>=2, var_count>=2, strain_count>=2,  GAorTC>=2)
                write.table(gene_list_top, paste(out_dir,"/","gene_list_top.txt",sep=""), sep="\t", quote=F, row.names=F)
            } else { write("No mutation on gene",paste(out_dir,"/","summary.txt",sep=""), append=T)
            }
            write.table(var_filtered_concise, paste(out_dir,"/","variant_filtered.txt",sep=""), sep="\t", quote=F, row.names=F)
            
            #vcf file for IGV viewer
            vcf_one_alt_idx = paste(vcf_one_alt$CHROM,vcf_one_alt$POS, sep=":")
            var_filtered_concise_idx = paste(var_filtered_concise$CHROM,var_filtered_concise$POS, sep=":")
            var_filtered_concise_idx = unique(var_filtered_concise_idx)
            idx = grep(paste(var_filtered_concise_idx,collapse = "|"),vcf_one_alt_idx)
            vcf_one_alt = vcf_one_alt[idx,]; names(vcf_one_alt)[1] = "#CHROM"
            vcf_version <- readLines(input_file,n = 1)
            writeLines(vcf_version,con =  paste(out_dir,"/","variants.vcf",sep=""))
            write.table(vcf_one_alt, paste(out_dir,"/","variants.vcf",sep=""), sep="\t", quote=F, row.names=F,append = T)
            
        }
        EndTime = Sys.time()
        print(EndTime-StartTime) # Running time
        write("", paste(out_dir, "/", "summary.txt",sep=""), append=T)
        write(paste("Start time = ", StartTime, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
        write(paste("End time = ", EndTime, sep=""),paste(out_dir,"/","summary.txt",sep=""), append=T)
        
    }   
}