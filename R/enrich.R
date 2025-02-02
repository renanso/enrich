#' enRich
#' This is the main function of the package.
#' @examples
#' # example code
#' #input values
#' ref= system.file("extdata", "genome.fasta", package = "enRich")# reference genome
#' file.copy(ref,"./") #saving reference genome in the directory
#' index=system.file("extdata", "genome.fasta.fai", package = "enRich")# reference genome index
#' file.copy(index,"./") #saving index in the current directory
#' bed_file= system.file("extdata", "bed_test.bed", package = "enRich") # bed file
#' file.copy(bed_file,"./") #savind bed file in the current directory
#' blast_db= "./db/genome" #blast database name
#' makeblastdb_path= "/path/to/blast/ncbi-blast-2.1X.0+/bin/makeblastdb" 
#' blastn_path= "/path/to/blast/ncbi-blast-2.1X.0+/bin/blastn" 
#' blast_out="./blast_result.txt"
#' frag_size=400 # frament size in bp for the region with the target SNP in the center
#' w_size= 120 #size of the probe
#' s_size= 1 #number of oligonucleotides to move
#' gc_min= 40 # minimum GC%
#' gc_max= 46 # maximum GC%
#' blast_db= "./db/genome" #blast database name
#' pid= 90 # percentage identity of blast result (default 90)
#' len= 120 # length of the match in the blast result (default 120)
#' blast_hit= 1 # maximum number of blast match for the probe (default 1)
#' tm_min= 70 # minimum Tm (default 70)
#' tm_max= 76 #maximum Tm (default 76)
#' HpTm= 60 #Hairpin maximum Tm (default 60)
#' HmTm= 60 # Heterodimer maximum Tm (default 60)
#' # Run the function
#' library(enRich)
#' enrich("genome.fasta", "bed_test.bed", frag_size, w_size, s_size, gc_min, gc_max, blast_db, makeblastdb_path, blastn_path, blast_out, pid, len, blast_hit,tm_min, tm_max, HpTm, HmTm)
#' @importFrom dplyr %>%
#' @export
enrich <- function(ref, bed_file, frag_size, w_size, s_size, gc_min, gc_max, 
                   blast_db, makeblastdb_path, blastn_path, blast_out, pid, len, blast_hit,tm_min, tm_max, HpTm, HmTm)

{
  message("Reading genome file")
  ## Read Genome file
  genome<- Biostrings::readDNAStringSet(ref, format="fasta",
                            nrec=-1L, skip=0L, seek.first.rec=FALSE,
                            use.names=TRUE, with.qualities=FALSE)
  seq_name <- names(genome)
  sequence <- paste(genome)
  df <- data.frame(seq_name, sequence)
  message("Done")
  colnames(df)<-c('qseqid','sequence')
  
  ## Read Bed file
  message("Reading bed file")
  bed<-read.table(bed_file)
  colnames(bed)<-c("chr","start","end","locus")
  message("Done")
  
  ## Cut targets
  message("Cutting fragments")
  bed$strand <- rep("+", nrow(bed))
  ranges <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns = T, start.field="start", end.field="end", seqnames.field = "chr")
  targets<- Biostrings::getSeq(genome, ranges)
  seq_name<- bed$locus
  sequence<- paste(targets)
  df1 <- data.frame(seq_name, sequence)
  colnames(df1)<-c('qseqid','sequence')
  write.csv(df1, "targets_seq.csv")
  message("Done")
  
  ## Slide window
  message("Running sliding window")
  df2<-df1 %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Windows = list(slide_window(sequence, w_size, s_size))) %>%
    tidyr::unnest(Windows) %>%
    dplyr::group_by(sequence, qseqid) %>%
    dplyr::mutate(Window_ID = paste(qseqid, row_number(), sep = "_")) %>%
    dplyr::ungroup() %>%
    dplyr::select(qseqid, Windows)
  names(df2$Windows)<-df2$qseqid
  df3<-stack(df2$Windows)
  
  ## add unique number to each probe
  df3$ind2<-paste0(df3$ind,"_",row.names(df3))
  df3<-df3[,c(1,3)]
  message("Done")
  
  ## GC filter
  message("Calculating GC%")
  colnames(df3)<-c('sequence', 'qseqid')
  df3$gc<-round(as.numeric(sapply(df3[,1],gc_fun)),1)
  message("Done")
  
  ##Apply filters
  message("GC% filter")
  df4<- df3 %>% dplyr::filter(gc >= gc_min, gc <= gc_max)
  message("Done")
  
  #plot gc content before and after
  pdf("gc_content.pdf")
  par(mfrow=c(1,2))
  hist(df3$gc, main = "Original GC %", breaks = 5)
  hist(df4$gc, main = "Filtered GC %", breaks = 5)
  dev.off()

  #probe candidates number
  message(paste("Number of GC filtered candidates:",nrow(df4)))
  
  ##changing table to fasta format
  df4[,2] <- paste0(">",df4[,2]) #add ">" to headers
  
  #bind rows of headers and seqs
  probes_fasta <- c(rbind(df4[,2], df4[,1]))
  probes_fasta[1:10]
  write.table(probes_fasta, "candidates_gc_filtered.fasta", row.names=FALSE,sep="\t", quote = FALSE, col.names = FALSE)
  
  ##BLAST
  
  message("Creating Blast database")
  
  # check makeblastdb installation
  Sys.which("makeblastdb")
  Sys.which("blastn")
  
  # Construct the makeblastdb command
  command1 <- paste0(makeblastdb_path,
                     " -in ", ref,
                     " -dbtype ", "nucl",
                     " -out ", blast_db,
                     " -parse_seqids ")
  # Run the makeblastdb command
  system(command1)
  message("Done")
  
  message("Running Blast")
  # Construct the BLASTN command
  command2 <- paste0(blastn_path,
                     " -db ", blast_db,
                     " -query ", "./candidates_gc_filtered.fasta",
                     " -out ", blast_out,
                     " -outfmt 6")
  # Run the BLASTN command
  system(command2)
  
  #add names to columns
  blast_result<-read.table("blast_result.txt")
  blast_names<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
  colnames(blast_result)<- blast_names
  
  message("Done") 
  
  message(paste("Size of BLAST object:",nrow(blast_result)))
  
  message("Filtering length and identity in the Blast result")
  
  ## Filter BLAST, Tm, hairpin and homodimer calculations
  ## Filter 1: % of identity and fragment length
  #probes binding in the correct place and with high identity to off targets will pass.
  
  probes1<- blast_result %>% dplyr::filter(length == len) %>% dplyr::filter(pident > pid)
  message(paste("Candidate probes here:",nrow(probes1)))
  
  ## check unique genes before next filter
  gene<-(stringr::str_split_fixed(probes1$qseqid, '\\|', 3))[,2]
  probes1$gene<-gene
  probes1$gene<-as.factor(probes1$gene)
  summary(probes1$gene)
  message("Done")
  
  ## Filter 2: Filter out probes that binds to other parts of the genome based on the blast_hit parameter
  message("Filtering probe candidates that binds to multiple sites in the genome")
  probes2<-probes1 %>%
    dplyr::group_by(qseqid) %>%
    dplyr::filter(n()<=blast_hit)
  summary(probes2$gene)
  message("Done")
  message(paste("Candidate probes here:",nrow(probes2)))
  
  ## Filter 3: Matching pattern to filter probes with specific match to target region. these filter are used to remove potential off-targets
  ## chromosome and position match test
  original_snp<-(stringr::str_split_fixed(probes2$qseqid, "\\|", 3))[,1]
  original_chr<-(stringr::str_split_fixed(original_snp, "_", 2))[,1]
  original_position<-as.numeric((stringr::str_split_fixed(original_snp, "_", 2))[,2])
  
  ## comparing chromosome from original SNP name and the match in the Blast search
  message("Checking if probe is in the right chromosome and right position")
  probes2$original_chr<- original_chr
  probes2$target_chr<-probes2$sseqid
  probes2$chr_match<- probes2$original_chr==probes2$target_chr
  
  ## comparing position from original SNP name and the match in the Blast search
  probes2$original_position<- original_position
  probes2$test<- probes2$original_position - probes2$sstart
  probes2$pos_match<- probes2$test < (frag_size/2) & probes2$test > -(frag_size/2) ##assuming the maximum distance base on the design 200bp up and down
  
  ## filter chromosome and positions
  probes3<-probes2 %>% dplyr::filter(chr_match == "TRUE") %>% dplyr::filter(pos_match == "TRUE")
  message(paste("Candidate probes here:",nrow(probes3)))
  
  ## drop the initial datasets to free up memory
  message("Done")
  rm(blast_result, probes1)
  gc()
  
  ## Thermodynamics
  #prepare data set adding sequence to blast results filtered
  dna <- Biostrings::readDNAStringSet('candidates_gc_filtered.fasta', format='fasta')
  seq_name = names(dna)
  sequence = paste(dna)
  df5 <- data.frame(seq_name, sequence)
  colnames(df5)<-c('qseqid','sequence')
  probes5 <- (merge(df5, probes3, by = 'qseqid'))
  
  ##Tm-iterate on the column (apply function to each row in a column)
  probes5$tm<-round(as.numeric(sapply(probes5[,2],calculate_tm)),1)
  
  # probes were already filtered for GC, but we will add the information as a column
  probes5$gc<-round(as.numeric(sapply(probes5[,2],gc_fun)),1)
  
  ## Apply Tm filters
  message("Applying filter for Tm")
  probes6<- probes5 %>% dplyr::filter(tm >= tm_min, tm <= tm_max)
  message("Done")
  message(paste("Candidate probes here:",nrow(probes6)))
  
  ## Histograms
  pdf("TM_gc_after_filter.pdf")
  par(mfrow=c(1,2))
  hist(probes6$gc, main = "GC content filtered")
  hist(probes6$tm, main = "Tm filtered")
  dev.off()
  
  ##secondary structure
  ##Note that the maximum length of ``seq`` is 60 bp. This is a cap suggested
  #by the Primer3 team as the longest reasonable sequence length for which a
  #two-state NN model produces reliable results.
  # We can still analyze the probes by creating fragments up to 60bp.
  #frag1 = first 60bp (1 to 60)
  #frag2 = last 60bp (61 to 120)
  #frag3 = first 30bp + last 30 bp (1 to 30 + 91 to 120)
  #frag4 = middle 60bp (31 to 90)
  #frag5 = first 30bp + 60-90bp (1 to 30 + 61 to 90)
  #frag6 = 30-60 + last 30 bp (31 to 60 + 91 to 120)
  probes6$frag1<-substr(probes6$sequence, 1, 60)
  probes6$frag2<-substr(probes6$sequence, 61, 120)
  probes6$frag3<-paste0(substr(probes6$sequence, 1, 30), substr(probes6$sequence, 91, 120))
  probes6$frag4<-substr(probes6$sequence, 31, 90)
  probes6$frag5<-paste0(substr(probes6$sequence, 1, 30), substr(probes6$sequence, 61, 90))
  probes6$frag6<-paste0(substr(probes6$sequence, 31, 60), substr(probes6$sequence, 91, 120))
  
  #fragment list
  fragments<- list(probes6$frag1,probes6$frag2,probes6$frag3,probes6$frag4,probes6$frag5,probes6$frag6)
  
  ##hairpin Tm calculation
  message("Hairpin Tm calculation")
  probes6_hairpin<-lapply(fragments,calculate_hairpin)
  probes6$hairpin_temp_frag1<-probes6_hairpin[[1]][["temp"]]
  probes6$hairpin_temp_frag2<-probes6_hairpin[[2]][["temp"]]
  probes6$hairpin_temp_frag3<-probes6_hairpin[[3]][["temp"]]
  probes6$hairpin_temp_frag4<-probes6_hairpin[[4]][["temp"]]
  probes6$hairpin_temp_frag5<-probes6_hairpin[[5]][["temp"]]
  probes6$hairpin_temp_frag6<-probes6_hairpin[[6]][["temp"]]
  message("Done")
  
  ## homodimer Tm calculation
  message("Homodimer Tm calculation")
  probes6_homodimer<-lapply(fragments,calculate_homodimer)
  probes6$homodimer_temp_frag1<-probes6_homodimer[[1]][["temp"]]
  probes6$homodimer_temp_frag2<-probes6_homodimer[[2]][["temp"]]
  probes6$homodimer_temp_frag3<-probes6_homodimer[[3]][["temp"]]
  probes6$homodimer_temp_frag4<-probes6_homodimer[[4]][["temp"]]
  probes6$homodimer_temp_frag5<-probes6_homodimer[[5]][["temp"]]
  probes6$homodimer_temp_frag6<-probes6_homodimer[[6]][["temp"]]
  message("Done")
  message("Hairpin and homodimer filters")
  
  ## filter for hairpin Tm
  probes7<- probes6 %>% dplyr::filter(hairpin_temp_frag1 < HpTm) %>%
    dplyr::filter(hairpin_temp_frag2 < HpTm) %>%
    dplyr::filter(hairpin_temp_frag3 < HpTm) %>%
    dplyr::filter(hairpin_temp_frag4 < HpTm) %>%
    dplyr::filter(hairpin_temp_frag5 < HpTm) %>%
    dplyr::filter(hairpin_temp_frag6 < HpTm)
  
  ## filter for homodimer Tm
  probes8<- probes7 %>% dplyr::filter(homodimer_temp_frag1 < HmTm) %>%
    dplyr::filter(homodimer_temp_frag2 < HmTm) %>%
    dplyr::filter(homodimer_temp_frag3 < HmTm) %>%
    dplyr::filter(homodimer_temp_frag4 < HmTm) %>%
    dplyr::filter(homodimer_temp_frag5 < HmTm) %>%
    dplyr::filter(homodimer_temp_frag6 < HmTm)
  message("Done")
  message("Sort probes by GC%")
  
  ## sort probes from highest GC
  probes9 <- probes8[order(probes8$gc, decreasing = TRUE),]
  message("Done")
  message("Keeping one probe per gene")
  
  ## Filter for one probe per gene
  probes10<- probes9 %>% dplyr::distinct(gene, .keep_all = TRUE)
  message("Done")
  
  #Visualize probe filters
  message("Plot probe filter")
  ## plot probe filtering
  vals<- c(nrow(df3), nrow(df4), nrow(probes2), nrow(probes3), nrow(probes6), nrow(probes9), nrow(probes10))
  filters<-c("1-all_candidates","2-gc_test","3-blast_test","4-position_test","5-tm_test","6-structure_test","7-unique_test")
  flow_data<- data.frame(cbind(filters,vals))
  flow_data$vals<-as.numeric(flow_data$vals)
  plot2<-ggplot(flow_data, aes(y=vals, x=filters)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=vals), vjust=-0.5, size = 5) +
  ggtitle("Probe Summary") +
  ylab("probes (log10)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20)) +
  scale_y_log10()
  pdf("plot_filters.pdf", width = 15, height = 10)
  print(plot2) 
  dev.off() 
  message("Done")
  
  ## Visualize probe positions
  message("Plot probe position")
  ##bins every 1M
  dist=1e6
  bins<-seq(from=1, to=0.8e11, by=dist)
  ##adding the bins to dataframe
  probes10$bins<-findInterval(probes10$sstart, bins)
  # extracting the info needed for plots
  df6<-probes10[,c(3,10,35)]
  df6$bins2<-paste0(df6$sseqid,"-",df6$bins)
  df7<-data.frame(table(df6$bins2))
  colnames(df7)<-c("bins2","freq")
  df8 <- merge(df6,df7,by="bins2")
  #keep unique bins
  df9<- df8 %>% dplyr::distinct(bins2, .keep_all = TRUE)
  df9$original_chr <- as.numeric(gsub('Chr', '', df9$sseqid))
  n_chr<- nrow(df)
  chr_len<- nchar(df$sequence)
  plot<- ggplot(df9,aes(x=sstart, y= original_chr)) +
         ggplot2::geom_count(aes(size=freq)) +
         ggplot2::scale_y_continuous(breaks = seq(1, n_chr, by = 1)) +
         ggplot2::theme_bw() +
         ggplot2::theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                         panel.border = element_blank(),panel.background = element_blank()) +
         ggplot2::annotate("segment", x = 1, xend = c(chr_len), y = seq(1, n_chr, by = 1),
                           yend = seq(1, n_chr, by = 1),colour = "black") +
         ggplot2::scale_size_continuous(name = paste("Probe density \nBin size=", dist, "bp")) +
         ggplot2::xlab("Position (bp)") + ggplot2::ylab("Chromosome") +
         ggplot2::theme(legend.position.inside =c(0.92, 0.92),legend.justification = c("right", "top"),
                        legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))
  pdf("plot_all_probes.pdf", width = 15, height = 10)
  print(plot)
  dev.off()
  message("Done")
  
  ## Formating the final files
  message("Saving final files")
  ## probes summary
  write.csv(probes10, "final_probes_summary.csv")
  ## fasta file
  names<-(stringr::str_split_fixed(probes10$qseqid, '\\|',3))[,c(1:3)]
  final_probes<- cbind(paste0(names[,1],"-",names[,2],"-",names[,3],"-Tm=", probes10$tm,"-GC=", probes10$gc), probes10$sequence)
  ## changing table to fasta format
  final_probes[,1] <- paste0(">",final_probes[,1]) #add ">" to headers
  ## bind rows of headers and seqs
  probes_fasta <- c(rbind(final_probes[,1], final_probes[,2]))
  probes_fasta[1:10]
  write.table(probes_fasta, "final_probes.fasta", row.names=FALSE,sep="\t", quote = FALSE, col.names = FALSE)
  message("Done")
  message("enRich is finalized")
  proc.time()
  
}
