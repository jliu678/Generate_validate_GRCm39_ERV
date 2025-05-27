library(data.table)

repmask2gtf <- function (name) {
  # browser()
  out_name <- sprintf("chr%1$s/chr%1$s.fa.out", name)
  gff_name <- sprintf("chr%1$s/chr%1$s.fa.out.gff", name)
  gtf_name <- sprintf("gtfs/chr%1$s.gtf", name)
  
  
  header <- c("bit score", "perc div.", "perc del.", "perc ins.", "chr", "begin", "end", "(left)",
              "side", "matching repeat", "repeat class/family", "rep begin", "rep end", "rep (left)",
              "ID", "inc higer score")
  
  data <- read.table(out_name, skip = 2, header = FALSE, fill = TRUE,
                     col.names = paste0("V",seq_len(16)))
  names(data) <- header
  gff <- read.table(gff_name, skip=2, header=F)
  
  gff <- gff[!grepl("Simple|Unknown", data$`repeat class/family`),]
  data <- data[!grepl("Simple|Unknown", data$`repeat class/family`),]
  
  threshold <- 750
  gff <- gff[(data$`bit score` > threshold),]
  data <- data[(data$`bit score` > threshold),]
  
  if(nrow(gff) == 0 || nrow(data) == 0) {
    cat('No Repeatmasker annotation after filtering for',name,'\n')
    return()
  }
  
  # gff <- gff[grepl("ERV", data$`repeat class/family`),]
  # data <- data[grepl("ERV", data$`repeat class/family`),]
  
  seqname <- data$chr
  source <- rep(c("RepeatMasker"), each=length(seqname))
  feature <- rep(c("mRNA"), each=length(seqname))
  start <- data$begin
  end <- data$end
  score <- data$`perc div.`
  strand <- gff$V7
  frame <- rep(c("."), each=length(seqname))
  gene_id <- sapply(
    1:length(seqname), 
    function(x) paste0("gene_id \"RepMasker", x, "\"", collapse='')
    )
  gene_name <- mapply(
    function(i, x, y, z) paste0("gene_name \"", x, ":", y, ":", z, ":", i, "\"", collapse=''),
    1:length(seqname), 
    data$chr, 
    data$`repeat class/family`, 
    data$`matching repeat`
    )
  ID <- sapply(
    1:length(seqname), 
    function(x) paste0("ID \"RepMasker", x, "\"", collapse='')
    )
  attr <- mapply(
    function(id, x, y) paste0(id, "; ", x, "; ", y),
    ID, 
    gene_id, 
    gene_name
    )
  
  gtf <- data.table(
    seqname=seqname,
    source=source,
    feature=feature,
    start=start,
    end=end,
    score=score,
    strand=strand,
    frame=frame,
    attribute=attr
  )
  write.table(gtf, file=gtf_name, quote=FALSE, sep='\t', col.names=F, row.names=F)
}

for (chr in list.dirs(path=".", full.names=TRUE, recursive=FALSE)){
  # browser()
  if (grepl("chr", chr, fixed=TRUE)){
    print(paste0("generating GTF for ", gsub("./chr", "", chr)))
    if (!file.exists(sprintf("gtfs/chr%1$s.gtf", gsub("./chr", "", chr)))){
      repmask2gtf(gsub("./chr", "", chr))
    }else{
      cat(sprintf("gtfs/chr%1$s.gtf", 
                  gsub("./chr", "", chr)), 
          'already exist\n')
    }
  }
}


change_feature_gtf <- function(gtf_tab) {
  gtf_tab$V3 <- rep("exon", nrow(gtf_tab))
  gtf_tab
}

rm_id_gtf <- function(gtf_tab) {
  gtf_tab$V9 <- lapply(gtf_tab$V9, function(x) gsub("(ID \"[A-Za-z0-9]*\"; )", "", x))
  gtf_tab
}

rename_gtf <- function(gtf_tab) {
  gtf_tab$V9 <- lapply(gtf_tab$V9, function(x) {
    new_id <- sprintf("gene_id \"RepMasker_%s\"; ", uuid::UUIDgenerate())
    gsub("(gene_id \"[A-Za-z0-9]*\"; )", new_id, x)
  })
  gtf_tab
}

mod_gtf <- function(chr_name,input_dir="./gtfs/",output_dir="./gtfs2/",skip=0) {
  gtf <- read.table(chr_name, header = FALSE, sep = "\t", quote = "")
  if (!('change_feature_gtf' %in% skip)){gtf <- change_feature_gtf(gtf)}
  if (!('rm_id_gtf' %in% skip)){gtf <- rm_id_gtf(gtf)}
  if (!('rename_gtf' %in% skip)){  gtf <- rename_gtf(gtf)}
  # browser()
  gtf <- as.data.frame(lapply(gtf,unlist))
  write.table(gtf, 
              file=gsub(input_dir, output_dir, chr_name), 
              quote=FALSE, 
              sep='\t', 
              col.names=F, 
              row.names=F)
}

for (chr in list.files(path="./gtfs", full.names=T, recursive=F)) {
  if (!file.exists(gsub("./gtfs/", "./gtfs2/", chr))){
    print(chr)
    mod_gtf(chr)
  }
}

mod_gtf_for_AGAT <- function(chr_name,input_dir="./gtfs/",output_dir="./gtfs2/") {
  gtf <- read.table(chr_name, header = FALSE, sep = "\t", quote = "")
  t_uuid<- stringr::str_extract(string = gtf$V9,
                                pattern = 'gene_id( "RepMasker_.*?"; )',
                                group = 1)
  gtf$V9<-paste0('ID',t_uuid,gtf$V9)
  gtf$V3 <- rep("mRNA", nrow(gtf))
  
  # browser()
  gtf <- as.data.frame(lapply(gtf,unlist))
  write.table(gtf, 
              file=gsub(input_dir, output_dir, chr_name), 
              quote=FALSE, 
              sep='\t', 
              col.names=F, 
              row.names=F)
}

# in order to keep same UUID,convert from gtfs2
dir.create("./gtfs_for_fasta_needed_by_salmon/")
for (chr in list.files(path="./gtfs2", full.names=T, recursive=F)) {
  if (!file.exists(gsub("./gtfs2/", "./gtfs_for_fasta_needed_by_salmon/", chr))){
    print(chr)
    mod_gtf_for_AGAT(chr,input_dir ="./gtfs2/",output_dir="./gtfs_for_fasta_needed_by_salmon/")
  }
}

# **next step is bash below, see ~/siyi_summer2023/erv/01compare_rm/agat_gtf2fa.sh**
# cd /PHShome/jn22/siyi_summer2023/erv/00repeatmasker
# CHRS=`ls gtfs_for_fasta_needed_by_salmon/*.gtf`
# for chr in $CHRS
# do
# echo "extracting .fa from ${chr}"
# echo "output file is ${chr%.*}.fa"
# agat_sp_extract_sequences.pl -g ${chr} -t mRNA --merge 0\
# -f ./GRCm39.primary_assembly.genome.fa -o ${chr%.*}.fa
# done

# cat *.fa > mouse_erv_combined_repeatMasker_on_m39
# then rename the mouse_erv_combined_repeatMasker_on_m39 by adding .fa