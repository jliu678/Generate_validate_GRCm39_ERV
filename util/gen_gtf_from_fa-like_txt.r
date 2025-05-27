data.table::fread('transgene.fa',header = F)->sq

sq[grepl(t(sq),pattern = ' ')] -> tgNames

paste0(t(sq),collapse = '')->sq1
strsplit(sq1,split = paste0(tgNames$V1,collapse = '|'))->seqs
tail(unlist(seqs), -1)->seqs

sapply(tgNames$V1, function(i) gsub("> ", "", i))->seqname
rep(c("custom_"), each=length(seqname))->source
rep(c("exon"), each=length(seqname))->feature
rep(c(1), each=length(seqname))->start
sapply(seqs, function(i) nchar(i))->end
rep(c("500"), each=length(seqname))->score
rep(c("+"), each=length(seqname))->strand # I think they are + right?
rep(c("."), each=length(seqname))->frame

sapply(1:length(seqname), function(x) paste0("gene_id \"custom_", x, "\"", collapse=''))->gene_id
sapply(seqname, function(x) paste0("gene_name \"", x, "\"", collapse=''))->gene_name

mapply(function(x, y) paste0(x, "; ", y), gene_id, gene_name) -> attr

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

# don't forget to remove the row names
write.table(gtf, file='transgene.gtf', quote=FALSE, sep='\t', col.names = F, row.names = F)
