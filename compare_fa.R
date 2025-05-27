# compare the outputs repeatmasker to that of gEVE

library('rBLAST')
library('Biostrings')
library('data.table')

library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]/2) #not to overload your computer
registerDoParallel(cl)

reptmask_chrY <- readDNAStringSet('repeatmasker/chrY.fa')
# reptmask_chrY_testing <- readDNAStringSet('repeatmasker/chrY_testing.fa')

bla <- blast(db='./blastdb/mus38_chrY/mus38_chrY', type = "blastn")

full_matches <- foreach(i=1:length(reptmask_chrY)) %dopar% {
  library('rBLAST')
  library('data.table')
  
  match_i = data.table(predict(bla, reptmask_chrY[i,], BLAST_args = "-perc_identity 99")) #calling a function
  #do other things if you want
  
  match_i #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
saveRDS(full_matches, file = "repeatmasker_mus38_full_named_matches.rds")


names(full_matches)<-reptmask_chrY@ranges@NAMES %>%
  gsub(pattern=' .*$',replacement='')


#stop cluster
stopCluster(cl)


gtflist <- read.table("repeatmasker/chrY.gtf")

erv_matches <- full_matches[grepl("ERV", gtflist$V16)]
length(erv_matches[lapply(erv_matches, nrow) > 0]) / length(erv_matches)
erv_matches_list <- unique(do.call(c, lapply(erv_matches, function (tbl) tbl$'SubjectID')))
length(erv_matches_list) / length(erv_matches)

line_matches <- full_matches[grepl("LINE", gtflist$V16)]
length(line_matches[lapply(line_matches, nrow) > 0]) / length(line_matches)
line_matches_list <- unique(do.call(c, lapply(line_matches, function (tbl) tbl$'SubjectID')))
length(line_matches_list) / length(line_matches)

other_matches <- full_matches[!grepl("ERV|LINE", gtflist$V16)]
length(other_matches[lapply(other_matches, nrow) > 0]) / length(other_matches)
other_matches_list <- unique(do.call(c, lapply(other_matches, function (tbl) tbl$'SubjectID')))
length(other_matches_list) / length(other_matches)

length(full_matches[lapply(full_matches, nrow) > 0])
full_matches_list <- unique(do.call(c, lapply(full_matches, function (tbl) tbl$'SubjectID')))

# list(predict(bla, reptmask_chrY[2,], BLAST_args = "-perc_identity 97"))
# cl <- list()
# for (i in 40:440){
#   cl <- append(cl,  data.table(predict(bl, reptmask_chrY[i,], BLAST_args = "-perc_identity 99")))
# }
# length(cl[lapply(lapply(cl, data.table), nrow) > 0])
