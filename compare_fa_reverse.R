library('rBLAST')
library('Biostrings')
library('data.table')

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]/3) #not to overload your computer
registerDoParallel(cl)

mus38_chrY <- readDNAStringSet('mus38_split/chrY.fa')
bla <- blast(db='./blastdb/rm_mus39_chrY/rm_mus39_chrY', type = "blastn")

rev_full_matches <- foreach(i=1:length(mus38_chrY)) %dopar% {
  library('rBLAST')
  library('data.table')
  
  match_i = data.table(predict(bla, mus38_chrY[i,], BLAST_args = "-perc_identity 99")) #calling a function
  #do other things if you want
  
  match_i #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}

length(rev_full_matches[lapply(rev_full_matches, nrow) > 0])
rev_full_matches_list <- unique(do.call(c, lapply(rev_full_matches, function (tbl) tbl$'SubjectID')))

mus38_matches_diff <- setdiff(rev_full_matches_list, full_matches_list)

if (identical(length(rev_full_matches),length(mus38_chrY))){
  names(rev_full_matches) <- mus38_chrY@ranges@NAMES
}

if (all(mus38_matches_diff %in% names(rev_full_matches))){
  rev_full_matches[mus38_matches_diff]->tt1
}
library(magrittr)

lapply(tt1, function(i){
  i$SubjectID
})->tt1a

lapply(tt1a,function(i){
  full_matches[i]
})->tt1b

lapply(tt1b,function(i){
  lapply(i,function(j){
    j$SubjectID
  })
})->tt1c

lapply(names(tt1c),function(i){
  lapply(tt1c[[i]],function(j){
    i %in% j
  })
}) %>%
setNames(names(tt1c)) %>%
unlist() %>%
any()





