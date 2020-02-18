library("Biostrings")
library("seqinr") 



## Load DNA sequences

fastaFile  = readDNAStringSet("E:/Fast/Maam Zoya/data and code/LncRNAsSmilarity/lncRNA_seq_clean.fa")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df_q <- data.frame(seq_name, sequence)
write.csv(df_q,'E:/Fast/Maam Zoya/data and code/LncRNAsSmilarity/lncRNA_seq_Query.csv')
# Query sequnces
#qry_file <- system.file("E:/Fast/Maam Zoya/data and code/LncRNAsSmilarity/lncRNA_seq_clean.fa", package = "seqinr")
#Qry_Seq <- read.fasta("E:/Fast/Maam Zoya/data and code/LncRNAsSmilarity/query_sample.fa" , seqtype="DNA", as.string="FALSE")
#Qry_Seq1 <- getSequence(Qry_Seq[[1]], as.string=TRUE)

#Target sequences 
#trg_Seq <- read.fasta("E:/Fast/Maam Zoya/data and code/LncRNAsSmilarity/target_sample.fa", seqtype="DNA", as.string="FALSE")
#trg_Seq1 <- getSequence(trg_Seq[[1]], as.string=TRUE)
fastaFile2  = readDNAStringSet("E:/Fast/Maam Zoya/data and code/LncRNAsSmilarity/NONCODEv5_human.fa")
seq_name2 = names(fastaFile2)
sequence2 = paste(fastaFile2)
df_t <- data.frame(seq_name2, sequence2)
write.csv(df_t,'E:/Fast/Maam Zoya/data and code/LncRNAsSmilarity/lncRNA_seq_target.csv')


data("BLOSUM62")


Smith_Waterman_Scores <- data.frame(Query_Seq_ID=as.character(),
                                    Query_Seq=character(), 
                                    Target_Seq_ID=as.character(),
                                    Target_Seq=character(),
                                    Score=as.numeric(),
                                    stringsAsFactors=FALSE) 
length(trg_Seq)
## Perform alignment
for (q in 1:length(Qry_Seq)) {
  # query seq
  #print(typeof(getSequence(Qry_Seq[[q]], as.string =FALSE)))
  qry <- getSequence(Qry_Seq[[q]], as.string =FALSE)
  query_seq<-DNAString(paste(qry,collapse=""))
  # query seq id
  query_seq_id<-getName(Qry_Seq[[q]])
  print("query Id:")
  print(q)
  for (t in 1:length(trg_Seq))
  {
    
    #print(query_seq)
    # target seq
    trg <- getSequence(trg_Seq[[t]], as.string =FALSE)
    target_seq<-DNAString(paste(trg,collapse=""))
    #print(target_seq)
    
    #target seq id
    target_seq_id<-getName(trg_Seq[[t]])
    
    pa <- pairwiseAlignment(query_seq, target_seq,substitutionMatrix=BLOSUM62,gapOpening = -10, gapExtension = -0.5,type="local",scoreOnly = TRUE)
    print("score:")
    print(pa)
    
    Smith_Waterman_Scores <- rbind(Smith_Waterman_Scores, c(query_seq_id,paste(qry,collapse=""),query_seq_id,paste(trg,collapse=""),pa))
    print("traget id:")
    print(t)
    
  }
 
   
  
}

#names(Smith_Waterman_Scores)[1] <- "First.Protein"
#names(Smith_Waterman_Scores)[2] <- "Second.Protein"
#names(Smith_Waterman_Scores)[3] <- "Score"
write.csv(Smith_Waterman_Scores,'E:\\MyData.csv')


### Normalize Smith-Waterman similarity scores
dt <- data.table(Smith_Waterman_Scores)
dt.lookup <- dt[First.Protein == Second.Protein]
setkey(dt,"First.Protein" )
setkey(dt.lookup,"First.Protein" )
colnames(dt.lookup) <- c("First.Protein","Second.Protein","Score1")
dt <- dt[dt.lookup]
setkey(dt,"Second.Protein" )
setkey(dt.lookup,"Second.Protein")
colnames(dt.lookup) <- c("First.Protein","Second.Protein","Score2")
dt <- dt[dt.lookup][
  , Normalized :=  Score / (sqrt(Score1) * sqrt(Score2))][
    , .(First.Protein, Second.Protein, Normalized)]
dt <- dt[order(dt$First.Protein),]

Smith_Waterman_Scores <- as.data.frame(dt)