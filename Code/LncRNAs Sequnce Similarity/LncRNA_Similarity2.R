SeqScore <- function()
{

library("Biostrings")
library("seqinr") 
data("BLOSUM100")
query_seq  <- readDNAStringSet("E:/Fast/Maam Zoya/data and code/LncRNAsSmilarity/query_sample.fa")

target_seq  <- readDNAStringSet("E:/Fast/Maam Zoya/data and code/LncRNAsSmilarity/target_sample.fa")

#Create empty list to fill with percent scores and matching sequences:
DList=NULL
QueSeqList = NULL
TotList = NULL

#initiating the counters:
i=1
j=1
c=0

#Perform alignment and generate percent identity scores:
while(i<=length(query_seq))
{
  while(j<=length(target_seq))
  {
    SeqAlign <- pairwiseAlignment(query_seq, target_seq[j],substitutionMatrix=BLOSUM100,gapOpening = -10, gapExtension = -0.5,type="local",scoreOnly = TRUE)
    PercAlign <- SeqAlign
    
      DList = c(DList, names(query_seq[i]), names(target_seq[j]))
      QueSeqList = c(QueSeqList, toString(target_seq[j]))
      c=c+1
    
  }
  i=i+1
  j=1 #to reset the inner while loop
}


unlist<-t(sapply(DList, unlist));
outputMatrix<-cbind(DList,QueSeqList)
outputMatrix<-as.matrix(outputMatrix, ncol=3)
write.csv(outputMatrix, "E:/Fast/Maam Zoya/data and code/LncRNAsSmilarity/outputMatrix.csv")
}