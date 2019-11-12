install.packages("Peptides")
library(Peptides)

cox1 <- read.fasta(file="R/win-library/3.6/seqinr/sequences/cox1.fasta", seqtype="AA")
aaComp(cox1[1]) #Compute the amino acid composition of a protein sequence for human sequence
aIndex(cox1[1]) #Computes de aliphatic index of the protein sequence
charge(cox1, pH=7, pKscale="Lehninger") #Predict the net charge of the protein
charge(seq="FLPVLAG", pH=7, pKscale="EMBOSS")
hydrophobicity(cox1) #Compute hydrophobicity

#A function for calculating hydrophobicity by chunks
slidingwindowhyrophobicity <- function(windowsize, inputseq)
{
  hpwindow <- seq(1, length(inputseq)-windowsize, by = windowsize) 
  # Find the length of the GCwindow 
  n <- length(hpwindow) 
  # Make a vector the same length but "blank" for us to fill 
  Chunks <- numeric(n) 
  for (i in 1:n) {
    chunk <- inputseq[hpwindow[i]:(hpwindow[i]+windowsize-1)] 
    chunkhp <- hydrophobicity(chunk) 
    print(chunkhp) 
    Chunks[i] <- chunkhp
  }
  plot(hpwindow,Chunks,type="b",xlab="aa start position",ylab="Hydrophobicity",main=paste("Hydrophobcity Plot with windowsize ", windowsize)) 
}

CloningProtein <- cox1[1][[1]]
slidingwindowhyrophobicity(20, CloningProtein)
