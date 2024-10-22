AAtoDNA <- function(inAA,species="hsa",fullOPT=FALSE,optIt=1){
  humanCodon <- read.table(header = TRUE, 
                           stringsAsFactors = FALSE, 
text="AA DNA
A   gcc
C	tgc
D	gac
E	gag
F	ttc
G	ggc
H	cac
I	atc
K	aag
L	ctg
M	atg
N	aac
P	ccc
Q	cag
R	cgg
S	agc
T	acc
V	gtg
W	tgg
Y	tac")

inAA <- toupper(as.character(inAA))
outDNA <- toupper(sedit(inAA,humanCodon$AA, humanCodon$DNA, wild.literal=FALSE))
outDNA <- DNAString(outDNA)
outDNA <- GeneCodon(as.character(outDNA) , organism = species)
}