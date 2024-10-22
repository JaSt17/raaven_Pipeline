GeneCodon <-
  function(seq,organism="ec",max=TRUE,scale=0.5,numcode= 1){
    if(sum(s2c(toupper(seq))%in%c("A","T","G","C")) != nchar(seq)){
      stop("The input sequence must be nucleotide sequence only containg 'ATGC'")
    }
    #data(wSet)
    load("input/wSet.rda")
    assign("translate",seqinr::translate,envir=.GlobalEnv)
    if(seq=="") return("")
    code_invert=hash()
    stardardSeq="GAGTCGTTCAGTAAACACTGTCAACGGCAGGCTCGAATATAACCCGAAACAACGTGGAATTTTCTCACCGTCCTTCCGTGCTCTCTGGGGGTTATGTAGCCAGCAACTGTAAGATACGGCAGCAGGAACTGAAAGTCAGTGATCGCGCGTGGATTACCTCATTATGACCTATTGCGCGATGCCGGTTCCATT"
    for(i in seq(1,192,3)){
      code_invert[[substr(stardardSeq,i,i+2)]]=translate(s2c(stardardSeq)[i:(i+2)],numcode = numcode)
    }
    code=invert(code_invert)
    codon_w <- hash(keys=names(wSet), values=wSet[row.names(wSet) == organism,] )
    amino=translate(s2c(seq),numcode = numcode)
    new_seq=c()
    if(max == TRUE){
      for(i in amino){
        w_value=sapply(code[[i]],function(x)codon_w[[tolower(x)]])
        optimal_codon=code[[i]][sort(w_value,index=TRUE,decreasing=TRUE)$ix[1]]
        new_seq=append(new_seq,optimal_codon)
      }
    } else {
      for(i in amino){
        code_num=length(code[[i]])
        w_value=sapply(code[[i]],function(x)codon_w[[tolower(x)]])
        if(code_num == 1){
          optimal_codon=code[[i]]
        } else {
          optimal_codon=code[[i]][sort(w_value,index=TRUE)$ix[sample(1:round(code_num*scale),1)]]
        }
        new_seq=append(new_seq,optimal_codon)
      }
    }
    return(paste(new_seq,collapse=""))
  }