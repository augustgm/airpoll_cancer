library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg19)
findReverseComplement = function(n_string){
  
  rc = c("A","T","C","G")
  names(rc) = c("T","A","G","C")
  
  n_string_split = unlist(strsplit(n_string,""))
  paste(rc[n_string_split],collapse="")
  
}

getGenomeContext = function(bs_genome){
  tri_context <- matrix(ncol = 64, nrow = length(BSgenome::seqnames(bs_genome)))
  
  i <- 1
  for(chr in BSgenome::seqnames(bs_genome)){
    tri_context[i,] <- Biostrings::trinucleotideFrequency(bs_genome[[chr]])
    i <- i+1
  }
  tri_context_bs_df <- data.frame(
    tri_context_name = rownames(as.data.frame(Biostrings::trinucleotideFrequency(bs_genome[[as.character(BSgenome::seqnames(bs_genome)[1])]]))),
    count=colSums(tri_context)
  ) %>%
    mutate(tri_context_name=as.character(tri_context_name)) %>%
    mutate(count=as.numeric(count)) %>%
    rowwise() %>%
    mutate(tri_context_name=ifelse(substr(tri_context_name,2,2)%in%c("A","G"),findReverseComplement(tri_context_name),tri_context_name)) %>%
    group_by(tri_context_name) %>%
    summarise(count=sum(count)) %>%
    mutate(freq=prop.table(count))
}


getMouseHumanGenomeContext = function(){
  print("Procesing BSgenome.Hsapiens.UCSC.hg19")
  human_genome_context_df = getGenomeContext(BSgenome.Hsapiens.UCSC.hg19)
  print("Procesing BSgenome.Mmusculus.UCSC.mm10")
  mouse_genome_context_df = getGenomeContext(BSgenome.Mmusculus.UCSC.mm10)
  print("Joining Tables")
  tri_context_df = full_join(mouse_genome_context_df,human_genome_context_df,by="tri_context_name",
                             suffix=c(".mm10",".hg19"))
}


plotTriNucleotideContext = function(sigs.input,title_text="Mutation Tri-Nucleotide Context"){
  color_bar = c("#28bdee","#000000","#e62a27","#cbcacb","#a2c966","#ecc8c5")
  
  sigs.input %>%
    ggplot(aes(x=variable,y=value,fill=nt_change)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=color_bar) +
    theme_cowplot() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=2),
      axis.title.x=element_blank()) +
    facet_wrap(~nt_change,nrow=1,scales="free_x") +
    ylab("Count") +
    ggtitle(title_text)
}


