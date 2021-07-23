
library(dplyr)
library(RColorBrewer)
library(readxl)
library(enrichR)
library(enrichplot)
library(ggplot2)

########################################################################## function ################
####################################################################################################
get_pathway <- function(deg_list, dbs, p_value, idx){
  
  enriched <- enrichr(deg_list, dbs)
  p_value <- 0.05
  
  KEGG_path <- list()
  GO_MF_path <- list()
  GO_CC_path <- list()
  GO_BP_path <- list()
  Biocarta_path <- list()
  Reactome_path <- list()
  
  if (idx == 1) {
    KEGG_path <- enriched$KEGG_2019_Human[which(enriched$KEGG_2019_Human$Adjusted.P.value < p_value),]
    GO_MF_path <- enriched$GO_Molecular_Function_2018[which(enriched$GO_Molecular_Function_2018$Adjusted.P.value < p_value),]
    GO_CC_path <- enriched$GO_Cellular_Component_2018[which(enriched$GO_Cellular_Component_2018$Adjusted.P.value < p_value),]
    GO_BP_path <- enriched$GO_Biological_Process_2018[which(enriched$GO_Biological_Process_2018$Adjusted.P.value < p_value),]
    Biocarta_path <- enriched$BioCarta_2016[which(enriched$BioCarta_2016$Adjusted.P.value < p_value),]
    Reactome_path <- enriched$Reactome_2016[which(enriched$Reactome_2016$Adjusted.P.value < p_value),]

  }
  
  if (idx == 2) {
    KEGG_path <- enriched$KEGG_2019_Human[which(enriched$KEGG_2019_Human$P.value < p_value),]
    GO_MF_path <- enriched$GO_Molecular_Function_2018[which(enriched$GO_Molecular_Function_2018$P.value < p_value),]
    GO_CC_path <- enriched$GO_Cellular_Component_2018[which(enriched$GO_Cellular_Component_2018$P.value < p_value),]
    GO_BP_path <- enriched$GO_Biological_Process_2018[which(enriched$GO_Biological_Process_2018$P.value < p_value),]
    Biocarta_path <- enriched$BioCarta_2016[which(enriched$BioCarta_2016$P.value < p_value),]
    Reactome_path <- enriched$Reactome_2016[which(enriched$Reactome_2016$P.value < p_value),]
  }

  
  if (dim(KEGG_path)[1]!=0) {
    KEGG_path$DB <- "KEGG"
  }
  
  if (dim(GO_MF_path)[1]!=0) {
    GO_MF_path$DB <- "GO_MF"
  }
  
  if (dim(GO_CC_path)[1]!=0) {
    GO_CC_path$DB <- "GO_CC"
  }
  
  if (dim(GO_BP_path)[1]!=0) {
    GO_BP_path$DB <- "GO_BP"
  }
  
  if (dim(Biocarta_path)[1]!=0) {
    Biocarta_path$DB <- "Biocarta"
  }
  
  if (dim(Reactome_path)[1]!=0) {
    Reactome_path$DB <- "Reactome"
  }
  

  enriched_path <- rbind(KEGG_path,  
                         GO_MF_path,              
                         GO_CC_path,              
                         GO_BP_path,
                         Biocarta_path,
                         Reactome_path)
                         
  enriched_path <- enriched_path[enriched_path$P.value < p_value,]
  return(enriched_path)
}



plot_bar <- function(up, down, n_list, title){
  
  
  if(dim(up)[1]!=0 & dim(down)[1]!=0){
    up$type <- "Factor 1"
    down$type <- "Factor 2"
    
    inter_term <- intersect(up$Term, down$Term)
      if (length(inter_term)!=0) { 
        up <- up[-match(inter_term, up$Term),]
        down <- down[-match(inter_term, down$Term),]
      }
    
    
    up <- up[order(up$Combined.Score), ]  # sort
    if(dim(up)[1] <= n_list){ up <- up }
    else{ up <- up [c(1:n_list),] }
    
    
    down <- down[order(down$Combined.Score), ]  # sort
    if (dim(down)[1] <= n_list) { down <- down}
    else{ down <- down [c(1:n_list),] }
    
    
    gos <- rbind(down, up)
    gos$Term <- factor(gos$Term, levels=gos$Term)
     
    ## Diverging Barcharts  
    ggplot(gos, aes(x=Term, y=Combined.Score , label=Combined.Score)) + 
      geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
      scale_fill_manual(
        values = c("#f8766d", "#00ba38"),
        name="Expression", 
        # values = c('up'='#f8766d', 'down'='#00ba38'),
        breaks = c("Factor 1", "Factor 2"),
        labels = c("Factor 1", "Factor 2")
                        ) + 
      
      labs(subtitle="Combined scores from EnrichR", 
           title= title) + 
      coord_flip()
    }else if(dim(down)[1]!=0 & dim(up)[1]==0){
        down$type <- "Factor 2"
        down <- down[order(down$Combined.Score), ]
        if (dim(down)[1] <= n_list) {
          down <- down
          }
        else{
          down <- down [c(1:n_list),]
        }
        gos <- down
        gos$Term <- factor(gos$Term, levels=gos$Term)
        
        # Diverging Barcharts
        ggplot(gos, aes(x=Term, y=Combined.Score , label=Combined.Score)) + 
          geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
          scale_fill_manual(name="Expression", 
                            values = c("Factor 2"="#00ba38")
                            ,labels = c("Factor 2")
                            ) + 
          labs(subtitle="Combined scores from EnrichR", 
               title= title) + 
          coord_flip()
  } else if(dim(up)[1]!=0 & dim(down)[1]==0){
    up$type <- "Factor 1"
    
    up <- up[order(up$Combined.Score), ]  # sort
    if (dim(up)[1] <= n_list) {
      up <- up
      }
    else{
      up <- up [c(1:n_list),]
    }
    
    gos <- up
    gos$Term <- factor(gos$Term, levels=gos$Term)
    
    # Diverging Barcharts
    ggplot(gos, aes(x=Term, y=Combined.Score , label=Combined.Score)) + 
      geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
      scale_fill_manual(name="Expression", 
                        values = c("Factor 1"="#f8766d")
                        ,labels = c("Factor 1")
      ) + 
      labs(subtitle="Combined scores from EnrichR", 
           title= title) + 
      coord_flip()
    
  } 
  
  
  else{
    print('no list')}

}



plot_bar_q <- function(up, down, n_list, title){
  
  
  if(dim(up)[1]!=0 & dim(down)[1]!=0){
    up$type <- "Factor 1"
    down$type <- "Factor 2"
    
    inter_term <- intersect(up$Term, down$Term)
    if (length(inter_term)!=0) { 
      up <- up[-match(inter_term, up$Term),]
      down <- down[-match(inter_term, down$Term),]
    }
    
    
    up <- up[order(-log10(up$Adjusted.P.value)), ]  # sort
    if(dim(up)[1] <= n_list){ up <- up }
    else{ up <- up [c(1:n_list),] }
    
    
    down <- down[order(-log10(down$Adjusted.P.value)), ]  # sort
    if (dim(down)[1] <= n_list) { down <- down}
    else{ down <- down [c(1:n_list),] }
    
    
    gos <- rbind(down, up)
    gos$Term <- factor(gos$Term, levels=gos$Term)
    
    ## Diverging Barcharts  
    ggplot(gos, aes(x=Term, y=-log10(Adjusted.P.value) , label=-log10(Adjusted.P.value))) + 
      geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
      scale_fill_manual(
        values = c("#f8766d", "#00ba38"),
        name="Expression", 
        # values = c('up'='#f8766d', 'down'='#00ba38'),
        breaks = c("Factor 1", "Factor 2"),
        labels = c("Factor 1", "Factor 2")
      ) + 
      
      labs(subtitle="Combined scores from EnrichR", 
           title= title) + 
      coord_flip()
  }else if(dim(down)[1]!=0 & dim(up)[1]==0){
    down$type <- "Factor 2"
    down <- down[order(-log10(down$Adjusted.P.value)), ]
    if (dim(down)[1] <= n_list) {
      down <- down
    }
    else{
      down <- down [c(1:n_list),]
    }
    gos <- down
    gos$Term <- factor(gos$Term, levels=gos$Term)
    
    # Diverging Barcharts
    ggplot(gos, aes(x=Term, y=-log10(Adjusted.P.value) , label=-log10(Adjusted.P.value))) + 
      geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
      scale_fill_manual(name="Expression", 
                        values = c("Factor 2"="#00ba38")
                        ,labels = c("Factor 2")
      ) + 
      labs(subtitle="Combined scores from EnrichR", 
           title= title) + 
      coord_flip()
  } else if(dim(up)[1]!=0 & dim(down)[1]==0){
    up$type <- "Factor 1"
    
    up <- up[order(-log10(up$Adjusted.P.value)), ]  # sort
    if (dim(up)[1] <= n_list) {
      up <- up
    }
    else{
      up <- up [c(1:n_list),]
    }
    
    gos <- up
    gos$Term <- factor(gos$Term, levels=gos$Term)
    
    # Diverging Barcharts
    ggplot(gos, aes(x=Term, y=-log10(Adjusted.P.value) , label=-log10(Adjusted.P.value))) + 
      geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
      scale_fill_manual(name="Expression", 
                        values = c("Factor 1"="#f8766d")
                        ,labels = c("Factor 1")
      ) + 
      labs(subtitle="Combined scores from EnrichR", 
           title= title) + 
      coord_flip()
    
  } 
  
  
  else{
    print('no list')}
  
}



