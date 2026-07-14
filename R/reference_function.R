################ REFERENCE STRUCTURE ANALYSIS


##HEAVY CHAIN WRANGLING

library(stringr)
library(dplyr)
#remove top and bottom extra

reference_analysis<- function(h_file, k_file) {
  h_file<- readLines(h_file)
  h_file_info<- h_file[5:33]
  h_file<- h_file[37:length(h_file)]
  h_file<- h_file[-((length(h_file)-7):length(h_file))]
  
  ###convert to dataframe, separate info into columns
  h_df<- do.call(rbind, strsplit(h_file, "  /"))
  h_df<- as.data.frame(h_df)
  
  #replace empty rows with content of previous row for Ab chains
  h_df$V1<- trimws(h_df$V1)
  h_df$V1[h_df$V1==""] <- NA
  for (i in 2:nrow(h_df)) {
    if (is.na(h_df$V1[i])) {
      h_df$V1[i]<- h_df$V1[i-1]
    }
  }
  
  #prepare ab resis dataframe
  h_abresis<- as.data.frame(str_split_fixed(h_df$V1, "/", 5))
  h_abresis$V4<- substr(h_abresis$V4, 1, 8)
  h_abresis$V4<- str_replace(h_abresis$V4, "\\(", "-")
  h_abresis$V5<- substr(h_abresis$V5, 1, 8)
  h_abresis$dist<- substr(h_df$V2, 28,32)
  h_abresis<- subset(h_abresis, select=c(V4,V5, dist))
  colnames(h_abresis)<- c("ab_resi", "ab_resi_atom", "dist")
  h_abresis$dist<- as.numeric(h_abresis$dist)
  
  #prepare ag resis dataframe
  h_agresis<- as.data.frame(str_split_fixed(h_df$V2, "/", 4))
  h_agresis$dist<- substr(h_agresis$V4, 13,16)
  h_agresis$V4<- substr(h_agresis$V4, 1,8)
  h_agresis$V3<- substr(h_agresis$V3, 1,8)
  h_agresis$V3<- str_replace(h_agresis$V3, "\\(", "-")
  h_agresis<- subset(h_agresis, select= c("V3", "V4", "dist"))
  colnames(h_agresis)<- c("ag_resi", "ag_resi_atom", "dist.")
  h_agresis$dist.<- as.numeric(h_agresis$dist.)
  
  #merge ab and ag resi information
  
  h_ab_ag<- cbind(h_abresis, h_agresis)
  h_contacts<- subset(h_ab_ag, select=c("ab_resi", "ag_resi", "dist"))
  
  
  
  ##LIGHT CHAIN WRANGLING
  
  #remove top and bottom extra
  k_file<- readLines(k_file)
  k_file_info<- k_file[5:33]
  k_file<- k_file[37:length(k_file)]
  k_file<- k_file[-((length(k_file)-7):length(k_file))]
  
  ###convert to dataframe, separate info into columns
  k_df<- do.call(rbind, strsplit(k_file, "  /"))
  k_df<- as.data.frame(k_df)
  
  #replace empty rows with content of previous row for Ab chains
  k_df$V1<- trimws(k_df$V1)
  k_df$V1[k_df$V1==""] <- NA
  for (i in 2:nrow(k_df)) {
    if (is.na(k_df$V1[i])) {
      k_df$V1[i]<- k_df$V1[i-1]
    }
  }
  
  #prepare ab resis dataframe
  k_abresis<- as.data.frame(str_split_fixed(k_df$V1, "/", 5))
  k_abresis$V4<- substr(k_abresis$V4, 1, 8)
  k_abresis$V4<- str_replace(k_abresis$V4, "\\(", "-")
  k_abresis$V5<- substr(k_abresis$V5, 1, 8)
  k_abresis<- subset(k_abresis, select=c(V4,V5))
  colnames(k_abresis)<- c("ab_resi", "ab_resi_atom")
  
  
  #prepare ag resis dataframe
  k_agresis<- as.data.frame(str_split_fixed(k_df$V2, "/", 4))
  k_agresis$dist<- substr(k_agresis$V4, 13,16)
  k_agresis$V4<- substr(k_agresis$V4, 1,8)
  k_agresis$V3<- substr(k_agresis$V3, 1,8)
  k_agresis$V3<- str_replace(k_agresis$V3, "\\(", "-")
  k_agresis<- subset(k_agresis, select=c("V3", "V4", "dist"))
  colnames(k_agresis)<- c("ag_resi", "ag_resi_atom", "dist")
  
  
  #merge ab and ag resi information
  k_ab_ag<- cbind(k_abresis, k_agresis)
  k_contacts<- subset(k_ab_ag, select=c("ab_resi", "ag_resi", "dist"))
  k_contacts$dist<- as.numeric(k_contacts$dist)
  
  ###now to get stats on the ab and ag resis....
  
  #heavy chain ab resis
  h_ab_table<- as.data.frame(table(h_contacts$ab_resi))
  colnames(h_ab_table)<- c("ab_resi", "count")
  h_ab_table<- merge(h_ab_table, aggregate(dist~ ab_resi, h_contacts, mean), by="ab_resi")
  h_ab_table$nod <- h_ab_table$count / h_ab_table$dist
  
  #ag resis interacting with heavy chain
  h_ag_table<- as.data.frame(table(h_contacts$ag_resi))
  colnames(h_ag_table)<- c("ag_resi", "count")
  h_ag_table<- merge(h_ag_table, aggregate(dist~ ag_resi, h_contacts, mean), by="ag_resi")
  h_ag_table$nod <- h_ag_table$count / h_ag_table$dist
  
  #light chain ab resis
  k_ab_table<- as.data.frame(table(k_contacts$ab_resi))
  colnames(k_ab_table)<- c("ab_resi", "count")
  k_ab_table<- merge(k_ab_table, aggregate(dist~ ab_resi, k_contacts, mean), by="ab_resi")
  k_ab_table$nod <- k_ab_table$count / k_ab_table$dist
  
  #ag resis interacting with light chain
  k_ag_table<- as.data.frame(table(k_contacts$ag_resi))
  colnames(k_ag_table)<- c("ag_resi", "count")
  k_ag_table<- merge(k_ag_table, aggregate(dist~ ag_resi, k_contacts, mean), by="ag_resi")
  k_ag_table$nod <- k_ag_table$count / k_ag_table$dist
  
  
  ##need to compute the quartiles new for the reference structure
  
  #first start by creating columns of all H and k ab resis or ag interacting residues
  ag_dist_quart<- rbind(h_ag_table, k_ag_table)
  ab_dist_quart<- rbind(h_ab_table, k_ab_table)
  
  ag_quart<- quantile(ag_dist_quart$nod, probs= c(1, 0.75, 0.5, 0.25))
  ab_quart<- quantile(ab_dist_quart$nod, probs= c(1, 0.75, 0.5, 0.25))
  
  
  ###grouping ab and ag residues based on their n/d ratios compared to reference quartiles (to be input at the beginning of the script)
  
  #Ab quartiles 
  #(Alphafold predicted paratope (AFPP) designated as Q0 or Q1 residues)
  for (i in 1: length(h_ab_table$nod)) {
    if (h_ab_table$nod[i] >= ab_quart[1]) {
      h_ab_table$Q[i]<- "Q0"
    }
    else if (h_ab_table$nod[i] >= ab_quart[2]) {
      h_ab_table$Q[i]<- "Q1"
    }
    else if (h_ab_table$nod[i] >= ab_quart[3]) {
      h_ab_table$Q[i]<- "Q2"
    }
    else if (h_ab_table$nod[i] >= ab_quart[4]) {
      h_ab_table$Q[i]<- "Q3"
    }
    else {
      h_ab_table$Q[i]<- "Q4"
    }
  }
  
  for (i in 1: length(k_ab_table$nod)) {
    if (k_ab_table$nod[i] >= ab_quart[1]) {
      k_ab_table$Q[i]<- "Q0"
    }
    else if (k_ab_table$nod[i] >= ab_quart[2]) {
      k_ab_table$Q[i]<- "Q1"
    }
    else if (k_ab_table$nod[i] >= ab_quart[3]) {
      k_ab_table$Q[i]<- "Q2"
    }
    else if (k_ab_table$nod[i] >= ab_quart[4]) {
      k_ab_table$Q[i]<- "Q3"
    }
    else {
      k_ab_table$Q[i]<- "Q4"
    }
  }
  
  
  #Ag quartiles
  #(Alphafold predicted epitope (AFPE) designated as Q0 or Q1 residues)
  for (i in 1: length(h_ag_table$nod)) {
    if (h_ag_table$nod[i] >= ag_quart[1]) {
      h_ag_table$Q[i]<- "Q0"
    }
    else if (h_ag_table$nod[i] >= ag_quart[2]) {
      h_ag_table$Q[i]<- "Q1"
    }
    else if (h_ag_table$nod[i] >= ag_quart[3]) {
      h_ag_table$Q[i]<- "Q2"
    }
    else if (h_ag_table$nod[i] >= ag_quart[4]) {
      h_ag_table$Q[i]<- "Q3"
    }
    else {
      h_ag_table$Q[i]<- "Q4"
    }
  }
  
  for (i in 1: length(k_ag_table$nod)) {
    if (k_ag_table$nod[i] >= ag_quart[1]) {
      k_ag_table$Q[i]<- "Q0"
    }
    else if (k_ag_table$nod[i] >= ag_quart[2]) {
      k_ag_table$Q[i]<- "Q1"
    }
    else if (k_ag_table$nod[i] >= ag_quart[3]) {
      k_ag_table$Q[i]<- "Q2"
    }
    else if (k_ag_table$nod[i] >= ag_quart[4]) {
      k_ag_table$Q[i]<- "Q3"
    }
    else {
      k_ag_table$Q[i]<- "Q4"
    }
  }
  
  ref_output<- list("Ag quartiles"=ag_quart, "Ab quartiles"=ab_quart, 
                    "Antibody H residues"= h_ab_table, "Antibody k residues"=k_ab_table, 
                    "Antigen residues H interacting"= h_ag_table, "Antigen residues k interacting"=k_ab_table)
  return(ref_output) 

}
