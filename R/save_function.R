############ SAVE FUNCTION

save_csv<- function(output, interaction_name){

  #first label every table with a H or k column to tell them apart
  h_ab_table<- cbind(chain="H", ouput[[3]])
  k_ab_table<- cbind(chain="k", output[[4]])
  h_ag_table<- cbind(int_chain="H", output[[5]])
  k_ag_table<- cbind(int_chain="k", output[[6]])

  #will make a mega table with all of these next to each other
  ab_output<- rbind(h_ab_table, k_ab_table)
  ag_output<- rbind(h_ag_table, k_ag_table)

  write.csv(ab_output, paste0(output_name, "_antibody_afpp.csv"))
  write.csv(ag_output, paste0(output_name, "_antigen_afpe.csv"))

  print(paste0("Files saved as:\n ", paste0(output_name, "_antibody_afpp.csv/n", output_name, "_antigen_afpe.csv")))
}
