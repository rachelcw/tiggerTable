#1. get only the gene instead of whole vcall (mutate func) 
# and group for each v_gene all v_calls (group by & summarise)
# transforms it into a vector (pull)
gene_allele <- filtered_data %>% mutate(v_gene = getGene(v_call, strip_d = F)) %>% group_by(v_gene) %>% summarise(v_call = paste0(unique(v_call), collapse = ","))
#new tigger - 
#allele_diff <- function(germs) {
#  germs <- lapply(germs, function(x) strsplit(x, '')[[1]])
#  germs_m <- t(sapply(germs, `length<-`, max(lengths(germs))))
#  setdiff_mat <- function(x) {
#    sum(!unique(x) %in% c('.', NA, "N"))#
#  }
#  idx = which(apply(germs_m, 2, setdiff_mat) > 1)
#  return(idx)
#}


# gm1 <- germline["IGHV4-34*12"]
# gm2 <- novel_list_df$novel_imgt[novel_list_df$polymorphism_call=="IGHV4-34*12_T170A"]
# print(gm2)
# gm2 <- novel_list_df$novel_imgt[novel_list_df$polymorphism_call=="IGHV4-34*12_T170A"][1] #first column
# print(gm2)
# gm3 <- germline["IGHV4-34*01"] # refrence
# 
# germs <- lapply(c(gm1,gm2,gm3), function(x) strsplit(x, '')[[1]])
# germs_m <- t(sapply(germs, `length<-`, max(lengths(germs))))