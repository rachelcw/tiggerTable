library(dplyr)
library(alakazam)
library(tigger)
library(stringr)
library(Biostrings)

germline <- readIgFasta("C:/Users/wrach/OneDrive - Bar-Ilan University/Documents/׳‘׳™׳•׳׳™׳ ׳₪׳•׳¨׳׳˜׳™׳§׳”/׳₪׳¨׳•׳™׳™׳§׳˜/projectPerson1/IGHV_gap_full.fasta")
# get germline reference
dataP1 <- read.table(file = "C:/Users/wrach/OneDrive - Bar-Ilan University/Documents/׳‘׳™׳•׳׳™׳ ׳₪׳•׳¨׳׳˜׳™׳§׳”/׳₪׳¨׳•׳™׳™׳§׳˜/projectPerson1/P1_I64_S1_collapsed.tsv", sep = '\t', header = TRUE)

#get rid of multiple assignment, starts without . - starts in first position, # no N at all, # consensus_count = 2
filtered_data  <- dataP1[!grepl(",", dataP1$v_call),]
filtered_data %>% .[grep("-",filtered_data$sequence_alignment,invert = T),]%>%.[grep("^[^\\.]",filtered_data$sequence_alignment),] %>% .[!grepl("N",filtered_data$sequence_alignment),] %>%.[filtered_data$consensus_count>=2,]

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
#filter data - add v_alleles and v_gene column
data <- filtered_data %>% mutate(v_alleles = getAllele(v_call, strip_d = F, first = T)) %>% 
  mutate(v_gene = getGene(v_call, strip_d = F)) 
genes <- unique(data$v_gene)
final_df <- c()
novel <- c()
novel_list <- list()
for (g in genes) {
  alleles <- unique(data$v_alleles[data$v_gene == g])
  allele_count <- c()
  # Allele count for final table
  for (a in alleles) {
    count <- sum(filtered_data$v_call == a)
    allele_count <- append(allele_count, count)
  }
  #all_novel <- data.frame(matrix(ncol = 0, nrow = 0)) # empty df
  novel_list[[g]] <- c()
  for (a in alleles) {
    sub <- data[data$v_gene == g, ] #data of g
    v_call_sub <- sub$v_call 
    sub$v_call <- a # dirisa test
    
    # run tigger on changed data and return findNovelAlleles answer in novel df
    novel <- findNovelAlleles(sub, germline[a] ,v_call="v_call", j_call="j_call",
                              junction="junction", junction_length="junction_length",
                              germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318)
    #all_novel <- dplyr::bind_rows(all_novel[[g]], novel)
    novel_list[[g]] <- dplyr::bind_rows(novel_list[[g]], novel)
    # input data to final_d
    ######################
    #swap novel
    found_novel <- c()
    for (n in novel$novel_imgt){ # seq in novels tigger
      for (a_new in unique(v_call_sub)[!unique(v_call_sub)%in%a]) {
        # n ==
        allele_name<- names(germline)[germline==n]
        if (length(allele_name)!=0){
          print("1")
          found_novel <- append(found_novel,allele_name)
        }else{
          # n shoter then germline seq
          if(grepl(n, germline[a_new])) {
            print("2")
            allele_name<- names(germline)[grep(n, germline[a_new])]
            found_novel <- c(found_novel,allele_name)
          }else{
            # n longer then germline seq
            if (grepl(germline[a_new], n)) {
              print("3")
              allele_name<- names(germline)[grep(germline[a_new], n)]
              found_novel <- c(found_novel,allele_name)
            }
            else {
              print("no if was good")
            }
          }
        }
      }
    }
    ######################## 
    # create conclusion final table
    new_row <- data.frame(g, paste(alleles, collapse = ","), paste(allele_count, collapse = ","), a, paste(novel$polymorphism_call, collapse = ","), paste(found_novel, collapse = ","), NA)
    names(new_row) <- c("Gene", "Allele", "Count", "Refrence", "Polymorphism_call", "Novel", "notes") #df column titles
    final_df <- dplyr::bind_rows(final_df, new_row) # add row
  }
}


library(data.table)
novel_list_df <- data.table::rbindlist(novel_list)

gm1 <- germline["IGHV4-34*12"]
gm2 <- novel_list_df$novel_imgt[novel_list_df$polymorphism_call=="IGHV4-34*12_T170A"]
gm2 <- novel_list_df$novel_imgt[novel_list_df$polymorphism_call=="IGHV4-34*12_T170A"][1] #first column
gm3 <- germline["IGHV4-34*01"] # refrence

germs <- lapply(c(gm1,gm2,gm3), function(x) strsplit(x, '')[[1]])
germs_m <- t(sapply(germs, `length<-`, max(lengths(germs))))
