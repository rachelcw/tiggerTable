library(dplyr)
library(alakazam)
library(tigger)
library(stringr)
library(Biostrings)
source('functions_tigger.R')


germline <- readIgFasta("C:/Users/wrach/OneDrive - Bar-Ilan University/Documents/ביואינפורמטיקה/פרוייקט/tiggerTable/IGHV_gap_full.fasta")
# get germline reference
dataP1 <- read.table(file = "C:/Users/wrach/OneDrive - Bar-Ilan University/Documents/ביואינפורמטיקה/פרוייקט/tiggerTable/P1_I64_S1_collapsed.tsv", sep = '\t', header = TRUE)
gene_allele_all <- dataP1 %>% mutate(v_gene = getGene(v_call, strip_d = F)) %>% group_by(v_gene) %>% summarise(v_call = paste0(unique(v_call), collapse = ","))


#get rid of multiple assignment, starts without . - starts in first position, # no N at all, # consensus_count = 2
filtered_data  <- dataP1[!grepl(",", dataP1$v_call),]
filtered_data  <-filtered_data %>% filter(!grepl("^[.]", sequence_alignment), consensus_count>=2) %>% rowwise() %>% mutate(v_seq = substr(sequence_alignment, 1, 318)) %>% filter(!grepl("N", v_seq), !grepl("-", v_seq))

gene_allele_filter <- filtered_data %>% mutate(v_gene = getGene(v_call, strip_d = F)) %>% group_by(v_gene) %>% summarise(v_call = paste0(unique(v_call), collapse = ","))
#filter data - add v_alleles and v_gene column
# data -> person
data <- filtered_data %>% mutate(v_alleles = getAllele(v_call, strip_d = F, first = T)) %>% mutate(v_gene = getGene(v_call, strip_d = F)) 
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
    allele_count <- c(allele_count, count)
  }
  
  novel_list[[g]] <- c()
  for (a in alleles) {
    sub <- data[data$v_gene == g, ] #data of g
    v_call_original <- sub$v_call 
    sub$v_call <- a # swap all the alleles of g to a
    # if (length(alleles)==1){
    #   #we have only one allele (homozygous)
    #   sub$v_call<-gsub(" ", "", paste(g, '*01'))
    #   
    # }else{
    #   sub$v_call <- a # swap all the alleles of g to a
    # }
    
    # run tigger on changed data and return findNovelAlleles answer in novel df
    novel <- findNovelAlleles(sub, germline[a] ,v_call="v_call", j_call="j_call",
                              junction="junction", junction_length="junction_length",
                              germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318)
    #all_novel <- dplyr::bind_rows(all_novel[[g]], novel)
    novel_list[[g]] <- dplyr::bind_rows(novel_list[[g]], novel)
    
    ######################
    #check why tigger didnt find- length#
    found_novel <- c()
    note<- c()
    for (n in novel$novel_imgt){ # seq in novels from tigger function
      for (a_new in unique(v_call_original)[!unique(v_call_original)%in%a]) {
        # n is the seq of the novel 
        allele_name<- names(germline)[germline==n]
        if (length(allele_name)!= 0){
          # n is found in the reference and we know the specific name of the allele
          found_novel <- c(found_novel,allele_name) 
          note<-c(note, "found")
        }else{
          ## n is didnt find in reference and we want to understand why-##
          # n shorter than germline seq
          if(grepl(n, germline[a_new])) {
            allele_name<- names(germline)[grep(n, germline[a_new])]
            found_novel <- c(found_novel,allele_name)
            note<-c(note, "shorter")
          }else{
            # n longer than germline seq
            if (grepl(germline[a_new], n)) {
              allele_name<- names(germline)[grep(germline[a_new], n)]
              found_novel <- c(found_novel,allele_name)
              note<-c(note, "longer")
            }
            else {
              note<-c(note, "Suspected")
              # we need to analyze the reason
            }
         }
       }
      }
    }
    ######################## 
    # create conclusion final table
    new_row <- data.frame(g, paste(alleles, collapse = ","), paste(allele_count, collapse = ","), a, paste(novel$polymorphism_call, collapse = ","), paste(found_novel, collapse = ","), paste(note, collapse = ","))
    names(new_row) <- c("Gene", "Allele", "Count", "Refrence", "Polymorphism_call", "Novel", "notes") #df column titles
    final_df <- dplyr::bind_rows(final_df, new_row) # add row
  }
}


#write.csv(final_df,"C:/Users/wrach/OneDrive - Bar-Ilan University/Documents/ביואינפורמטיקה/פרוייקט/tiggerTable/final_df20_02.csv", row.names = FALSE)

library(data.table)
novel_list_df <- data.table::rbindlist(novel_list) #df with the results from findNovelAllele
#write.csv(novel_list_df,"C:/Users/wrach/OneDrive - Bar-Ilan University/Documents/ביואינפורמטיקה/פרוייקט/tiggerTable/novel_list_df20_02.csv", row.names = FALSE)



# check<-data[data$v_call == "IGHV4-4*02",]
# check<-check$sequence_alignment[1]
# check
# gm1 <- germline["IGHV4-4*02"]
# gm2 <- check
# print(gm2)
# gm2 <- novel_list_df$novel_imgt[novel_list_df$polymorphism_call=="IGHV4-34*12_T170A"][1] #first column
# print(gm2)
# gm3 <- germline["IGHV4-34*01"] # refrence
# 
# germs <- lapply(c(gm1,gm2), function(x) strsplit(x, '')[[1]])
# germs_m <- t(sapply(germs, `length<-`, max(lengths(germs))))
