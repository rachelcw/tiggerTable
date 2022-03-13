test_gene <- function(allele_name, refrence_allele, data, germline) {
  final_df_gene <- c()
  #novel <- c()
  #for (g in genes) {
  allele_count <- c()
  # Allele count for final table
  #if allele != null
  #alleles <- unique(data$v_alleles[data$v_gene == gene_name])
  #for (a in alleles) {
  #  count <- sum(filtered_data$v_call == a)
  #  allele_count <- c(allele_count, count)
  #}
  #for (a in alleles) {
  gene_name <- substr(allele_name,1,nchar(allele_name)-3)
  
  test_y <- data.frame()
  
  sub <- data[data$v_gene == gene_name, ] #data of g
  v_call_original <- sub$v_call 
  sub$v_call <- refrence_allele # swap all the alleles of g to a
  
  # run tigger on changed data and return findNovelAlleles answer in novel df
  # def 0.125
  novel_0 <- findNovelAlleles(sub, germline[refrence_allele], v_call="v_call", j_call="j_call",
                              junction="junction", junction_length="junction_length",
                              germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318)
  test_y <- rbind(test_y, novel_0)
  # 0.1
  novel_1 <- findNovelAlleles(sub, germline[refrence_allele] ,v_call="v_call", j_call="j_call",
                              junction="junction", junction_length="junction_length",
                              germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318, y_intercept = 0.1)
  test_y <- rbind(test_y, novel_1)
  # 0.09
  novel_09 <- findNovelAlleles(sub, germline[refrence_allele] ,v_call="v_call", j_call="j_call",
                               junction="junction", junction_length="junction_length",
                               germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318, y_intercept = 0.09)
  test_y <- rbind(test_y, novel_09)
  # 0.08
  novel_08 <- findNovelAlleles(sub, germline[refrence_allele] ,v_call="v_call", j_call="j_call",
                               junction="junction", junction_length="junction_length",
                               germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318, y_intercept = 0.08)
  test_y <- rbind(test_y, novel_08)
  # 0.07
  novel_07 <- findNovelAlleles(sub, germline[refrence_allele] ,v_call="v_call", j_call="j_call",
                               junction="junction", junction_length="junction_length",
                               germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318, y_intercept = 0.07)
  test_y <- rbind(test_y, novel_07)
  # 0.06
  novel_06 <- findNovelAlleles(sub, germline[refrence_allele] ,v_call="v_call", j_call="j_call",
                               junction="junction", junction_length="junction_length",
                               germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318, y_intercept = 0.06)
  test_y <- rbind(test_y, novel_06)
  # 0.05
  novel_05 <- findNovelAlleles(sub, germline[refrence_allele] ,v_call="v_call", j_call="j_call",
                               junction="junction", junction_length="junction_length",
                               germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318, y_intercept = 0.05)
  test_y <- rbind(test_y, novel_05)
  # 0.04
  novel_04 <- findNovelAlleles(sub, germline[refrence_allele] ,v_call="v_call", j_call="j_call",
                               junction="junction", junction_length="junction_length",
                               germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318, y_intercept = 0.04)
  test_y <- rbind(test_y, novel_04)
  # 0.03
  novel_03 <- findNovelAlleles(sub, germline[refrence_allele] ,v_call="v_call", j_call="j_call",
                               junction="junction", junction_length="junction_length",
                               germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318, y_intercept = 0.03)
  test_y <- rbind(test_y, novel_03)
  # 0.02
  novel_02 <- findNovelAlleles(sub, germline[refrence_allele] ,v_call="v_call", j_call="j_call",
                               junction="junction", junction_length="junction_length",
                               germline_min=1, seq="sequence_alignment", min_seqs = 1, pos_range = 1:318,  y_intercept = 0.02)
  test_y <- rbind(test_y, novel_02)
  View(test_y)
  }
  
###################################################################################
# # vec of strings racheli
germline <- readIgFasta("C:/Users/wrach/OneDrive - Bar-Ilan University/Documents/ביואינפורמטיקה/פרוייקט/tiggerTable/IGHV_gap_full.fasta")
# # get germline reference
dataP1 <- read.table(file = "C:/Users/wrach/OneDrive - Bar-Ilan University/Documents/ביואינפורמטיקה/פרוייקט/tiggerTable/P1_I64_S1_collapsed.tsv", sep = '\t', header = TRUE)

# vec of strings eden
#germline <- readIgFasta("C:/Users/Eden/OneDrive - Bar-Ilan University/Desktop/tigger/IGHV_gap_full.fasta")
# get germline reference
#dataP1 <- read.table(file = "C:/Users/Eden/OneDrive - Bar-Ilan University/Desktop/tigger/P1_I64_S1_collapsed.tsv", sep = '\t', header = TRUE)
#get rid of multiple assignment, starts without . - starts in first position, # no N at all, # consensus_count = 2
filtered_data  <- dataP1[!grepl(",", dataP1$v_call),]
filtered_data  <-filtered_data %>% filter(!grepl("^[.]", sequence_alignment), consensus_count>=2) %>% rowwise() %>% mutate(v_seq = substr(sequence_alignment, 1, 318)) %>% filter(!grepl("N", v_seq), !grepl("-", v_seq))
#filter data - add v_alleles and v_gene column
# data -> person
data <- filtered_data %>% mutate(v_alleles = getAllele(v_call, strip_d = F, first = T)) %>% mutate(v_gene = getGene(v_call, strip_d = F)) 
###################################################################################

allele_name <- "IGHV4-34*02"
refrence_allele <- "IGHV4-34*01"
# test func
test_gene(allele_name, refrence_allele, data = data, germline = germline)