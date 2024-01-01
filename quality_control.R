library(tidyverse)
library(Biostrings)
library(AhoCorasickTrie)
library(cowplot)

work_path <- '/home/projects2/kvs_chunxu_work/gpmdb/'
data_path <- '~/gpmdb/quanti/pro_iso_quantify/'

# Initiate the mapping coverage data_frame which consists of three columns:
# gpm_sample, human_mapping, mouse_mapping
mapping_coverage <- data.frame(
  gpm_sample = readLines(paste0(data_path,'human_samples_batch5.txt'),n = 1000),
  human_mapping = NA,
  mouse_mapping = NA
)

# Mapping peptides to a reference protein sequence database
for(i in 1:nrow(mapping_coverage)){
  sample_name <- mapping_coverage$gpm_sample[i]
  mapping_coverage$human_mapping[i] <- mappingPercent(sample_name,'human')
  mapping_coverage$mouse_mapping[i] <- mappingPercent(sample_name,'mouse')
}

mappingPercent <- function(sample_name,species){
  mapPeptides <- function(peptides,reference) {
    
    ### Test input
    peptides <- as.character(peptides)
    stopifnot(is.character(reference))
    
    ### Helper function
    mfUnlistAndExtract <- function(aList, pattern) {
      tmp <- unlist(aList)
      tmp <- tmp[str_detect(names(tmp), pattern = pattern)]
      return(tmp)
    }
    
    ### Map and format
    AhoCorasickTrie::AhoCorasickSearch(
      keywords = peptides,
      text = reference,
      alphabet = "aminoacid"
    ) %>% 
      enframe(
        name = 'protein',
        value = 'list'
      ) %>% 
      rowwise() %>% 
      reframe(
        protein = protein[1],
        peptide  = mfUnlistAndExtract(list, pattern ='Keyword')
        #position = mfUnlistAndExtract(list, pattern ='Offset'),
      ) %>%  
      return()
  }
  
  ## human or mouse
  if(species == 'human'){
    load(file = paste0(data_path,'01_human_gencode.v44.protein.Rdata')) # gencodeAa, gencode44annotation, gen44pc
    gencodeAaChr <- as.character(gencodeAa)
  } else if(species == 'mouse'){
    load(file = paste0(data_path,'02_mouse_gencode.v33.protein.Rdata')) # gencodeAa, gencode44annotation, gen44pc
    gencodeAaChr <- as.character(gencodeAa)
  }
  
  gpm_peptides <- read_delim(paste0(work_path,sample_name,paste0('/',sample_name,'_peptides.xml')),progress = FALSE)
  mapPepTable <- mapPeptides(peptides = gpm_peptides$sequence,
                             reference = gencodeAaChr)
  mapCoveragePercent <- length(unique(mapPepTable$peptide)) / length(unique(gpm_peptides$sequence))
  mapCoveragePercent <- sprintf("%.2f%%", mapCoveragePercent * 100)
  print(mapCoveragePercent)
  return(mapCoveragePercent)
}

# visualization
# histogram & dot figure
mapping_coverage$human_mapping <- as.numeric(gsub("%", "", mapping_coverage$human_mapping)) / 100
mapping_coverage$mouse_mapping <- as.numeric(gsub("%", "", mapping_coverage$mouse_mapping)) / 100
p1 <- ggplot(mapping_coverage, aes(x=human_mapping)) +
  geom_histogram(binwidth = 0.005, fill="#AEE0F2",color = '#BAB2B4') +
  #geom_vline(xintercept = 0.95, color = "red") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal() +
  labs(title="Human Mapping Histogram", x="Human Mapping Percentage", y="Count")
p2 <- ggplot(mapping_coverage, aes(x=mouse_mapping)) +
  geom_histogram(binwidth = 0.005, fill="#DD9EB8", color="#BAB2B4") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal() +
  labs(title="Mouse Mapping Histogram", x="Mouse Mapping Percentage", y="Count")
mapping_coverage$isHuman <- ifelse(mapping_coverage$human_mapping >= 0.95, "Y", "N")
mapping_coverage$isHuman <- as.factor(mapping_coverage$isHuman)
mapping_coverage <- mapping_coverage %>% filter(!is.na(isHuman))
p3 <- ggplot(mapping_coverage, aes(x=human_mapping, y=mouse_mapping, color = isHuman)) +
  geom_point() +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal() +
  labs(title="Human vs. Mouse Mapping", x="Human Mapping Percentage", y="Mouse Mapping Percentage", color = 'Mapping Coverage (Human) > 95%')


