library(rvest)
library(tidyverse)

##Deletions

#get the data
dels <- read_html("http://www.umd.be/DMD/4DACTION/W_DMDT1/3") %>% 
  html_nodes('tr') %>% 
  html_text() %>%
  gsub("\\r", ";", .) %>% 
  gsub("\\s+;", ";", .) %>%
  gsub(";\\s+", ";", .) %>%
  gsub("\\#", "n", .)  %>%
  gsub(";$", "", .) %>%
  trimws() %>%
  .[grepl("^p|^Protein",.)] %>% 
  read_delim(., delim = ";", trim_ws = T) %>% 
  rename(Records = `n records`)

dels$Records <- as.integer(dels$Records)

dels$exons <- gsub("Large rearrangementDeletion from exon ", "", dels$Rearrangement)

dels <- dels %>% separate(exons, c('start','end'), sep = ' to ')
dels <- dels %>% mutate(type = 'Deletion')

n_patients_dels <- sum(dels$Records)
allexons_dels <- mapply(getexons, dels$start, dels$end, dels$Records)
exoncounts_dels <- tibble(exon = 1:79, count = as.integer(table(unlist(allexons_dels))), pct_patients = as.integer(table(unlist(allexons_dels))) / n_patients_dels * 100, type = 'Deletion')

#add frameshift info
fr_del <- dels %>% filter(`Mutation type` == 'Fr.')
allexons_dels_fr <- mapply(getexons, fr_del$start, fr_del$end, fr_del$Records)
fr.table <- table(unlist(allexons_dels_fr))
fr_del_counts <- tibble(exon = names(fr.table), count = fr.table, frame = 'Frameshift', type = 'Deletion', type_frame = 'Frameshift Deletion')

infr_del <- dels %>% filter(`Mutation type` == 'InF')
allexons_dels_infr <- mapply(getexons, infr_del$start, infr_del$end, infr_del$Records)
infr.table <- table(unlist(allexons_dels_infr))
infr_del_counts <- tibble(exon = names(infr.table), count = infr.table, frame = 'In-frame', type = 'Deletion', type_frame = 'In-Frame Deletion')


##Duplications
dups <- read_html("http://www.umd.be/DMD/4DACTION/W_DMDT1/4") %>% 
  html_nodes('tr') %>% 
  html_text() %>%
  gsub("\\r", ";", .) %>% 
  gsub("\\s+;", ";", .) %>%
  gsub(";\\s+", ";", .) %>%
  gsub("\\#", "n", .)  %>%
  gsub(";$", "", .) %>%
  trimws() %>%
  .[grepl("^p|^Protein",.)] %>% 
  read_delim(., delim = ";", trim_ws = T) %>% 
  rename(Records = `n records`)

dups$Records <- as.integer(dups$Records)
dups$exons <- gsub("Large .* exon ", "", dups$Rearrangement)

dups <- dups %>% separate(exons, c('start','end'), sep = ' to ') 
dups <- dups %>% filter(str_detect(Rearrangement,'Duplication'))
dups <- dups %>% mutate(end = if_else(str_detect(Rearrangement,'5\''), 1L, as.integer(end)))

dups <- dups %>% mutate(type = 'Duplication')

n_patients_dups <- sum(dups$Records)
allexons_dups <- mapply(getexons, dups$start, dups$end, dups$Records)
exoncounts_dups <- tibble(exon = 1:79, count = as.integer(table(unlist(allexons_dups))), pct_patients = as.integer(table(unlist(allexons_dups))) / n_patients_dups * 100, type = 'Duplication')

#add frameshift info
fr_dup <- dups %>% filter(`Mutation type` == 'Fr.')
allexons_dups_fr <- mapply(getexons, fr_dup$start, fr_dup$end, fr_dup$Records)
fr_dup.table <- table(unlist(allexons_dups_fr))
fr_dup_counts <- tibble(exon = names(fr_dup.table), count = fr_dup.table, frame = 'Frameshift', type = 'Duplication', type_frame = 'Frameshift Duplication')

infr_dup <- dups %>% filter(`Mutation type` == 'InF')
allexons_dups_infr <- mapply(getexons, infr_dup$start, infr_dup$end, infr_dup$Records)
infr_dup.table <- table(unlist(allexons_dups_infr))
infr_dup_counts <- tibble(exon = names(infr_dup.table), count = infr_dup.table, frame = 'In-frame', type = 'Duplication', type_frame = 'In-Frame Duplication')



## merge
deldup <- bind_rows(exoncounts_dels, exoncounts_dups)

n_total <- n_patients_dels + n_patients_dups
deldup <- deldup %>% mutate(pct_from_total = count / n_total * 100)


deldup_frame <- bind_rows(fr_del_counts, infr_del_counts, fr_dup_counts, infr_dup_counts)
deldup_frame$exon <- as.integer(deldup_frame$exon)
n_total_fr <- sum(fr_del$Records) + sum(fr_dup$Records) + sum(infr_del$Records) + sum(infr_dup$Records)
deldup_frame <- deldup_frame %>% mutate(pct_from_total = count / n_total * 100)

deldup_frame$type_frame <- factor(deldup_frame$type_frame, 
                                  levels = c("Frameshift Deletion", "In-Frame Deletion",
                                             "Frameshift Duplication", "In-Frame Duplication"))


#plot
ggplot(deldup, aes(exon, pct_from_total, col = type, fill = type)) +
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  ylab("Pct of patients with\ndeleted/duplicated exon") + 
  theme(legend.position = c(0.2,0.83), legend.title = element_blank())

ggsave('UMD_DMD_dels_dups_stacked.png', width = unit(4, 'in'), height  = unit(3, 'in'))


ggplot(deldup_frame, aes(exon, pct_from_total, col = type_frame, fill = type_frame)) +
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  ylab("Pct of patients with\ndeleted/duplicated exon") + 
  theme(legend.position = c(0.2,0.8), 
        legend.text = element_text(size = 4),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank()) +
  scale_fill_manual(values=c("brown2", "coral", "navy", 'cyan2')) +
  scale_color_manual(values=c("brown2",  "coral", "navy", 'cyan2'))

ggsave('UMD_DMD_dels_dups_frames_stacked.png', width = unit(4, 'in'), height  = unit(3, 'in'))

  
