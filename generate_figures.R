### SUMMARY STATS

library(readr)
library(dplyr)
library(purrr)

## Figure 3: BWA-MEM vs Bowtie 2 default alignments:

bw_filename42 <- paste("~/COLLEGE/winter2018/bio465/project/bwamem/miseq/depth/19_1.5_42.depth.10000window",
                       sep = '')
bw_filename <- paste("~/COLLEGE/winter2018/bio465/project/bwamem/miseq/depth/19_1.5.depth.10000window",
                     sep = '')

bow_filename42 <- paste("~/COLLEGE/winter2018/bio465/project/bowtie2/miseq/depth/20_0_15_2_42.depth.10000window", 
                        sep = '')
bow_filename <- paste("~/COLLEGE/winter2018/bio465/project/bowtie2/miseq/depth/20_0_15_2.depth.10000window",
                      sep = '')

bwa_42 <- read_tsv(bw_filename42)
bwa <- read_tsv(bw_filename)

bow_42 <- read_tsv(bow_filename42)
bow <- read_tsv(bow_filename)

#z <- merge(bwa,bow, by="pos")
z <- merge(bwa_42,bow_42, by="pos")

ggplot(z,aes(pos)) +
  geom_line(aes(y = depth.y, colour = "Bowtie 2"), alpha = 0.95) +
  geom_line(aes(y = depth.x, colour = "BWA-MEM"), alpha = 0.95) +
  labs(x = "Position per 10,000 bp", y = "High Qual Reads") +
  theme(text = element_text(size=18), legend.position = "bottom") +
  scale_colour_manual(name = NULL, values = c("#00B6EB", "#00BA38"))

cor(z$depth.x, z$depth.y)


#### Figure 4: PICARD GRAPHS
setwd("~/COLLEGE/winter2018/bio465/project/bowtie2/miseq/picard_stats/")

data <- 
  do.call("rbind", lapply(list.files(pattern="*_0_15_2.picard_stats"), function(fn) 
    data.frame(Filename=fn, read_tsv(fn, comment = '#', skip = 1) )
  ))

x <- subset(data, CATEGORY!="PAIR")
z <- extract(x, Filename, c("seedSize"), "(\\d+)*", remove=T, convert = T)


ggplot(z, aes(x = seedSize, y = PF_READS_ALIGNED, colour=CATEGORY)) +
  geom_line(size = 2) +
  labs(x = "Seed Size", y = "Reads Aligned") +
  theme(text = element_text(size=18), axis.text.x = element_text(hjust = 1, size = 18)) +
  guides(colour = F)

ggplot(z, aes(x = seedSize, y = PF_MISMATCH_RATE, colour=CATEGORY)) +
  geom_line(size = 2) +
  labs(x = "Seed Size", y = "Mismatch Rate") +
  theme(text = element_text(size=18), axis.text.x = element_text(hjust = 1, size = 18)) +
  guides(colour = F)

ggplot(z, aes(x = seedSize, y = PF_INDEL_RATE, colour=CATEGORY)) +
  geom_line(size = 2) +
  labs(x = "Seed Size", y = "Indel Rate") +
  theme(text = element_text(size=18), axis.text.x = element_text(hjust = 1, size = 18)) +
  guides(colour = F)

ggplot(z, aes(x = seedSize, y = PF_READS_IMPROPER_PAIRS, colour=CATEGORY)) +
  geom_line(size = 2) +
  labs(x = "Seed Size", y = "Improper Pairs") +
  theme(text = element_text(size=18), axis.text.x = element_text(hjust = 1, size = 18)) +
  guides(colour = F)

  # GENERATE LABELS
  #scale_colour_discrete(name = NULL, labels = c("First in Pair", "Second in Pair"))


#### Figure 5: MAPQ

setwd("~/COLLEGE/winter2018/bio465/project/bowtie2/miseq/mapq/")

data_ign <- 
  do.call("rbind", 
          lapply(list.files(pattern="*_0_15_2_ignore_quals.mapquniform"),
                 function(fn) 
                  data.frame(Filename=fn, data.frame(t(read_tsv(fn, col_types = "_n"))))
                 )
  )
z_ign <- extract(data_ign, Filename, c("seedSize"), "(\\d+)*", remove=T, convert = T)


data <- 
  do.call("rbind", 
          lapply(list.files(pattern="*_0_15_2.mapquniform"),
                 function(fn) 
                   data.frame(Filename=fn, data.frame(t(read_tsv(fn, col_types = "_n"))))
          )
  )
z <- extract(data, Filename, c("seedSize"), "(\\d+)*", remove=T, convert = T)

m <- merge(z, z_ign, by = "seedSize")

## FIGURE 5: PERCENT OF MM READS
ggplot(m, aes(x = seedSize)) +
  geom_line(aes(y = (X1.x + X2.x)/X44.x, colour = 'ignore'), size = 2, show.legend = TRUE) +
  geom_line(aes(y = (X1.y + X2.y)/X44.y, colour = 'no ignore'), size = 2, show.legend = TRUE) +
  labs(title = "Percent Reads that Multi-Mapped (0 or 1)",x = "Seed Size", y = "Percent Reads") +
  theme(text = element_text(size=18), axis.text.x = element_text(hjust = 1, size = 18)) +
  scale_colour_manual(name = NULL, values = c("#F8766D", "#619CFF"), breaks=c("no ignore", 'ignore'), labels = c( "Ignore Read Quality", "Use Read Quality")) +
  guides(colour = F)

## PERCENT OF Q42 READS
ggplot(m, aes(x = seedSize)) +
  geom_line(aes(y = (X43.x)/X44.x, colour = 'ignore'), size = 2, show.legend = TRUE) +
  geom_line(aes(y = (X43.y)/X44.y, colour = 'no ignore'), size = 2, show.legend = TRUE) +
  labs(title = "Percent Reads with Perfect Quality (42)",x = "Seed Size", y = "Percent Reads") +
  theme(text = element_text(size=18), axis.text.x = element_text(hjust = 1, size = 18)) +
  scale_colour_manual(name = NULL, values = c("#F8766D", "#619CFF"), breaks=c("no ignore", 'ignore'), labels = c( "Ignore Read Quality", "Use Read Quality")) +
  guides(colour = F)


#   colour x y PANEL group
# 1 #F8766D 1 1     1     1
# 2 #B79F00 2 2     1     2
# 3 #00BA38 3 3     1     3
# 4 #00BFC4 4 4     1     4
# 5 #619CFF 5 5     1     5
# 6 #F564E3 6 6     1     6

