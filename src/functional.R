# library(dplyr)
# library(glue)
library(reutils)
# fname <- "~/git/research/metagenome/gene_coding/lactobacillus_johnsonii/other/synonymous.tab"
# synonymous <- read.table(fname, sep="\t", header = FALSE)
# synonymous_noref <- filter(synonymous, V2 != "REF") %>% group_by(V3) %>% summarise(samples=paste0(V2, collapse = ","))
# nrow(synonymous_noref)

#clostridium_sporogenes  enterobacteria_phage_lambda  enterococcus_faecalis lactobacillus_johnsonii
strain <- "clostridium_sporogenes"
fname <- glue("~/git/research/metagenome/gene_coding/{strain}/treatment/nonsynonymous.tab")
nonsynonymous <- read.table(fname, sep="\t", header = FALSE)
nonsynonymous_noref <- filter(nonsynonymous, V2 != "REF") %>% group_by(V3) %>% summarise(samples=paste0(V2, collapse = ","))
nrow(nonsynonymous_noref) - length(intersect(as.vector(synonymous_noref$V3), as.vector(nonsynonymous_noref$V3)))

positive <- filter(nonsynonymous_noref, !grepl(pattern="N", samples))
s <- strsplit(positive$samples, split = ",")
positive <- data.frame(V1 = rep(positive$V3, sapply(s, length)), sample = unlist(s))

negative <- filter(nonsynonymous_noref, !grepl(pattern="P", samples))
s <- strsplit(negative$samples, split = ",")
negative <- data.frame(V1 = rep(negative$V3, sapply(s, length)), sample = unlist(s))

columns <- c("QueryID", "SubjectID", "PercIdentity",
             "AlignmentLength", "Mismatches", "GapOpens",
             "QueryStart", "QueryEnd", "SubjectStart",
             "SubjectEnd", "Evalue", "BitScore")

data <- negative

pos <- data[1, 1]
sample <- data[1, 2]
df <- read.table(file = glue("~/git/research/metagenome/gene_coding/{strain}/treatment/{sample}/treatment_{pos}_blast.txt"),
                 sep = "\t", comment.char = "#")
colnames(df) <- columns

for (i in 2:nrow(data)) {
  pos <- data[i, 1]
  sample <- data[i, 2]
  tmp <- read.table(file = glue("~/git/research/metagenome/gene_coding/{strain}/treatment/{sample}/treatment_{pos}_blast.txt"),
                    sep = "\t", comment.char = "#")
  colnames(tmp) <- columns
  df <- rbind(df, tmp)
}

write.table(x = df,
            file = glue("~/git/research/metagenome/gene_coding/{strain}/onlynegative_nonsynonymous.tab"), sep="\t", row.names = FALSE)


# library(ggplot2)
# x <- c(1,2,3,4,5,6,7)
# y <- c(1,1,2,2,1,1,1)
# w <- c(10,2,3,1,7,2,1)
# label = c("h1=10", "d1=2", "d2=3","h2=1","d3=7","h3=2","d4=1")
# distribution <- as.factor(c("holes","dirts","dirts","holes","dirts","holes","dirts"))
# 
# df <- data.frame(x=x,y=y,w=w,distribution=distribution)
# 
# ggplot(df, aes(x=x, y=y, color=distribution, label=label))+ geom_text(vjust = 0, nudge_y = 0.1)+ geom_point(size=w+3) + 
#   scale_x_continuous(limits = c(0,8), breaks=seq(0,8,1)) + scale_y_continuous(limits = c(0,2.5), breaks=seq(0,2,0.5))

# clostridium_sporogenes  enterobacteria_phage_lambda  enterococcus_faecalis lactobacillus_johnsonii
# clostridium_sporogenes GCF_000960175.1_ASM96017v1
# enterococcus_faecalis GCF_902161805.1_25426_7_320
# lactobacillus_johnsonii GCF_003316915.1_ASM331691v1

strain = "lactobacillus_johnsonii"
type = "treatment"
group = "positive"
assembly = "GCF_003316915.1"

df = read.csv(glue("~/git/research/metagenome/gene_coding/{strain}_{type}_only{group}_nonsynonymous_annotation.tab"), sep="\t")
df = df[df$Assembly == assembly,]

dir.create(file.path(glue("/tmp/{strain}/{type}/{group}/")), showWarnings = TRUE, recursive = TRUE)

for (i in df$Protein) {
  print(i)
  p = efetch(i, db="protein", rettype = "gp", retmode = "text", outfile = glue("/tmp/{strain}/{type}/{group}/{i}.gb"))
}



library(clusterProfiler)
library(dplyr)
library(glue)
search_kegg_organism("Clostridium sporogenes", by="scientific_name")
search_kegg_organism("Enterococcus faecalis", by="scientific_name")
search_kegg_organism("Lactobacillus johnsonii", by="scientific_name")

# kegg_code                    scientific_name common_name
# 2907       efa         Enterococcus faecalis V583        <NA>
# 2908       efl           Enterococcus faecalis 62        <NA>
# 2909       efi        Enterococcus faecalis OG1RF        <NA>
# 2910       efd          Enterococcus faecalis D32        <NA>
# 2911       efs Enterococcus faecalis Symbioflor 1        <NA>
# 2912       efn        Enterococcus faecalis DENG1        <NA>
# 2913       efq   Enterococcus faecalis ATCC 29212        <NA>

# kegg_code                  scientific_name common_name
# 2828       ljo  Lactobacillus johnsonii NCC 533        <NA>
# 2829       ljf   Lactobacillus johnsonii FI9785        <NA>
# 2830       ljh Lactobacillus johnsonii DPC 6026        <NA>
# 2831       ljn     Lactobacillus johnsonii N6.2        <NA>

############################################ Clostridium sporogenes
strain = "cld"
clostridium_sporogenes = read.csv(glue("http://rest.kegg.jp/list/{strain}"), header = FALSE, sep = "\t")
colnames(clostridium_sporogenes) = c("KEEGId", "Desc")

# other
other_genes_negative = c("galU", "gltA", "trpS")
other_genes_positive = c("efp", "trpS", "rlmN", "ade", "purE", "glmS", "flhA", "gcvT", "dhaK", "gltA",
                         "fliY", "thiF", "buk", "bcmB", "gap", "murJ", "rapZ", "spoIIR", "prdB", "yjjW",
                         "prdB", "nadA", "dapD")
# treatment
treatment_genes_positive = c("trpS", "rlmN", "rsgA", "purE", "glmS", "flhA", "gcvT", "dhaK", "gltA", "fliY",
                             "minD", "pepV", "thiF", "buk", "bcmB", "cobO", "fliW", "murJ", "murJ", "rapZ", 
                             "yjjW", "nadA")

genes = treatment_genes_positive

df = filter(clostridium_sporogenes, grepl(pattern=genes[1], Desc))
for (i in 2:length(genes)) {
  tmp = filter(clostridium_sporogenes, grepl(pattern=genes[i], Desc))
  if (nrow(tmp) == 0) {
    cat(genes[i])
    cat("\n")
  }
  df = rbind(df, tmp)
}
lst = sapply(as.vector(unique(df$KEEGId)), function(x) strsplit(x, ":")[[1]], USE.NAMES = FALSE)[2,]
kk <- enrichKEGG(gene         = lst,
                 organism     = strain,
                 pvalueCutoff = 0.05)
head(kk)

write.csv(kk, file=glue("/tmp/clostridium_sporogenes_{strain}_other_genes_negative_kegg.csv"), row.names = FALSE)
write.csv(kk, file=glue("/tmp/clostridium_sporogenes_{strain}_other_genes_positive_kegg.csv"), row.names = FALSE)


############################################ Enterococcus faecalis
for (strain in c("efa", "efl", "efi", "efd", "efs", "efn", "efq")) {
  enterococcus_faecalis = read.csv(glue("http://rest.kegg.jp/list/{strain}"), header = FALSE, sep = "\t")
  colnames(enterococcus_faecalis) = c("KEEGId", "Desc")
  
  
  other_genes_negative = c("gshAB")
  other_genes_positive = c("dnaE", "tsaD", "moaA", "budA", "mgtE", "lepA", "rfbA", "cydC", "nrdD", "dfrE", "fss1")
  
  treatment_genes_negative = c("gshAB", "recQ", "addA", "nifJ")
  treatment_genes_positive = c("ezrA", "pta", "dnaE", "walK", "mutS", "recX", "purB", "rimM", "mgtE", "cydC", "ruvB", "pheA", "fss1")
  
  for (grp in c("other_genes_negative", "other_genes_positive", "treatment_genes_negative", "treatment_genes_positive")) {
    genes = get(grp)
    df = filter(enterococcus_faecalis, grepl(pattern=genes[1], Desc))
    if (nrow(df) > 0) {
      if (nrow(df) > 1) {
        for (i in 2:length(genes)) {
          tmp = filter(enterococcus_faecalis, grepl(pattern=genes[i], Desc))
          if (nrow(tmp) == 0) {
            cat(genes[i])
            cat("\n")
          }
          df = rbind(df, tmp)
        }
      }
      lst = sapply(as.vector(unique(df$KEEGId)), function(x) strsplit(x, ":")[[1]], USE.NAMES = FALSE)[2,]
      kk <- enrichKEGG(gene         = lst,
                       organism     = strain,
                       pvalueCutoff = 0.05)
      if (nrow(kk) > 0 && !is.null(kk)) {
        write.csv(kk, file=glue("/tmp/enterococcus_faecalis_{strain}_{grp}_kegg.csv"), row.names = FALSE)
        
      }
      
    }
  }
}


#################################################### Lactobacillus johnsonii
for (strain in c("ljo", "ljf", "ljh", "ljn")) {
  lactobacillus_johnsonii = read.csv(glue("http://rest.kegg.jp/list/{strain}"), header = FALSE, sep = "\t")
  colnames(lactobacillus_johnsonii) = c("KEEGId", "Desc")
  
  
  other_genes_negative = c("gap", "glyQ", "nusG", "pepV", "infB", "miaA", "recF", "rpoC", "yjeM", "topA", "xerC", "manX")
  other_genes_positive = c("tuf", "phoU", "pyk", "greA", "recN", "rpoC", "spxA")
  
  treatment_genes_negative = c("infB", "miaA", "recF", "xerC")
  treatment_genes_positive = c("tuf", "pyk", "greA", "recN")
  
  for (grp in c("other_genes_negative", "other_genes_positive", "treatment_genes_negative", "treatment_genes_positive")) {
    genes = get(grp)
    df = filter(lactobacillus_johnsonii, grepl(pattern=genes[1], Desc))
    if (nrow(df) > 0) {
      if (nrow(df) > 1) {
        for (i in 2:length(genes)) {
          tmp = filter(lactobacillus_johnsonii, grepl(pattern=genes[i], Desc))
          if (nrow(tmp) == 0) {
            cat(genes[i])
            cat("\n")
          }
          df = rbind(df, tmp)
        }
      }
      lst = sapply(as.vector(unique(df$KEEGId)), function(x) strsplit(x, ":")[[1]], USE.NAMES = FALSE)[2,]
      kk <- enrichKEGG(gene         = lst,
                       organism     = strain,
                       pvalueCutoff = 0.05)
      if (nrow(kk) > 0 && !is.null(kk)) {
        write.csv(kk, file=glue("/tmp/lactobacillus_johnsonii_{strain}_{grp}_kegg.csv"), row.names = FALSE)
      }
    }
  }
}