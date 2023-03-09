df <- read.csv("data-raw/Ssclerotiorum_FPKM.csv")

rownames(df) <- df$gene
df$gene <- NULL

sclero <- as.matrix(df)

usethis::use_data(sclero)
