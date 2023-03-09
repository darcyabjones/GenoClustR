df <- read.csv("data-raw/RNAseq_from_Sucher_v2.csv")

rownames(df) <- df$gene
df$gene <- NULL
df

sucher <- as.matrix(df)

usethis::use_data(sucher)
