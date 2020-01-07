########
library(readxl)
library(readr)
library(xlsx)


suppTables = list.files("suppTables", full = TRUE)

harmList = lapply(harmFiles, read.csv,as.is=TRUE, row.names=1)
for(i in seq(along=harmList)) {
	write.xlsx(harmList[[i]], file = "supp_tables/Extended_Data_Figure3-4.xlsx",
		sheetName = names(harmList)[i],append=TRUE)
}