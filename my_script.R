install.packages(c("ggplot2", "tibble"))
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BioManager")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
BiocManager::install(c("rtracklayer","GenomicRanges"))
BiocManager::available("TxDb.Hsapiens")
ls()
rm(list = ls())
ls()
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("IRanges")
library(IRanges)
rng <- IRanges(start=4, end=13)
print(rng)
x <- IRanges(start=c(4,7,2,20), end=c(13,7,5,23))
print(x)
names(x) <- letters[1:4]
end(x) <- end(x) +4
start(x) <5
x[start(x)<5]
x <- IRanges(start = c(40, 80), end = c(67, 114))
x +4L
flank(x, width = 7)
flank(x, width = 7, start = FALSE)
library(ggplot2)
library(IRanges)
set.seed(0)
alns <- IRanges(start=sample(seq_len(50), 20), width = 5)
head (alns, 4)
reduce(alns)
gaps(alns)
qry <- IRanges(start=c(1,26,19,11,21,7), end=c(16,30,19,15,24,8), names=letters[1:6])
sbj <- IRanges(start = c(1,19,10), end = c(5,29,16), naems=letters[24:26])
hts <- findOverlaps(ary, sbj)
hts <- findOverlaps(qry, sbj)
hts
hts_within <- findOverlaps(qry, sbj, type="widthin")
hts_within <- findOverlaps(qry, sbj, type="within")
hts_within
BiocManager::install(“GenomicFeatures”)
library(BiocManager)
library(GenomicRanges)
library(ggplot2)
library(tibble)
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")
library(GenomicFeatures)
BiocManager::available("TxDb.Hsapiens")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
mm_genes <- genes(txdb)
head(mm_genes)
length(mm_genes)
mm_exons <- exons(txdb)
length(mm_exons)
chr1_genes <- mm_genes[seqnames(mm_genes) == "chr1"]
head(chr1_genes)
length(chr1_genes$gene_id)
git add analysis.R
fname <-file.choose(airway_cilData.csv)
fname <-file.choose(airway_cilData.csv)
fname <- file.choose()
fname
colData <- read.csv(fname, row.names = 1)
str(colData)
fname <- file.choose()
counts <- read.csv(fname, row.names = 1)
counts <- as.matrix(counts)
dim(counts)
head(counts)
library("SummarizedExoeriment")
library("SummarizedExperiment")
se <- SummarizedExperiment(assay=counts, colData = colData)
dim(assay(se))
subset(se, , dex == "trt")
colData(se)
colSums((assay(se)))
se$lib.size <- colSums(assay(se))
colData(se)
required_pkgs = c("TCGAbiolinks","GEOquery","GenomicDataCommons","limma","curatedTC
GAData","recount“,"curatedMetagenomicData","phyloseq","HMP16SData",
required_pkgs = c(
"TCGAbiolinks",
"GEOquery",
"GenomicDataCommons",
"limma",
"curatedTCGAData",
"recount",
"curatedMetagenomicData",
"phyloseq",
"HMP16SData"
)
BiocManager::install(required_pkgs)
library(GEOquery)
gse = getGEO("GSE103512")[[1]]
se = as(gse, "SummarizedExperiment")
colData(se)
head(colData(se))
names(colData)
names(colData(se))
dim(colData(se))
with(colData(se), table('cancer.type.ch1', 'normal.ch1'))
sds = apply(assay(se),1,sd)
dat = assay(se)[order(sds,decreasing = TRUE)[1:500],]
mdsvals = cmdscale(dist(t(dat)))
mdsvals = as.data.frame(mdsvals)
mdsvals$Type=factor(colData(se)[,'cancer.type.ch1'])
head(mdsvals)
mdsvals$Normal = factor(colData(se)[,'normal.ch1'])
head(mdsvals)
library(ggplot2)
ggplot(mdsvals, aes(x=v1,y=v2,shape = Normal,color=Type)) +
geom_point(alpha=0.6) + theme(text=element_text(size18))
geom_point(alpha=0.6) + theme(text=element_text(size=18))
ggplot(mdsvals, aes(x=v1,y=v2,shape = Normal,color=Type))
ggplot(mdsvals, aes(x=v1,y=v2,shape = Normal,color=Type)) +
geom_point(alpha=0.6) + theme(text=element_text(size=18))
head(mdsvals)
str(mdsvals)
ggplot(mdsvals, aes(x=V1,y=V2,shape = Normal,color=Type)) +
geom_point(alpha=0.6) + theme(text=element_text(size=18))
BiocManager::install("RNAseq123")
dir.create("RNAseq")
setwd("RNAseq")
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar"
,mode="wb")
utils
utils::untar("GSE63310_RAW.tar", exdir=".")
files < c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt",
"GSM1545538_purep53.txt",
"GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3. txt",
"GSM1545541_JMS8-4.txt",
"GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7C.txt"
"GSM1545545_JMS9-P8c.txt")
files <- c(
"GSM1545535_10_6_5_11.txt",
"GSM1545536_9_6_5_11.txt",
"GSM1545538_purep53.txt",
"GSM1545539_JMS8-2.txt",
"GSM1545540_JMS8-3.txt", # 공백 제거
"GSM1545541_JMS8-4.txt",
"GSM1545542_JMS8-5.txt",
"GSM1545544_JMS9-P7C.txt", # 쉼표 추가
"GSM1545545_JMS9-P8c.txt"
)
for (i in paste(files, ".gz", seq=""))
R.utils::gunzip(i, overwrite=TRUE)
files <- c(
"GSM1545535_10_6_5_11.txt",
"GSM1545536_9_6_5_11.txt",
"GSM1545538_purep53.txt",
"GSM1545539_JMS8-2.txt",
"GSM1545540_JMS8-3.txt",
"GSM1545541_JMS8-4.txt",
"GSM1545542_JMS8-5.txt",
"GSM1545544_JMS9-P7C.txt",
"GSM1545545_JMS9-P8c.txt"
)
gz_files <- paste0(files, ".gz")
for (i in gz_files) {
if (file.exists(i)) {
R.utils::gunzip(i, overwrite = TRUE)
} else {
cat("File not found:", i, "\n")
}
}
read.delim(file.path(".", files[1]), nrow=5)
library(RNAseq123)
x <- readDGE(file.path(".", files), column=c(1,3))
class(x)
str(x)
samplenames <- substring(colnames(x), 12, nchar(columes(x)))
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames
colnames(x)<-samplenames
group<-as.factor(c("LP","ML","Basal","Basal","ML","LP","Basal","ML","LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004", "L006", "L008"),c(3,4,2)))
x$samples$lane<-lane
x$samples
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns = c("SYMBOL", "TXCHROM"), keytype = "ENTREZID")
head(genes)
x$genes <- genes
str(x)
savehistory("~/Medical_Genetics/RNAseq/GES.Rhistory")
install.packages("dplyr", dependencies = TRUE)
library(gapminder)
library(dplyr)
glimpse(gapminder)
korea_data <- gapminder %>% filter(country == "Korea, Rep.")
print(korea_data)
ggplot(korea_data, aes(x = year, y = lifeExp)) + geom_line() + geom_point() + ggtitle("대한민국의 기대수명 변화") + xlab("연도") + ylab("기대수명")
install.packages("ggplot2")
library(ggplot2)
ggplot(korea_data, aes(x = year, y = lifeExp)) + geom_line() + geom_point() + ggtitle("대한민국의 기대수명 변화") + xlab("연도") + ylab("기대수명")
install.packages("showtext")
library(showtext)
font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()
ggplot(korea_data, aes(x = year, y = lifeExp)) +
geom_line() +
geom_point() +
ggtitle("대한민국의 기대수명 변화") +
xlab("연도") +
ylab("기대수명")
year_20007 <- gapminder %>%> filter(year == 2007)
year_2007 <- gapminder %>% filter(year == 2007)
glimpse(gapminder)
library(gapminder)
library(dplyr)
year_2007 <- gapminder %>% filter(year == 2007)
gglot(year_2007)
ggplot(year_2007, aes(x = gdpPercap, y = lifeExp)) +
geom_point(aes(color = continent, size = pop)) +
scale_x_log10() +
ggtitle("2007년 연도별 GDP와 기대수명") +
xlab("1인당 GDP (로그 스케일)") +
ylab("기대수명")
continent_lifeExp
continent_lifeExp <- gapminder %>% group_by(continent) %>% summerise(mean_lifeExp = mean(lifeExp))
continent_lifeExp <- gapminder %>% group_by(continent) %>% summarise(mean_lifeExp = mean(lifeExp))
continent_lifeExp
africa_data <- gapminder %>% filter(continent == "Africa")
ggplot(africa_data)
ggplot(africa_data, aes(x=year, y= lifeExp, color = country)) + geom_line() + ggtitle("아프리카 대륙의 기대 수명 변화") + xlab("연도")+ ylab("기대수명")
year_2007
year_2007 <- gapminder %>% filter(year == 2007)
highest_lifeExp <- year_2007 %>% filter(lifeExp == max(lifeExp))
print(highest_lifeExp)
lowest_lifeExp <- year_2007 %>% filter(lifeExp == min(lifeExp))
print(lowest_lifeExp)
continent_gdp <- year_2007 %>% group_by(continent) %>% summarise(mean_gdpPercap = mean(gdpPercap))
continent_gdp
ggplot(continent_gdp, aes(x = continent, y = mean_gdpPercap)) +
geom_bar(stat = "identity", fill = "skyblue") +
ggtitle("2007년 대륙별 평균 GDP") +
xlab("대륙") +
ylab("평균 1인당 GDP")
highest_lifeExp <- gapminder %>% filter(lifeExp == max(lifeExp))
huest_lifeExp$country
highest_lifeExp$country
highest_lifeExp$gdpPercap
ggplot(year_2007, aes(x=continent, y=lifeExp)) + geom_boxplot(fill = "lightblue") + ggtitle("2007년 대륙별 기대수명의 분포") + xlab("대륙") + ylab("기대수명")
continent_c_count <- gapminder %>% group_by(continent) %>% summarise(c_count = n_distinct(country))
continent_c_count
continent_country<- gapminder %>% distinct(continent, country) %>% arrange(continent, country)
print(continent_country, n=Inf)
view(continent_country)
View(continent_country)
library(IRanges)
ranges<-IRanges(start = c(5,15,25,35,45), end = c(10,20,30,40,50))
width(ranges)
filtered_ranges <- ranges[start(ranges) >= 20]
end(ranges) <- end(ranges) +5
if (!requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")
}
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicRanges")
library(GenomicRanges)
gr <- GRanges(
seqnames = c("chr1", "chr1", "chr2", "chr2"),
ranges = IRanges(start = c(100, 200, 300, 400), end = c(150, 250, 350, 450)),
strand = c("+", "+", "-", "-")
)
BiocManager::install("GenomicRanges", force = TRUE)
BiocManager::install(c("GenomicRanges", "GenomeInfoDb", "IRanges"), force = TRUE)
library(GenomicRanges)
gr <- GRanges(
seqnames = c("chr1", "chr1", "chr2", "chr2"),
ranges = IRanges(start = c(100, 200, 300, 400), end = c(150, 250, 350, 450)),
strand = c("+", "+", "-", "-")
)
print(gr)
chr1_ranges <- gr[seqnames(gr) == "chr1"]
negative_srand_width(gr[stand(gr) == "-"])
negative_strand_width <- width(gr[strand(gr) == "-"])
total_negative_width <- sum(negative_strand_width)
gr$gene_id <- c("geneA", "geneB", "geneC", "geneD")
print(gr)
print(chr1_ranges)
print(total_negative_width)
library(GenomicRanges)
gr_disease <- GRanges(
seqnames = c("chr1", "chr1", "chr2"),
ranges = IRanges(start = c(100, 300, 500), end = c(200, 400, 600)),
strand = c("+", "+", "-")
)
gr2 <- GRanges(
seqnames = c("chr1", "chr1", "chr2", "chr1", "chr1"),
ranges = IRanges(start = c(150, 350, 450, 100, 120), end = c(250, 450, 550, 200, 480)),
strand = c("+", "+", "-", "+", "+")
)
overlap_counts <- countOverlaps(gr_disease, gr2)
print(overlap_counts)
most_overlapped <- gr_disease[overlap_counts == max(overlap_counts)]
print(most_overlapped)
gr_disease[overlap_counts==max(overlap_counts)]
savehistory(file = "my_script.R")
