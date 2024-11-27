> library(TxDb.Hsapiens.UCSC.hg38.knownGene)
> library(TxDb.Mmusculus.UCSC.mm10.ensGene)
> txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
> mm_genes <- genes(txdb)
> head(mm_genes)
GRanges object with 6 ranges and 1 metadata column:
  seqnames              ranges strand |
  <Rle>           <IRanges>  <Rle> |
  ENSMUSG00000000001     chr3 108107280-108146146      - |
  ENSMUSG00000000003     chrX   77837901-77853623      - |
  ENSMUSG00000000028    chr16   18780447-18811987      - |
  ENSMUSG00000000031     chr7 142575529-142578143      - |
  ENSMUSG00000000037     chrX 161117193-161258213      + |
  ENSMUSG00000000049    chr11 108343354-108414396      + |
  gene_id
<character>
  ENSMUSG00000000001 ENSMUSG00000000001
ENSMUSG00000000003 ENSMUSG00000000003
ENSMUSG00000000028 ENSMUSG00000000028
ENSMUSG00000000031 ENSMUSG00000000031
ENSMUSG00000000037 ENSMUSG00000000037
ENSMUSG00000000049 ENSMUSG00000000049
-------
  seqinfo: 66 sequences (1 circular) from mm10 genome
> length(mm_genes)
[1] 39017
> mm_exons <- exons(txdb)
> length(mm_exons)
[1] 348801
> chr1_genes <- mm_genes[seqnames(mm_genes) == "chr1"]
> head(chr1_genes)
GRanges object with 6 ranges and 1 metadata column:
  seqnames              ranges strand |
  <Rle>           <IRanges>  <Rle> |
  ENSMUSG00000000544     chr1 166130238-166166510      + |
  ENSMUSG00000000817     chr1 161780691-161788495      - |
  ENSMUSG00000001138     chr1   36511867-36528237      + |
  ENSMUSG00000001143     chr1   36422065-36445271      - |
  ENSMUSG00000001305     chr1 186720978-186749358      - |
  ENSMUSG00000001674     chr1 121553835-121567989      - |
  gene_id
<character>
  ENSMUSG00000000544 ENSMUSG00000000544
ENSMUSG00000000817 ENSMUSG00000000817
ENSMUSG00000001138 ENSMUSG00000001138
ENSMUSG00000001143 ENSMUSG00000001143
ENSMUSG00000001305 ENSMUSG00000001305
ENSMUSG00000001674 ENSMUSG00000001674
-------
  seqinfo: 66 sequences (1 circular) from mm10 genome
> length(chr1_genes$gene_id)
[1] 2027