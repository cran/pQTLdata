## ----setup, include=FALSE-----------------------------------------------------
set.seed(0)
knitr::opts_chunk$set(
  out.extra = 'style="display:block; margin: auto"',
  fig.align = "center",
  fig.path = "pQTLdata/",
  collapse = TRUE,
  comment = "#>",
  dev = "png")

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
pkgs <- c("dplyr", "grid", "EnsDb.Hsapiens.v75", "ensembldb", "IRanges", "knitr", "org.Hs.eg.db", "S4Vectors", "VennDiagram")
for (p in pkgs) if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
    if (!requireNamespace(p)) warning(paste0("This vignette needs package `", p, "'; please install"))
}
invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))

## ----caprion, echo=FALSE------------------------------------------------------
accession_list <- c("P04745", "Q9GZN8", "P0C0L5", "Q96N11", "P48960",
                    "Q9NUQ9", "P04062", "P69905", "P62805", "O14745",
                    "P23381", "P54577")
updated_list <- dplyr::filter(pQTLdata::caprion,Accession %in% accession_list) |>
                dplyr::ungroup() |>
                dplyr::select(Gene,Gene.orig,Protein,Accession,ensGenes,chr,start,end)
knitr::kable(updated_list,caption="Updated information on Caprion")

## ----overlap1, echo=FALSE, fig.cap="Overlap of Olink Explore HT, SomaScan V4.1 and Caprion", fig.height=9, fig.width=9, messages=FALSE, results="hide"----
suppressMessages(library(dplyr))
suppressMessages(library(VennDiagram))
uniprot <- list(Olink=pull(pQTLdata::Olink_Explore_HT,UniProt.ID),
                SomaScan=pull(pQTLdata::SomaScanV4.1,UniProt.ID),
                Caprion=pull(pQTLdata::caprion,Accession))
lapply(uniprot,head)
olink_somascan_caprion <- VennDiagram::venn.diagram(uniprot,filename = NULL,disable.logging = TRUE,
                                                    cex = 2.5, cat.cex = 2.5, cat.pos = c(0,0,180),
                                                    height=8,width=8,units="in")
grid.newpage()
grid.draw(olink_somascan_caprion)

## ----overlap2, echo=FALSE, fig.cap="Overlap of Olink Explore HT, SomaScan 11K and Caprion", fig.height=9, fig.width=9, messages=FALSE, results="hide"----
suppressMessages(library(dplyr))
suppressMessages(library(VennDiagram))
uniprot <- list(Olink=pull(pQTLdata::Olink_Explore_HT,UniProt.ID),
                SomaScan=pull(pQTLdata::SomaScan11k,UniProt.ID),
                Caprion=pull(pQTLdata::caprion,Accession))
lapply(uniprot,head)
olink_somascan_caprion <- VennDiagram::venn.diagram(uniprot,filename = NULL,disable.logging = TRUE,
                                                    cex = 2.5, cat.cex = 2.5, cat.pos = c(0,0,180),
                                                    height=8,width=8,units="in")
grid.newpage()
grid.draw(olink_somascan_caprion)

## ----inf1, echo=FALSE---------------------------------------------------------
knitr::kable(pQTLdata::inf1,caption="Olink/inflammation panel")

## ----refdata, warnings=FALSE--------------------------------------------------
ensembldb::metadata(EnsDb.Hsapiens.v75)
ensembldb::keytypes(EnsDb.Hsapiens.v75)
exon_info <- ensembldb::exons(EnsDb.Hsapiens.v75)
gene_info <- ensembldb::genes(EnsDb.Hsapiens.v75)
transcript_info <- ensembldb::transcripts(EnsDb.Hsapiens.v75)
colnames(S4Vectors::mcols(gene_info))
colnames(S4Vectors::mcols(transcript_info))
overlaps <- IRanges::findOverlaps(transcript_info, gene_info)
overlapping_transcripts <- transcript_info[queryHits(overlaps)]
overlapping_genes <- gene_info[subjectHits(overlaps)]
overlap_data <- data.frame(
  transcript_id = mcols(overlapping_transcripts)$tx_id,
  gene_id = S4Vectors::mcols(overlapping_genes)$gene_id,
  gene_name = S4Vectors::mcols(overlapping_genes)$gene_name,
  start = pmax(start(overlapping_transcripts), start(overlapping_genes)),
  end = pmin(end(overlapping_transcripts), end(overlapping_genes))
)
gene_symbols <- c("BRCA1", "TP53")
gene_data <- subset(overlap_data,gene_name%in%gene_symbols,select=-c(gene_id,gene_name))
cols <- c("UNIPROTID","PROTEINID","GENEID","GENENAME","SEQNAME","TXID")
info <- ensembldb::select(EnsDb.Hsapiens.v75, keys = gene_symbols, 
                          columns = cols, keytype = "SYMBOL") |>
        dplyr::left_join(head(gene_data,15),by=c('TXID'='transcript_id')) |>
        subset(!is.na(UNIPROTID)&!is.na(start)&!is.na(end))
knitr::kable(head(info,15),caption="Annotation for BRCA1 and TP53")
keytypes(org.Hs.eg.db)
uniprot_ids <- ensembldb::select(org.Hs.eg.db, keys = gene_symbols, columns = "UNIPROT", keytype = "SYMBOL")

