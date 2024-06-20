
library(tidyverse)
#library(gprofiler2)
#library(UpSetR)
library(biomaRt)
library(clusterProfiler)
library(GO.db)

select <- dplyr::select

sheep112 <- list(name="sheep112", dataset="oarambouillet_gene_ensembl", version="112", kegg="oas")

mouse112 <- list(name="mouse93", dataset="mmusculus_gene_ensembl", version="112", kegg="mmu")


# Allow slow biomaRt queries
options(timeout=3600)


cache <- function(name_parts, value) {
    name <- paste(name_parts,collapse="_")
    dir.create("output/cache", recursive=TRUE, showWarnings=FALSE)
    filename <- paste0("output/cache/",name,".rds")
    if (!file.exists(filename)) {
        saveRDS(value, filename)
    }
    readRDS(filename)
}

get_ensembl <- function(spec) {
    ensembl <- useEnsembl(biomart="genes", version=spec$version)
    ensembl <- useDataset(dataset=spec$dataset, mart=ensembl)
    ensembl
}

#listFilters(get_ensembl()) |> View()
#listAttributes(get_ensembl()) |> View()

get_gene_info <- function(spec) {
    cache(c("gene_info",spec$name), {
        getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', 'gene_biotype', 'entrezgene_id'),
            mart = get_ensembl(spec))
    }) |> 
    mutate(entrezgene_id=as.character(entrezgene_id))
}

get_gene_go <- function(spec) { 
    cache(c("gene_go",spec$name), {
        getBM(attributes = c('ensembl_gene_id', 'go_id'),
            mart = get_ensembl(spec),
            curl = curl::new_handle(timeout_ms=3600000))
    })
}

get_ontologies <- function(spec) {
    gene_info <- get_gene_info(spec)
    gene_go <- get_gene_go(spec)
    
    go_terms <- AnnotationDbi::select(GO.db, keys=keys(GO.db), columns=columns(GO.db))
    
    go_ancestors <- c(as.list(GOBPANCESTOR), as.list(GOMFANCESTOR), as.list(GOCCANCESTOR))
    go_ancestors <- tibble(term=names(go_ancestors), ancestor=go_ancestors) |>
        unnest(ancestor) |>
        filter(ancestor != "all") |>
        bind_rows(tibble(term=keys(GO.db), ancestor=keys(GO.db)))
    #nrow(go_ancestors)
    #nrow(distinct(go_ancestors))
    
    # There seem to be a few missing terms from GO.db. Too bad.
    gene_goall <- gene_go |>
        filter(go_id != "") |>
        inner_join(
            select(go_ancestors, go_id=term, ancestor),
            by="go_id", relationship="many-to-many") |>
        #select(gsid=go_id, gene=ensembl_gene_id)
        select(gsid=ancestor, gene=ensembl_gene_id)
    
    gene2name <- select(gene_info, gene=ensembl_gene_id, name=external_gene_name)
    
    get_gson_go <- function(part) {
        go_terms_wanted <- filter(go_terms, ONTOLOGY==part)
        gson::gson(
            gsid2gene=filter(gene_goall, gsid %in% go_terms_wanted$GOID) |> distinct(),
            gsid2name=select(go_terms_wanted, gsid=GOID, name=TERM),
            gene2name=gene2name)
    }
    
    gson_bp <- get_gson_go("BP")
    gson_mf <- get_gson_go("MF")        
    gson_cc <- get_gson_go("CC")
    
    gson_kegg_entrez <- gson_KEGG(spec$kegg)
    kegg_cats <- clusterProfiler:::kegg_category_data()
    gson_kegg <- gson::gson(
        gsid2gene=
            inner_join(gson_kegg_entrez@gsid2gene,
                    select(gene_info, gene=entrezgene_id, ensembl_gene_id), 
                    by="gene", relationship="many-to-many", na_matches="never") |>
            select(gsid, gene=ensembl_gene_id) |>
            distinct(),
        gsid2name=
            gson_kegg_entrez@gsid2name |>
            left_join(mutate(kegg_cats,gsid=paste0(spec$kegg,id)),by="gsid") |>
            mutate(name=paste0(name.x," (",category,")")) |>
            select(gsid,name),
        gene2name=gene2name)
    
    # TODO: reactome
    # - Seems to have sheep data, no matching ENSEMBL ids.
    
    gson_list <- list(
        "GO:BP" = gson_bp,
        "GO:MF" = gson_mf,
        "GO:CC" = gson_cc,
        "KEGG"  = gson_kegg)
    
    gson_list
}


gson_to_matrix <- function(gson, genes, min_size=10) {
    mat <- split(gson@gsid2gene$gene, gson@gsid2gene$gsid) |>
        map(\(set) genes %in% set)
    mat <- do.call(rbind, mat)
    sizes <- rowSums(mat)
    keep <- sizes >= min_size
    mat <- mat[keep,,drop=FALSE]
    sizes <- sizes[keep]
    
    mat <- mat - rowMeans(mat)
    mat <- mat / sqrt(rowSums(mat*mat))
    
    info <- tibble(
            gsid = rownames(mat),
            size = sizes) |>
        left_join(gson@gsid2name, by="gsid")
    
    list(
        matrix=mat,
        info=info)
}