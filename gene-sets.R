
library(tidyverse)
#library(gprofiler2)
#library(UpSetR)
library(biomaRt)
library(clusterProfiler)
library(GO.db)
library(withr)

select <- dplyr::select

sheep112 <- list(name="sheep112", dataset="oarambouillet_gene_ensembl", version="112", kegg="oas")

mouse112 <- list(name="mouse112", dataset="mmusculus_gene_ensembl", version="112", kegg="mmu")
mouse110 <- list(name="mouse110", dataset="mmusculus_gene_ensembl", version="110", kegg="mmu")


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
    
    gson_kegg_entrez <- cache(c("kegg",spec$name), 
        gson_KEGG(spec$kegg)
    )
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


gson_to_matrix <- function(gson, genes, min_size=10, competitive=FALSE) {
    #message("split")
    #mat <- split(gson@gsid2gene$gene, gson@gsid2gene$gsid) |>
    #    map(\(set) genes %in% set)
    #message("bind")
    #mat <- do.call(rbind, mat)
    
    message("table")
    df <- gson@gsid2gene
    mat <- table(
        factor(df$gsid, unique(df$gsid)),
        factor(df$gene, genes))
    
    message("filter, etc")
    sizes <- rowSums(mat)
    keep <- sizes >= min_size
    mat <- mat[keep,,drop=FALSE]
    sizes <- sizes[keep]
    
    #Competitive version?
    if (competitive) {
        mat <- mat - rowMeans(mat)
    }
    
    mat <- mat / sqrt(rowSums(mat*mat))
    
    info <- tibble(
            gsid = rownames(mat),
            size = sizes) |>
        left_join(gson@gsid2name, by="gsid")
    
    list(
        matrix=mat,
        info=info)
}

grouped_description <- function(score, name, decimal_places=1) {
    fmt <- paste0("%.",decimal_places,"f")
    score <- sprintf(fmt, score)
    
    result <- c()
    current <- ""
    for(i in seq_len(length(score))) {
        if (current != score[i]) {
            result[length(result)+1] <- paste0("(",score[i],")")
            current <- score[i]
        }
        result[length(result)+1] <- name[i]
    }
    
    paste(result, collapse=" ")
}


# TODO: This can be made *much* faster and more memory efficient.
correlation_enrichment_slow <- function(gson, scores, up=TRUE, min_size=2, p_cutoff=0.05, n_top=20, n_desc=20, n_samples=1000, decimal_places=1, seed=1) {
    #message("mat")
    mat <- gson_to_matrix(gson, names(scores), min_size=min_size, competitive=TRUE)
    
    s <- scores
    if (!up) s <- -s
    s <- s-mean(s)
    s <- s/sqrt(sum(s*s))
    r <- (mat$matrix %*% s)[,1]
    
    if (p_cutoff >= 1.0 || n_samples < 1) {
        cutoff <- -Inf
    } else {
        #message("null sample")
        #TODO: this is actuall K Nearest Neighbors
        #dist <- with_seed(seed, replicate(n_samples, max(mat$matrix %*% sample(s))))
        samples <- with_seed(seed, replicate(n_samples, sample(s)))
        
        # Alternative method:
        #message("nn2")
        #result <- RANN::nn2(mat$matrix, t(samples))
        #dist <- (result$nn.dists[,1]**2-2) * -0.5
        
        #message("null multiply and max")
        dist <- apply(mat$matrix %*% samples, 2, max)
        
        cutoff <- quantile(dist, 1-p_cutoff) |> as.vector()
    }
    
    result <- mutate(mat$info, r=.env$r) |> 
        #filter(r >= cutoff) |> 
        arrange(-r) |> 
        slice_head(n=n_top)
    
    #message("desc")
    result$top_genes <- map_chr(result$gsid,\(gsid) {
        genes <- filter(gson@gsid2gene,gsid==.env$gsid) |> 
            left_join(gson@gene2name,by="gene") |>
            left_join(enframe(scores,"gene","score"), by="gene") |>
            arrange(score * ifelse(up,-1,1)) |>
            slice_head(n=n_desc)
        grouped_description(genes$score, genes$name, decimal_places)
    })
    
    result <- relocate(result, name, r, size, top_genes)
    
    list(result=result, cutoff=cutoff)
}





correlation_enrichment <- function(gson, scores, up=TRUE, min_size=2, p_cutoff=0.05, n_top=20, n_desc=20, n_samples=1000, decimal_places=1, seed=1) {
    # scores should be named.
    # Convert to indices into scores vector
    df <- gson@gsid2gene |>
        mutate(gene = match(gene, names(scores))) |>
        filter(!is.na(gene))
    
    # Convert into list of vectors of indices
    # Discard sets that are too small
    sets <- split(df$gene, df$gsid)
    n1 <- map_dbl(sets, length)
    keep <- n1 >= min_size
    sets <- sets[keep]
    n1 <- n1[keep]
    
    s <- scores
    if (!up) s <- -s
    s <- s-mean(s)
    s <- s/sqrt(sum(s*s))
    n <- length(s)
    
    # See correlation_enrichment_notes.txt
    n0 <- n - n1
    b <- 1/sqrt(n1*n1/n0+n1)
    a <- -n1/n0 * b
    scale <- b-a
    
    # Function to calculate correlations
    get_r <- function(s) map_dbl(sets,\(idx) sum(s[idx])) * scale
    
    r <- get_r(s)
    
    if (p_cutoff >= 1.0 || n_samples < 1) {
        cutoff <- -Inf
    } else {
        dist <- with_seed(seed, replicate(n_samples, max(get_r(sample(s)))))
        cutoff <- quantile(dist, 1-p_cutoff) |> as.vector()
    }
    
    result <- tibble(
            gsid = names(sets),
            size = n1,
            r = r) |>
        left_join(gson@gsid2name, by="gsid") |>
        #filter(r >= cutoff) |> 
        arrange(-r) |> 
        slice_head(n=n_top)
    
    result$top_genes <- map_chr(result$gsid,\(gsid) {
        genes <- filter(gson@gsid2gene,gsid==.env$gsid) |> 
            inner_join(gson@gene2name,by="gene") |>
            inner_join(enframe(scores,"gene","score"), by="gene") |>
            arrange(score * ifelse(up,-1,1)) |>
            slice_head(n=n_desc)
        grouped_description(genes$score, genes$name, decimal_places)
    })
    
    result <- relocate(result, name, r, size, top_genes)
    
    list(result=result, cutoff=cutoff)
}

