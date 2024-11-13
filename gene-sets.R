
library(tidyverse)
#library(gprofiler2)
#library(UpSetR)
library(biomaRt)
library(clusterProfiler)
library(GO.db)
library(withr)

select <- dplyr::select

sheep112 <- list(name="sheep112", dataset="oarambouillet_gene_ensembl", version="112", kegg="oas")

mouse112 <- list(name="mouse112", dataset="mmusculus_gene_ensembl", version="112", kegg="mmu", msigdb="EH8299")
mouse110 <- list(name="mouse110", dataset="mmusculus_gene_ensembl", version="110", kegg="mmu", msigdb="EH8299")

# Check for new msigdb:
# library(ExperimentHub)
# eh <- ExperimentHub()
# res <- query(eh , 'msigdb')



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

# Might poo in the namespace?
get_msigdb <- function(spec, collection="", subcollection="") {
    cache(c("msigdb",spec$name,collection,gsub(":","",subcollection)), {
        thing <- ExperimentHub::ExperimentHub()[[spec$msigdb]]
        
        if (collection != "")
            thing <- subsetCollection(thing, collection)
        if (subcollection != "")
            thing <- subsetCollection(thing, subcollection=subcollection)
        
        sets <- geneIds(thing)
        all_symbols <- unique(unlist(sets))
        info <- get_gene_info(spec)
        good_ensembl <- info |> 
            select(gene=ensembl_gene_id, name=external_gene_name) |> 
            distinct() |>
            filter(n()==1,.by=gene)
        gene2name <- tibble(name=all_symbols) |> 
            left_join(good_ensembl,by="name") |>
            mutate(gene=ifelse(is.na(gene),name,gene)) |>
            select(gene,name)
        gsid2gene <- sets |>
            enframe("gsid","name") |> 
            unnest("name") |>
            inner_join(gene2name,"name",relationship="many-to-many") |>
            select(gsid,gene)
        gsid2name <- purrr::map_chr(thing, \(i) paste(
                #i|>collectionType()|>bcCategory()|>replace_na(""),
                #i|>collectionType()|>bcSubCategory()|>replace_na(""),
                #i|>setName(),
                ifelse(i@shortDescription=="",setName(i),i@shortDescription),
                sep="/")) |> 
            enframe("gsid","name")
        gson::gson(
            gsid2gene=gsid2gene,
            gsid2name=gsid2name,
            gene2name=gene2name)
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
    
    if (!is.null(spec$msigdb)) {
#        gson_list[["MSigDB"]] <- get_msigdb(spec)
#        gson_list[["MSigDB_H"]] <- get_msigdb(spec, "h")
#        gson_list[["MSigDB_C2"]] <- get_msigdb(spec, "c2")
        gson_list[["MSigDB_REACTOME"]] <- get_msigdb(spec, subcollection="CP:REACTOME")
        gson_list[["MSigDB_WIKIPATHWAYS"]] <- get_msigdb(spec, subcollection="CP:WIKIPATHWAYS")
#        gson_list[["MSigDB_IMMUNESIGDB"]] <- get_msigdb(spec, subcollection="IMMUNESIGDB")
#        gson_list[["MSigDB_C8"]] <- get_msigdb(spec, "c8")
    }
    
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

grouped_description <- function(score, name, round_to=0.1) {
    #fmt <- paste0("%.",decimal_places,"f")
    #score <- sprintf(fmt, score)
    score <- as.character( round(score/round_to)*round_to )
    
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
correlation_enrichment_slow <- function(gson, scores, up=TRUE, min_size=2, p_cutoff=0.05, n_top=20, n_desc=20, n_samples=1000, round_to=0.1, seed=1) {
    #message("mat")
    mat <- gson_to_matrix(gson, names(scores), min_size=min_size, competitive=TRUE)
    
    s <- scores
    if (!up) s <- -s
    s <- s-mean(s)
    s <- s/sqrt(sum(s*s))
    r <- (mat$matrix %*% s)[,1]

#Fixme
#    if (p_cutoff >= 1.0 || n_samples < 1) {
#        cutoff <- -Inf
#    } else {
#        #message("null sample")
#        #TODO: this is actuall K Nearest Neighbors
#        #dist <- with_seed(seed, replicate(n_samples, max(mat$matrix %*% sample(s))))
#        samples <- with_seed(seed, replicate(n_samples, sample(s)))
#        
#        # Alternative method:
#        #message("nn2")
#        #result <- RANN::nn2(mat$matrix, t(samples))
#        #dist <- (result$nn.dists[,1]**2-2) * -0.5
#        
#        #message("null multiply and max")
#        dist <- apply(mat$matrix %*% samples, 2, max)
#        
#        cutoff <- quantile(dist, 1-p_cutoff) |> as.vector()
#    }
    
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
        grouped_description(genes$score, genes$name, round_to)
    })
    
    result <- relocate(result, name, r, size, top_genes)
    
    list(result=result) #, cutoff=cutoff)
}





correlation_enrichment <- function(gson, scores, up=TRUE, min_size=2, rank_based=FALSE, n_top=20, n_desc=20, n_samples=999, round_to=0.1, seed=1) {
    scores <- na.omit(scores)
    
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
    if (rank_based) s <- rank(s)
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
        grouped_description(genes$score, genes$name, round_to)
    })
    
    result <- relocate(result, name, r, size, top_genes)
    
    if (n_samples > 0) {
        dist <- with_seed(seed, replicate(n_samples, get_r(sample(s))))
        dist <- cbind(dist, r) # Important. Makes it conservative and valid!
        #max_dist <- with_seed(seed, replicate(n_samples, max(get_r(sample(s)))))
        max_dist <- apply(dist,2,max)
        
        # PFER = Per-Family Error-Rate = E(V). Sometimes written PFE
        # expected false discoveries at cutoff value r
        # Similar to BLAST E-value
        # You can achieve your desired PFER by choosing rows where the reported PFER is lower than the desired value.
        # PFER >= FWER
        result$FWER <- 0.0
        result$PFER <- 0.0
        result$FDR <- 0.0
        #col_beta <- ceiling((1-fdr_beta)*ncol(dist)) 
        for(i in seq_len(nrow(result))) {
            #active <- !(names(sets) %in% result$gsid[seq_len(i-1)])
            #max_dist <- apply(dist, 2, \(col) max(col[active]))
            #print(sum(active))
            #print(summary(max_dist))
            
            # The p-value at which we would reject the null hypothesis
            # because the true r is among p*(n_samples+1) out of n_samples+1.
            #fwer <- max(fwer, (sum(max_dist >= result$r[i])+1)/(n_samples+1))
            #fwer <- (sum(max_dist >= result$r[i])+1)/(n_samples+1)
            result$FWER[i] <- sum(max_dist >= result$r[i]) / ncol(dist)
            
            result$PFER[i] <- sum(dist >= result$r[i]) / ncol(dist)
            
            # Storey and Tibshirinai (2003) approximation. We conservatively let pi_0 = 1.
            #result$FDR[i] <- result$PFER[i] / i
            each_R <- sort(colSums(dist >= result$r[i]))
            result$FDR[i] <- mean(each_R / pmax(each_R,i))
            
            
            # YB(1999) upper bound, beta=0.05 method for FDR
            # https://doi.org/10.1016/S0378-3758(99)00041-5
            #if (i > R_beta)
            #    result$FDR[i] <- mean(each_R / (each_R + i - R_beta))
            #else
            #    result$FDR[i] <- result$FWER[i]
            # Seems to behave oddly for i=1.
            
            #each_R <- sort(colSums(dist >= result$r[i]))
            #result$FDR[i] <- mean(each_R / pmax(each_R,i))
            
            #More conservatively:
            #R_beta <- each_R[col_beta]
            #print(c(R_beta, mean(each_R)+2*sd(each_R)))
            # Let's 95% of the time under-estimate the number of discoveries
            #result$FDR[i] <- mean(each_R / pmax(1,each_R,i-pmax(0,R_beta-each_R)))
        }
        
        result$FDR <- result$FDR |> rev() |> cummin() |> rev()
    }
    
    result
}



# ClusterProfiler, same interface as correlation_enrichment.
# Scores should be centered on zero!
gsea_enrichment <- function(gson, scores, up=TRUE, min_size=2, n_top=20, n_desc=20, round_to=0.1, seed=1, exponent=1, ...) {
    scores <- na.omit(scores)
    scores <- sort(scores, decreasing=TRUE)
    
    output <- with_seed(seed, clusterProfiler::GSEA(
        scores, 
        exponent=exponent, 
        scoreType=ifelse(up,"pos","neg"), 
        minGSSize=min_size, 
        maxGSSize=Inf, 
        pvalueCutoff=1.0, 
        gson=gson,
        ...))
    
    result <- output@result |> 
        arrange(NES * ifelse(up,-1,1)) |> 
        slice_head(n=n_top) |>
        rownames_to_column("gsid") |>
        select(gsid, name=Description, NES, size=setSize, p.adjust, enrichmentScore, leading_edge, core_enrichment) |>
        as_tibble()
    
    result$top_genes <- map_chr(result$gsid,\(gsid) {
        genes <- filter(gson@gsid2gene,gsid==.env$gsid) |> 
            inner_join(gson@gene2name,by="gene") |>
            inner_join(enframe(scores,"gene","score"), by="gene") |>
            arrange(score * ifelse(up,-1,1)) |>
            slice_head(n=n_desc)
        grouped_description(genes$score, genes$name, round_to)
    })
    
    result$core_genes <- map_chr(result$core_enrichment,\(core) {
        core <- core |> str_split_1("/")
        if (!up) core <- rev(core)
        core_names <- gson@gene2name$name[ match(core, gson@gene2name$gene) ]
        grouped_description(scores[core], core_names, round_to)
    })
    
    result <- select(result, name, enrichmentScore, NES, FDR=p.adjust, size, core_genes, top_genes, gsid)
    
    result
}



gson_using_names <- function(gson) {
    all_names <- unique(gson@gene2name$name)
    mapping <- gson@gene2name$name
    names(mapping) <- gson@gene2name$gene
    
    gson@gene2name <- tibble(gene=all_names, name=all_names)
    gson@gsid2gene <- gson@gsid2gene |>
        mutate(gene=.env$mapping[gene]) |>
        distinct()
    
    gson
}


# Matrix has genes as rows and *something* as columns
matrix_enrichment <- function(mat, gson_list, enricher=correlation_enrichment, round_to=0.1, ...) {
    mat <- as.matrix(mat)
    has_colnames <- !is.null(colnames(mat))
    if (!has_colnames) {
        colnames(mat) <- paste0("S", seq_len(ncol(mat)))
    }
    
    results <- expand_grid(
        score_set=colnames(mat), 
        ontology=names(gson_list), 
        up=c(TRUE,FALSE))
    
    results$result <- furrr::future_pmap(results, .progress=TRUE, \(score_set,ontology,up) {
        enricher(
            gson_list[[ontology]], 
            mat[, score_set], 
            up=up,
            round_to=round_to,
            ...)
    })
    
    list(
        results=results,
        score_sets=mat,
        round_to=round_to)
}


matrix_enrichment_report <- function(result) {
    df <- result$results
    df$html <- pmap(df, \(score_set, ontology, up, result) {
        result$FWER <- NULL
        result$PFER <- NULL
        table <- gt::gt(result) |>
            gt::fmt_number(columns=r, decimals=3) |>
            gt::fmt_number(columns=FDR, decimals=3) |>
            gt::tab_style("font-size:50%;max-width:20em;word-break:break-all;", location=gt::cells_body(columns=gsid))
        
        bslib::nav_panel(
            title = paste(ontology, ifelse(up,"up","down")), 
            table)
    })
    
    df <- nest(df, .by="score_set")
    df$body <- pmap(df, \(score_set, data) {
        scores <- result$score_sets[,score_set] |> na.omit()
        desc_up <- scores |>
                sort(decreasing=TRUE) |>
                head(20) |>
                f(grouped_description(.,names(.),result$round_to))
        desc_down <- scores |>
                sort() |>
                head(20) |>
                f(grouped_description(.,names(.),result$round_to))
        htmltools::div(
            htmltools::p(htmltools::tag("b", "Top up:"), desc_up),
            htmltools::p(htmltools::tag("b", "Top down:"), desc_down),
            do.call(bslib::navset_tab, data$html)
        )
    })
    df$panel <- pmap(df, \(score_set, data, body) {
        bslib::nav_panel(
            title=score_set,
            body)
    })
    
    if (nrow(df) == 1) {
        return(bslib::page(df$body[[1]]))
    }
    
    do.call(bslib::navset_tab, df$panel)
}