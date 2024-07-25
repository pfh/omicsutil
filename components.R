

choose_lib_size <- function(counts_vec, target, per) {
    counts_vec <- counts_vec[ counts_vec != 0 ]
    guess <- log2(sum(counts_vec))
    f <- function(x) {
        sum(log1p(counts_vec * (per/2**x))) / log(2) - target
    }
    
    result <- uniroot(f, lower=guess-1, upper=guess+1, extendInt="yes")
    # TODO: check/report convergence
    
    2**result$root
} 

superb_log_counts <- function(counts, per=1e5, in_at_least=1) {
    real_lib_sizes <- colSums(counts)
    mat <- log1p(sweep(counts,2,per/real_lib_sizes,"*"))/log(2)
    target <- mean( colSums(mat) )  # min?
    
    effective_lib_sizes <- map_dbl(seq_len(ncol(counts)), \(i)
        choose_lib_size(counts[,i], target, per))
    
    mat <- log1p(sweep(counts,2,per/effective_lib_sizes,"*"))/log(2)
    keep <- rowSums(mat>=1) >= in_at_least
    table(keep)
    mat <- mat[keep,,drop=FALSE]
    list(
        matrix=mat,
        target_total=target,
        info=tibble(
            sample=colnames(counts), 
            real_lib_size=real_lib_sizes, 
            effective_lib_size=effective_lib_sizes))
}


nice_log_counts <- function(counts, per=1e5, in_at_least=1) {
    lib_sizes <- colSums(counts) * calcNormFactors(counts)
    mat <- log1p(sweep(counts,2,per/lib_sizes,"*"))/log(2)

    keep <- rowSums(mat>=1) >= in_at_least
    table(keep)
    mat <- mat[keep,,drop=FALSE]
    mat
}


nice_varimax <- function(mat, n=2) {
    pca <- prcomp(t(mat))
    
    loadings <- pca$rotation[,1:n,drop=FALSE]
    scores <- pca$x[,1:n,drop=FALSE]
    
    if (n == 1) { 
        rotation <- matrix(1)
    } else {
        rotation <- varimax(loadings, normalize=FALSE)$rotmat
    }
    
    # Flip components so loadings have positive outliers
    temp <- loadings %*% rotation
    flips <- ifelse(colSums(sign(temp)*(abs(temp)**4)) < 0, -1, 1)
    rotation <- sweep(rotation, 2, flips, '*')
    
    # Order by mean absolute scores
    reorder <- (scores %*% rotation) |> abs() |> colMeans() |> order() |> rev()
    rotation <- rotation[, reorder, drop=FALSE]
        
    colnames(rotation) <- paste0("VM_",seq_len(ncol(rotation)))
    
    vm_scores <- scores %*% rotation 
    vm_loadings <- loadings %*% rotation
    
    list(n=n, scores=vm_scores, loadings=vm_loadings, pca=pca)
}