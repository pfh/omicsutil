
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
    
    # Flip components so loadings have positive skew
    flips <- ifelse(colSums((loadings %*% rotation) ** 3) < 0, -1, 1)
    rotation <- sweep(rotation, 2, flips, '*')
    
    # Order by mean absolute scores
    reorder <- (scores %*% rotation) |> abs() |> colMeans() |> order() |> rev()
    rotation <- rotation[, reorder, drop=FALSE]
        
    colnames(rotation) <- paste0("VM_",seq_len(ncol(rotation)))
    
    vm_scores <- scores %*% rotation 
    vm_loadings <- loadings %*% rotation
    
    list(n=n, scores=vm_scores, loadings=vm_loadings, pca=pca)
}