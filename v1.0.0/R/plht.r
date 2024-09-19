#' @title Using PLHT to measure the postnatal development level -- plht(input matrix)
#'
#' @description Please check https://github.com/Ruismart/PLHT for details of the method and code example to plot the result;
#'
#' @param mat a normalized expression matrix:
#'                for bulk RNAseq, should be log2(TPM/CPM/TMM+1);
#'                for singlecell RNAseq, should be assays-'RNA/SCT'-data.
#' @param design A design table of meta.data containing several columns as:
#'                   sample - condition - replicate - tissue - age - batch - ...
#'                   (first three are necessary)
#'                   (design$sample must match colnames(mat))
#' @param stsp To specify standard sample condition in design$condition,
#'                 the final trajectory.score would set 0 for the median(/mean) of this condition.
#' @param stsp.type set "median"(default) or "mean" value for stsp
#'
#' @param scale Whether to do zscore-transformation for the input data matrix.
#'                  For singlecell, recommend to use @scale.data instead, then set scale=FALSE;
#'                  or use another function:plht.seur which is modified for seurat object as input.
#'
#' @param mat.coef.EpC coefficient matrix for PCs(principal components) of Epithelium, built-in data
#' @param mat.coef.Fib coefficient matrix for PCs(principal components) of Fibroblast, built-in data
#'
#' @param features.overlap TRUE/FALSE: TRUE(default) to include overlapped features of EpC/Fib in the calculation.
#'
#' @return A dataframe: sample - Tc.EpC - Tc.Fib  - condition - replicate - (other design columns);
#'             print information for the input data: EpC/Fib features   xxxx, expressed   xxxx
#'
#' @export
#' @examples #
#'
plht <- function(mat,
                 design,
                 stsp,
                 stsp.type="median",
                 scale=TRUE,
                 mat.coef.EpC = mat.coef.EpC,
                 mat.coef.Fib = mat.coef.Fib,
                 features.overlap = TRUE){
  #
  mat <- as.matrix(mat)
  design <- data.frame(design)
  rownames(design) <- design$sample

  #
  zscore_mat <- function(mat){
    Mean <- rowMeans(mat)
    Sd <- apply(mat,1,stats::sd)
    mat <- sweep(mat,1,Mean,FUN='-')
    mat <- sweep(mat,1,Sd,FUN='/')
    return(mat)
  }

  #
  if(scale==TRUE){
    mat <- mat[rowSums(mat)>0,]
    test.z <- zscore_mat(mat)
    test.z <- as.matrix(test.z)
  }else{
    test.z <- mat
  }

  ##   new PC results =  t(scale.data) %*% mat.coef
  ##
  ##   formula.EpC:   Tc = -1.0230*(Y-Y0) - 0.1259*(Z-Z0) + 2.6426    (Y: PC3; Z: PC4)
  ##   formula.Fib:   Tc  =  X - X0    (X: -PC1)

  feature1.tmp <- intersect(rownames(test.z), rownames(mat.coef.EpC))
  feature2.tmp <- intersect(rownames(test.z), rownames(mat.coef.Fib))

  overlap.raw <- intersect(rownames(mat.coef.EpC), rownames(mat.coef.Fib))
  overlap.tmp <- intersect(feature1.tmp, feature2.tmp)

  if(features.overlap==TRUE){
    feature1 <- feature1.tmp
    feature2 <- feature2.tmp
  }else{
    feature1 <- setdiff(feature1.tmp, overlap.raw)
    feature2 <- setdiff(feature2.tmp, overlap.raw)
  }

  overlap.new <- intersect(feature1, feature2)

  PC.EpC <- t(test.z[feature1,]) %*% mat.coef.EpC[feature1,]
  PC.Fib <- t(test.z[feature2,]) %*% mat.coef.Fib[feature2,]

  # Tc.EpC
  #Y0 <- median(PC.EpC[design$sample[design$condition %in% c(stsp)],"PC_3"])
  #Z0 <- median(PC.EpC[design$sample[design$condition %in% c(stsp)],"PC_4"])
  Y0 <- get(stsp.type)(PC.EpC[design$sample[design$condition %in% c(stsp)],"PC_3"])
  Z0 <- get(stsp.type)(PC.EpC[design$sample[design$condition %in% c(stsp)],"PC_4"])

  Tc.EpC = -1.023*(PC.EpC[,"PC_3"]-Y0) - 0.1259*(PC.EpC[,"PC_4"]-Z0) + 2.6426

  # Tc.Fib
  #X0 <- -median(PC.Fib[design$sample[design$condition %in% c(stsp)],"PC_1"])
  X0 <- -get(stsp.type)(PC.Fib[design$sample[design$condition %in% c(stsp)],"PC_1"])

  Tc.Fib <- -PC.Fib[,"PC_1"]-X0

  #
  result <- data.frame(sample = colnames(test.z),
                       Tc.EpC = Tc.EpC,
                       Tc.Fib = Tc.Fib,
                       design[colnames(test.z),grep("sample",colnames(design),invert=T)]
  )
  #
  cat(paste0("plht: \nEpC features  ",dim(mat.coef.EpC)[1],"\nexpressed  ",length(feature1),
             "\nFib features  ",dim(mat.coef.Fib)[1],"\nexpressed  ",length(feature2),

             "\n\noverlapped features  ",length(overlap.raw),
             "\nexpressed  ",length(overlap.tmp),
             "\nused  ",length(overlap.new)))

  return(result)


  # end
}




