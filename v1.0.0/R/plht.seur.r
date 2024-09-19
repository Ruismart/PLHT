#' @title Using PLHT to measure the postnatal development level -- plht.seur(input seurat.object)
#'
#' @description Please check https://github.com/Ruismart/PLHT for details of the method and code example to plot the result;
#'
#' @param seur.obj a well-normalized seurat object containing at least two conditions
#'                     which are associated with postnatal development
#'
#' @param assay To select assays-'RNA/SCT'-scale.data for calculation, default is "RNA"
#'
#' @param condition the condition column in seur.obj@meta.data
#'
#' @param stsp To specify standard sample condition in seur.obj@meta.data$condition,
#'                 the final trajectory.score would set 0 for the median(/mean) of this condition.
#' @param stsp.type set "median"(default) or "mean" value for stsp
#'
#' @param mat.coef.EpC coefficient matrix for PCs(principal components) of Epithelium, built-in data
#' @param mat.coef.Fib coefficient matrix for PCs(principal components) of Fibroblast, built-in data
#'
#' @param features.overlap TRUE/FALSE: TRUE(default) to include overlapped features of EpC/Fib in the calculation.
#'
#' @return same seurat object with Tc.EpC/Tc.Fib added as new columns in meta.data
#'
#' @export
#'
#' @examples #
#'
plht.seur <- function(seur.obj,
                      assay="RNA",
                      condition="condition",
                      stsp="CTL|WT|D0|P0",
                      stsp.type="median",
                      mat.coef.EpC = mat.coef.EpC,
                      mat.coef.Fib = mat.coef.Fib,
                      features.overlap = TRUE){
  #
  # mat <- as.matrix(seur.obj@assays[[assay]]@data)
  test.z <- as.matrix(seur.obj@assays[[assay]]@scale.data)


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

  names.stsp <- rownames(seur.obj@meta.data)[grep(stsp,seur.obj@meta.data[,condition])]

  # Tc.EpC
  Y0 <- get(stsp.type)(PC.EpC[names.stsp,"PC_3"])
  Z0 <- get(stsp.type)(PC.EpC[names.stsp,"PC_4"])

  Tc.EpC = -1.023*(PC.EpC[,"PC_3"]-Y0) - 0.1259*(PC.EpC[,"PC_4"]-Z0) + 2.6426

  # Tc.Fib
  X0 <- -get(stsp.type)(PC.Fib[names.stsp,"PC_1"])

  Tc.Fib <- -PC.Fib[,"PC_1"]-X0

  #
  seur.obj$Tc.EpC <- Tc.EpC
  seur.obj$Tc.Fib <- Tc.Fib


  #
  cat(paste0("plht: \nEpC features  ",dim(mat.coef.EpC)[1],"\nexpressed  ",length(feature1),
             "\nFib features  ",dim(mat.coef.Fib)[1],"\nexpressed  ",length(feature2),

             "\n\noverlapped features  ",length(overlap.raw),
             "\nexpressed  ",length(overlap.tmp),
             "\nused  ",length(overlap.new)))

  return(seur.obj)

}











