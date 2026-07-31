dnabin_to_phydat_collapsing <- function(ali)
{
  mat=as.character(as.matrix(ali))
  # Collapse non-polymorphic sites
  for (j in 1:ncol(mat)) {
    col=mat[, j]
    bases=col[col %in% c('A','C','G','T','a','c','g','t')]
    uniq=unique(bases)
    if (length(uniq)==1) mat[, j]=uniq
    if (length(uniq)==0) mat[, j]='A'
  }
  pali=phyDat(mat, type = "DNA")
  return(pali)
}
