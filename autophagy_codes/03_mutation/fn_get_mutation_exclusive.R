#####################################
# input matrix .data to get mutation pattern of genes in this data.
# by using cometExactTest, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4531541/

library(cometExactTest)

fn_get_mutation_pattern <- function(.data){
  # .data <- mut_test_data
  k = nrow(.data) #for now do test for a pair of genes
  
  #create a grid of binary matrix for k genes
  
  grid.mat = t(expand.grid(rep(list(0:1), k)))

  
  #colllapse grid and get all the levels (all posiible combinations)
  lvls = names(table(apply(grid.mat, 2, paste, collapse = '')))
   
  #convert various Variant_Classification class codes to binary (1 = mutated; 0 = nonmutated)
  .data[.data > 0] = 1

  #colllapse grid and get all the levels (all posiible combinations)
  mat.collapse = data.frame(table(apply(.data, 2, paste, collapse = '')))
  
  #check if for any missing combinations
  lvls.missing = lvls[!lvls %in% mat.collapse[,1]]
  if(length(lvls.missing) > 0){
    mat.collapse = rbind(mat.collapse, data.frame(Var1 = lvls.missing, Freq = 0))
  }
  
  #reorder
  mat.collapse = mat.collapse[order(mat.collapse$Var1),]
  
  #run comet exact test for significance
  pval = cometExactTest::comet_exact_test(tbl = as.numeric(x = mat.collapse$Freq), mutmatplot = FALSE)
  
  pval <- signif(pval, 3) #three decimal points
  
  # export the results
  res <- t(data.frame(value = c(as.numeric(mat.collapse$Freq), pval, 1 - pval)))
  colnames(res) <- c(as.character(mat.collapse$Var1),"exclusive_pvalue", "mutual_pvalue")
  dplyr::as.tbl(as.data.frame(res))
}
