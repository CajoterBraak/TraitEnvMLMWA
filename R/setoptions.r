
setoptions4pmax_test <- function(verbose = TRUE, FUN_test_statistics, with.MLM= FALSE, nrepet = 49,
                                 level = 'sites', perm.mat.spe = 0 , perm.mat.sit = 0, print = 10,
                                 falseSize = FALSE, alpha = 0.05, nmax = 6, how_to_permute=list(sites=permute::how(),species=permute::how()), filek = -1){
  # set options for use in PermutationTest_r_c_max
  # perm.mat.spe =0 : random sample; else number of species,
  # if perm.mat.spe or sit >= nmax : random sampling only
  # or if the number of systematic samples;
  # Beware f_allPerms is for one-dim trait and env as it removes the mirrored permutations. Not such a restriction for nmax 6
  # as with only 5 species or sites, more than one trait or env is not worthwile
  if (perm.mat.spe > 0 & perm.mat.spe < nmax) perm.mat.spe <- f_allPerms(perm.mat.spe)
  if (perm.mat.sit > 0 & perm.mat.sit < nmax) perm.mat.sit <- f_allPerms(perm.mat.sit)
  if (filek ==0 ) filek <- sample(1000,1)
  how_to_permute[[1]]$nperm <- nrepet
  how_to_permute[[2]]$nperm <- nrepet
  list(verbose = verbose, FUN_test_statistics = FUN_test_statistics, with.MLM= with.MLM,
       nrepet= nrepet, level = level,
       perm.mat.spe=perm.mat.spe, perm.mat.sit= perm.mat.sit, how_to_permute=how_to_permute,
       print = print,
       falseSize = falseSize, alpha = alpha, filek = filek )
}


setoptions4MLM<- function(   formula = y~ poly(trait,2) + poly(env,2) + trait:env + (1+trait|site)+ (1+env|species),
                             family = "nbinom2", K = 0,
                             model.names= c("MLM3"),
                             library = "glmmTMB", nAGQ=0,
                             estimation = FALSE, invert_p_value = TRUE, output = "short",
                             with.LRT= FALSE, LRT.only = FALSE, verbose = FALSE){
  # LRT: use the LRT test  instead of the Wald test
  # if LRT == FALSE and with.LRT == TRUE then LRT is an extra test statistic; name it in model.names!!!
  if (LRT.only) with.LRT <- TRUE
  list(formula = formula,  family =  family, model.names = model.names, library=library, nAGQ=nAGQ,
       estimation = estimation, invert_p_value = invert_p_value, K= K, output = output,
       with.LRT= with.LRT, LRT.only =LRT.only,  verbose = verbose)
}

setoptions4WA<- function(weighing = 0, level = "both", invert_p_value = TRUE, fast = TRUE, type= 'joint', wsvd=TRUE,  verbose = FALSE){
  list(weighing = weighing, level = level, invert_p_value = invert_p_value, fast = fast, type = type, wsvd= wsvd, verbose = verbose)
}
