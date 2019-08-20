# get Aravo data set ----------------------------------------------------------
data("aravo",  package = "ade4")
Y <- aravo$spe
SLA <- aravo$traits$SLA
Snow <- aravo$env$Snow
nrepet <- 19 # change to e.g. 499 or 999
result <- CWMSNC_regressions(Snow, Y, SLA, weighing = "N2", nrepet = nrepet)
names(result)
result$p_values
summary(result)
plot(result)


Snow <- aravo$env$Snow
Spread <- log(aravo$traits$Spread)
result <- CWMSNC_regressions(Snow, Y, Spread, weighing = "N2", nrepet = nrepet)
result$p_values
summary(result)
# in the one E vector - one T vector case, you can choose your own labels
plot(result, trait = "log Spread", env = "Snow")

# untransformed analysis of all pairs. See TutorialWA_Aravo_Multi.r
#                                      for transformations that appear useful
#  data frames E,L,T
result <- CWMSNC_regressions(aravo$env, aravo$spe, aravo$traits, weighing = "N2", nrepet = nrepet)
#result$p_values # contains all pairwise p-values
#                  for site-based, species-based and max-based permutations
#result$wFC      # contains all pairwise weighted fourth-corner correlations
# (site-based, species-based, and signed min-based)
summary(result, type = "max", p_value_adjust_method = "fdr", significance_level = 0.05)
# in the  E or T matrix or data frame case,
# names shoud refer to names of variables or labels of factors
# The next statement generates the available valid names.
(nam.list <- plot(result, trait = "", env = "" ))
plot(result, trait = "Spread", env = "Snow" )
## All plots
# for (trait in nam.list$trait.names){
#   for (env in nam.list$env.names){
#     print(plot(result, trait = trait, env = env ))
#   }
# }

