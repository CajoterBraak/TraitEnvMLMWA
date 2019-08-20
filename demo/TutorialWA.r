# Aravo data set ----------------------------------------------------------
data("aravo",  package = "ade4")
Y <- aravo$spe
SLA <- aravo$traits$SLA
Snow <- aravo$env$Snow

result <- CWMSNC_regressions(Snow, Y, SLA, weighing = "N2", nrepet = 99)
result$p_values
summary(result)
plot(result)

# untransformed analysis of all pairs. See TutorialWA_Aravo_Multi.r
#                                      for transformations that appear useful
#  data frames E,L,T
result <- CWMSNC_regressions(aravo$env, aravo$spe, aravo$traits, weighing = "N2", nrepet = 99)
#result$p_values # contains all pairwise p-values
#                  for site-based, species-based and max-based permutations
#result$wFC      # contains all pairwise weighted fourth-corner correlations
# (site-based, species-based, and signed min-based)
summary(result, type = "max", p_value_adjust_method = "fdr", significance_level = 0.05)
# in the  E or T matrix or data frame case, you cannot choose your own labels;
# names shoud refer to names of variables or labels of factors
(nam.list <- plot(result, trait = "", env = "" ))
plot(result, trait = "Spread", env = "Snow" )
## All plots
# for (trait in nam.list$trait.names){
#   for (env in nam.list$env.names){
#     print(plot(result, trait = trait, env = env ))
#   }
# }

# get Revisit data ---------------------------------------------------
data("Revisit")
## adapt dataframe dat for WA-based analyses
n_sites <- with(Revisit, nlevels(factor(site)))
n_species <- with(Revisit, nlevels(factor(species)))
species <- Revisit$species[seq(from  = 1, by = n_sites, length.out = n_species)]
sites <- Revisit$site[1:n_sites]
trait <- Revisit$trait[seq(from  = 1, by = n_sites, length.out = n_species)]
env <- Revisit$env[1:n_sites]
Y <- matrix(Revisit$value, nrow = n_sites,ncol = n_species,
            dimnames = list(sites= sites,species=species))
# #Alternative using the function dat4MLM2TE_obj
# TE_obj <- dat4MLM2TE_obj(Revisit)
# trait <- TE_obj$T; env <- TE_obj$E; Y <- TE_obj$L
#  CWM/SNC regressions ------------------------------------------
set.seed(1231)
result <- CWMSNC_regressions(env, Y, trait, weighing = "N2", nrepet = 99)
summary(result)
plot(result)



# Comparison of two versions of the max test------------------------------------
obj <- make_obj_for_traitenv(env,Y, trait, cutoff=0)
set.seed(1231)
# slow
system.time(aa<- WA_p_max(obj, nrepet = 99, fast = FALSE))
round(aa$p_values,4)

set.seed(1231)
# default
system.time(bb<- WA_p_max(obj, nrepet = 99, fast = TRUE))
round(bb$p_values,4)
# illustrating their equality
or1<- order(aa$sim.row)
or2<- order(bb$sim.row)
all.equal(or1,or2)
or1<- order(aa$sim.col)
or2<- order(bb$sim.col)
all.equal(or1,or2)


# fast, all three weighings   --------------------------------------------------

system.time(bb<- WA_p_max(obj, weighing = "N2", nrepet = 99))
round(bb$p_values,5)

system.time(FC<- WA_p_max(obj, weighing = "FC", nrepet = 99))
round(FC$p_values,5)

system.time(cc<- WA_p_max(obj, weighing = "unw", nrepet = 99))
round(cc$p_values,5)

# comparion with fourtcorner in ade4  ---------------------------------------

ade4::fourthcorner(data.frame(env),as.data.frame(Y),data.frame(trait), nrepet = 99)
result <- CWMSNC_regressions(env, Y, trait, weighing = "FC", nrepet = 99)
summary(result,digits = 4)



