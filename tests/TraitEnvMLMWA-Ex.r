library(TraitEnvMLMWA)
# Aravo data set ----------------------------------------------------------
data("aravo",  package = "ade4")
Y <- aravo$spe
SLA <- aravo$traits$SLA
Snow <- aravo$env$Snow

result <- CWMSNC_regressions(Snow, Y, SLA, weighing = "N2", nrepet = 999)
result$p_values
summary(result)
#plot(result)
plot(result, trait = "SLA", env  = "Snow")

Spread <- log(aravo$traits$Spread)
result_SLA_Spread <- CWMSNC_regressions(Snow, Y, Spread, weighing = "N2", nrepet = 999)
summary(result_SLA_Spread)


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


#  CWM/SNC regressions ------------------------------------------
set.seed(1231)
result <- CWMSNC_regressions(env, Y, trait, weighing = "N2", nrepet = 499)
summary(result)
plot(result)

# Comparison of two versions of the max test------------------------------------
obj <- make_obj_for_traitenv(env,Y, trait, cutoff=0)
set.seed(1231)
# slow
system.time(aa<- WA_p_max(obj, nrepet =499, fast = FALSE))
round(aa$p_values,4)

set.seed(1231)
# default
system.time(bb<- WA_p_max(obj, nrepet =499, fast = TRUE))
round(bb$p_values,4)
# illustrating their equality
or1<- order(aa$sim.row)
or2<- order(bb$sim.row)
all.equal(or1,or2)
or1<- order(aa$sim.col)
or2<- order(bb$sim.col)
all.equal(or1,or2)


# fast, all three weighings   --------------------------------------------------

system.time(bb<- WA_p_max(obj, weighing = "N2", nrepet = 999))
round(bb$p_values,5)

system.time(FC<- WA_p_max(obj, weighing = "FC", nrepet = 999))
round(FC$p_values,5)

system.time(cc<- WA_p_max(obj, weighing = "unw", nrepet = 999))
round(cc$p_values,5)

# comparion with fourtcorner in ade4  ---------------------------------------

ade4::fourthcorner(data.frame(env),as.data.frame(Y),data.frame(trait), nrepet = 999)
result <- CWMSNC_regressions(env, Y, trait, weighing = "FC", nrepet = 999)
summary(result,digits = 4)



