# get Aravo data set ----------------------------------------------------------
data("aravo",  package = "ade4")
Y <- aravo$spe
# Aravo data set: multi-trait multi-env test -----------------------------------
summary(aravo$env)
#E <- aravo$env[,-c(5,6)] # without Snow, without ZoogD
E <- aravo$env[,-c(5)] # without Snow, without ZoogD
# Slope is positive and very skew; therefore we logtransform, but there are zeroes
# add minimum non-zero value before taking logs
a <- min(E[E[,"Slope"]>0, "Slope"])
E[,"Slope"] <- log(E[,"Slope"]+ a)
summary(E)
T <- aravo$traits #[, -6 ] # without SLA
# all traits are skew; therefore we logtransform
T[,c(1,3)] <- log(T[,c(1,3)]+ 1) # traits with zeroes
T[,-c(1,3)] <- log(T[,-c(1,3)]) # traits without zeroes
summary(T)

set.seed(1231)
# score_test = TRUE gives the score test which
# takes intra-trait and intra-env correlations into account
#          If weighing == "FC", it give the test used
#              in double constrained correspondence analysis  )
# score_test = FALSE ignores such correlations.
#              If weighing == "FC", it give the gives the test used in RLQ

E<- aravo$env[,c(1,6,2)]
T <- aravo$traits[,c(2,6)]
nrepet <- 99 # change to e.g. 499 or 999
cc<- WA_p_max(E, Y, T, nrepet =nrepet, weighing = "N2", type = "onebyone")
cc$p_values



dd<- CWMSNC_regressions(E, Y, T, nrepet =nrepet, weighing = "N2")


aa<- WA_p_max(E, Y, T, nrepet =nrepet, weighing = "N2", score_test = TRUE)
aa$p_values

# Explicit conversion to dummy variables as performed in WA_p_max  -------
# convert any factor into dummy variables
E <-  model.matrix(as.formula(paste("~ 0", paste(names(E),  collapse= "+"), sep= "+")), data = E)
summary(E)
# convert any factor into dummy variables
T <-  model.matrix(as.formula(paste("~ 0", paste(names(T),  collapse= "+"), sep= "+")), data = T)
summary(T)
obj <- make_obj_for_traitenv(E,Y, T, cutoff=0)
set.seed(1231)

bb<- WA_p_max(obj, nrepet =nrepet, weighing = "N2", score_test = TRUE)
bb$p_values









