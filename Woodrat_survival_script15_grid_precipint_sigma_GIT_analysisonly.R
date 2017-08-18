#######################################
#Wood Rat Survival by Genotype        #
#E. Hunter, November 2016		  #
#######################################

rm(list=ls())

ELIZABETHOFFICE = F
ELIZABETH = T
KEVINOFFICE = F

###################################
##########  SET WORKING DIRECTORY


if(ELIZABETHOFFICE) rootDir <- "C:\\Users\\ehunter\\Dropbox\\Nevada\\Woodrat Project\\"
if(ELIZABETH) rootDir <- "~/Dropbox/Nevada/Woodrat Project/"
if(KEVINOFFICE) rootDir <- "E:\\Dropbox\\Woodrat Project\\"

ScriptDir <- paste(rootDir,"Rcode",sep="")
DataDir <- paste(rootDir,"Data",sep="")
FiguresDir <- paste(rootDir,"RawFigures",sep="")
BUGSDir <- paste(rootDir,"BUGS",sep="")

setwd(DataDir)

getwd()

##################################
##########  LOAD PACKAGES

library(R2jags)
library(coda)

###################################
##########   READ IN DATA

load(file="WoodratSurvivalData_20170817.RData")

####Data Descriptions

# D2 - [k,k] array; Squared distance between quadrant centers (quandrant center defined as the average x and y coordinates of all houses in a quadrant)
# eff.per - [k, t] array; Number of trap-nights trapped in each quadrant for each period
# firstQuad - [i] array; Initial quadrant where each individual was first trapped
# firsts - [i] array; The first period each individual was captured (when the enter the population)
# ind.quad.per - [i,t,k] array; Capture history of each individual, period, and quadrant
# isFu - [i] array; Logical indicating whether an individual is of the N. fuscipes type
# isHyb - [i] array; Logical indicating whether an individual is of the hybrid type
# isJuv - [i, t] array; Logical indicating whether an individual is a juvenile in a given period
# isMa - [i] array; Logical indicating whether an individual is of the N. macrotis type
# isMale - [i] array; Logical indicating whether an individual is a male
# mast.std - [i] array; standardized covariate of acorn mast abundance from previous year
# nindG - numeric; total number of individuals ever captured
# nperiodsG - numeric; total number of periods
# nquad - numeric; total number of quadrants
# perioddurG - [t-1] array; duration (in fraction of a year) between current period and the previous (from median date to median date)
# precipspring.std - numeric; standardized covariate of spring precipitation
# precipwinter.std - numeric; standardized covariate of winter precipitation
# prop.per.Fu - [t,k] array; relative proportion of N. fuscipes (compared to other types) in each period and quadrant (some values were interpolated when a quadrant was not surveyed, this was done to improve model convergence for movement probabilities and does not affect survival estimates)
# prop.per.hyb - [t,k] array; relative proportion of hybrids (compared to other types) in each period and quadrant (some values were interpolated when a quadrant was not surveyed, this was done to improve model convergence for movement probabilities and does not affect survival estimates)
# prop.per.Ma - [t,k] array; relative proportion of N. macrotis (compared to other types) in each period and quadrant (some values were interpolated when a quadrant was not surveyed, this was done to improve model convergence for movement probabilities and does not affect survival estimates)
# trapped - [t,k] array; which qudrants were sampled in each period
# trapped.num - [t] array; number of quadrants sampled in each period


##########################
#######    CREATE MODEL (JAGS)
setwd(BUGSDir)

filename <- "WoodRatModel_secr_v2_20170814"
model.file <- paste(filename, ".bug", sep="")

sink(model.file)

cat("

# Modified spatial mark recapture model
# E Hunter and K Shoemaker Nov 2016 

model{
  ############
  #OBSERVATION MODEL
  
  for(i in 1:n.ind){
    for(t in (firsts[i]+1):n.periods){
      #Loop only through quads that are open in each period
      for(r in 1:trapped.num[t]){                    
      	ind.quad.per[i,trapped[t,r],t] ~ dbin(lambda[i,t,trapped[t,r]],eff.per[trapped[t,r],t])
      	lambda[i,t,trapped[t,r]] <- lambda0.base * exp(-1*(D2[myQuad[i,t],trapped[t,r]])/(2*sigma2[i])) * alive[i,t] + indeff[i]   #total exposure as function of distance
					#EAH: 4/25 added individual effect to capture prob.
      }
    }
	sigma2[i] <- pow(sigma + maleSigma*isMale[i], 2)
  }
  
  ###########
  #QUAD MODEL

  #probability of moving from quad to new quad
  for(thisquad in 1:n.quad){
    for(toquad in 1:n.quad){
      isthisquad[thisquad,toquad] <- step(step(thisquad-toquad)+step(toquad-thisquad)-2)      # zero probability of moving to the same quad- if you move you really move!
      quadprob[thisquad,toquad] <-  (1-isthisquad[thisquad,toquad]) * (1/(n.quad-1)) 
    }
  }

  for (i in 1:n.ind){
    for(t in (firsts[i]+1):n.periods){
    	myQuad[i,t] <-  move[i,t]*newQuadCand_cat[i,t] + (1-move[i,t])*myQuad[i,t-1]  # which quad does it live in now?   
    	move[i,t] ~ dbern(transitionQuad[i,t])  							# latent variable: did it move to a new quad?
      transitionQuad[i,t] <- pow(p.move, period[t-1])  					# Probability of transitioning to a new quad
      newQuadCand_cat[i,t] ~ dcat(quadprob[myQuad[i,t-1],])   	
    } 
    myQuad[i,firsts[i]] <- firstQuad[i]  #~ dcat(firstquadProb[i,])
  }
  
  ############
  #SURVIVAL MODEL
  #Latent variable: living or dead
  for(i in 1:n.ind){
    for(t in (firsts[i]+1):n.periods){	
    	alive[i,t] ~ dbern(mualive[i,t])	                 # latent variable: is it still alive?
    	mualive[i,t] <- alive[i,(t-1)] * transition[i,t]     # probability of each ind being alive after transitioning
    	transition[i,t] <- pow(phi[i,t-1], period[t-1])      # Probability of transitioning to this period (based on period duration)
    }
  	alive[i,firsts[i]] ~ dbern(1)	#At first capture, individual is alive  
  }
  
  for(i in 1:n.ind){
    for(t in 1:(firsts[i]-1)){
      alive[i,t] ~ dbern(1)
    }
  }
  
  for(i in 1:n.ind){
    for(t in firsts[i]:(n.periods-1)){
    	logit(phi[i,t]) <- mu.phi[i,t]                       # expected survival rate for this individual at this time on the basis of covariates
    	mu.phi[i,t] <- logit.phi + maleEff*isMale[i] + juvEff*isJuv[i,t] + HybEff*isHyb[i] + MaEff*isMa[i] 	
		+ mastEff*mast[t] + precipWHyb*precipW[t]*isHyb[i] + precipSHyb*precipS[t]*isHyb[i]
		+ precipWFu*precipW[t]*isFu[i] + precipSFu*precipS[t]*isFu[i] + precipWMa*precipW[t]*isMa[i] + precipSMa*precipS[t]*isMa[i]
		+ FuNeighEff*Fu.prob[t,myQuad[i,t]]*isFu[i] + MaNeighEff*Ma.prob[t,myQuad[i,t]]*isMa[i] + HybNeighEff*Hyb.prob[t,myQuad[i,t]]*isHyb[i]
		+ MaleHybEff*isMale[i]*isHyb[i] + MaleMaEff*isMale[i]*isMa[i] + JuvHybEff*isJuv[i,t]*isHyb[i] + JuvMaEff*isJuv[i,t]*isMa[i]
		
		
    }
  }
  
  ############
  #PRIORS
  #Survival
  logit.phi <- log(phi0/(1-phi0))#survival on logit scale
  phi0 ~ dunif(0, 1)#mean survival rate (adult, female, Fu genetic type)
  
  juvEff ~ dunif(-4,4)#beta term for juvenile effect on survival
  maleEff ~ dunif(-4,4)#beta term for sex effect on survival
  HybEff ~ dunif(-4,4)#beta terms for effect of genotype on survival
  MaEff ~ dunif(-4,4)

  precipWinterEff ~ dunif(-4,4) #Effect of standardized winter precipitation on survival
  precipSpringEff ~ dunif(-4,4) #Effect of standardized spring precipication on survival
  mastEff ~ dunif(-4,4)   #Effect of standardized mast on survival
  precipWFu ~ dunif(-4,4) #Interaction of winter precipitation and being Fu
  precipSFu ~ dunif(-4,4) #Interaction of spring precipitation and being Fu
  precipWMa ~ dunif(-4,4) #Interaction of winter precipitation and being Ma
  precipSMa ~ dunif(-4,4) #Interaction of spring precipitation and being Ma
  precipWHyb ~ dunif(-4,4) #Interaction of winter precipitation and being hybrid
  precipSHyb ~ dunif(-4,4) #Interaction of spring precipitation and being hybrid

  #Interaction terms
  MaleHybEff ~ dunif(-4,4)
  MaleMaEff ~ dunif(-4,4)
  JuvHybEff ~ dunif(-4,4)
  JuvMaEff ~ dunif(-4,4)

  #Effects of being near like individuals
  HybNeighEff ~ dunif(-4,4)
  MaNeighEff ~ dunif(-4,4)
  FuNeighEff ~ dunif(-4,4)
  
  #Some individuals are not sexed or genotyped, so need priors for those missing data
  for(i in 1:n.ind){
    isMale[i] ~ dbern(male.prob)
    isHyb[i] ~ dbern(Hyb.prob[firsts[i], firstQuad[i]])
    isMa[i] ~ dbern(Ma.prob[firsts[i], firstQuad[i]])
    isFu[i] ~ dbern(Fu.prob[firsts[i], firstQuad[i]])
  }
  
  male.prob ~ dunif(0,1)   # probability of being male (account for unknown sex)
  
  #Detection    
  lambda0.base ~ dunif(0,1)    # mean times resident indiv captured at its given quad per independent survey event
  #sigma2 <- pow(sigma,2)
  sigma ~ dunif(0, 500)   # drop off in detection with distance
  maleSigma ~ dunif(0, 60) #Effect of being male on sigma
  
  #Individual effect on capture prob.
  sd.ind ~ dunif(0.0001, 1)
  prec.ind <- 1/pow(sd.ind,2)
  for (i in 1:n.ind){
	indeff[i] ~ dnorm(0, prec.ind)
  }

  
  #Movement among quads
  p.move ~ dunif(0, 1)


} #END MODEL


", fill=TRUE)

sink()


##########
# PACKAGE DATA FOR JAGS
##########

data <- list(
  D2 = D2,
  ind.quad.per=ind.quad.per, # MAIN DATA
  firsts=firsts, 
  n.ind=nindG, 
  n.periods=nperiodsG, 
  period=perioddurG, 
	eff.per=eff.per,
  isMale = isMale,
	isJuv=isJuv, 
	isHyb=isHyb,
	isFu=isFu,
	isMa=isMa,
	n.quad=nquad, 
  precipW = precipwinter.std,
  precipS = precipspring.std,
  mast = mast.std,
  trapped=trapped,
  trapped.num=trapped.num,
  firstQuad = firstQuad,
  Hyb.prob = prop.per.hyb,
  Fu.prob = prop.per.Fu,
  Ma.prob = prop.per.Ma
)

Z <- array(0, dim=c(nindG, nperiodsG))  
Z[] <- 1

inits <- function() {
		list(
		lambda0.base = runif(1, 0.01, 0.05),
  	sigma = runif(1, 200, 300 ),
  	alive = Z,
		male.prob = runif(1, 0.4, 0.6),
		maleEff = runif(1, -0.5, -0.5),
  	juvEff = runif(1, -1, 1),
  	MaleHybEff = runif(1, -1, 1),
  	MaleMaEff = runif(1, -1, 1),
  	JuvHybEff = runif(1, -1, 1),
  	JuvMaEff = runif(1, -1, 1),
		HybNeighEff = runif(1, -1, 1),
		MaNeighEff = runif(1, -1, 1),
		FuNeighEff = runif(1, -1, 1),
		precipWFu = runif(1, -1, 1),
		precipSFu = runif(1, -1, 1),
		precipWMa = runif(1, -1, 1),
		precipSMa = runif(1, -1, 1),
		precipWHyb = runif(1, -1, 1),
		precipSHyb = runif(1, -1, 1),
		mastEff = runif(1, -1, 1),
		HybEff = runif(1, -1, 1),
		MaEff = runif(1, -1, 1),
		p.move = runif(1,0,0.01),
		sd.ind = runif(1, 0.0001, 0.001),
		maleSigma = runif(1, 0, 10)
		)
}
inits()

params <- c("phi0", "sigma", "male.prob", "maleEff", "lambda0.base", "p.move", 
	"juvEff", "JuvHybEff", "JuvMaEff", "MaleHybEff", "MaleMaEff", "sd.ind", "maleSigma",
	"HybEff", "MaEff", "precipWFu", "precipSFu", "precipWMa", "precipSMa", "precipWHyb", "precipSHyb",
	"mastEff", "MaNeighEff", "HybNeighEff", "FuNeighEff")  


	Mod <- jags(data=data, inits=inits, parameters.to.save=params, model.file=model.file,
		n.iter=20000, n.burnin=10000, n.thin=5, n.chains=3)		


Mod.mcmc <- as.mcmc(Mod)
Mod.mcmc.list <- mcmc.list(Mod.mcmc)
heidel.diag(Mod.mcmc.list)
gelman.diag(Mod.mcmc.list)
plot(Mod.mcmc, ask=T)

#