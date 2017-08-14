#######################################
#Wood Rat Survival by Genotype        #
#E. Hunter, November 2016		  #
#######################################

rm(list=ls())

ELIZABETHOFFICE = T
ELIZABETH = F
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

library(lubridate)
library(rgeos)
library(rgdal)
library(raster)
library(tidyr)
library(R2jags)
library(coda)
library(spatstat)

###################################
##########   READ IN DATA

setwd(DataDir)
rat <- read.csv("CR_all_12dec2016.csv", header=TRUE)

	#head(rat)
	#nrow(rat) #13825 observations

effort <- read.csv("CR_TrapsAll_11dec2016.csv", header=TRUE)
effort$DataType <- as.character(effort$DataType)
effort$DATE <- dmy(effort$Date)

############
#  Process the capture and effort data 

rat$DATE <- dmy(rat$Date1900)
rat <- rat[order(rat$DATE),]
rat$qFu <- as.numeric(as.character(rat$qFu))
##REMOVE individuals X-422, X423, X424, X425 from analysis (only caught in last trapping period)
##REMOVE individuals: 830, 831, 457, 649R, unk-22may-AduF-dead, unk-2jul09  (only caught once in a trap that is not in effort file)
rat <- rat[rat$PJM.ID!="X-422" & rat$PJM.ID!="X-423" & rat$PJM.ID!="X-424" & rat$PJM.ID!="X-425" & rat$PJM.ID!="830"
	 & rat$PJM.ID!="831" & rat$PJM.ID!="457" & rat$PJM.ID!="649R" & rat$PJM.ID!="unk-22may-AduF-dead" & rat$PJM.ID!="unk-2jul09", ]

#TRAP NAME (Trap Line + Trap) (EDIT: PETER HAS CREATED TRAPID'S WITH CONSISTENT UTMS)
effort$TrapName <- as.character(effort$TrapID)		#paste(effort$Line, effort$Trap)
#Any non-matching UTMs?
trap.names <- unique(effort$TrapName)
prob.houses <- NULL
for (i in trap.names){
	temp <- effort[effort$TrapName==i,]
	UTM1 <- round(temp$UTM_N[1])
	prob.houses <- c(prob.houses, ifelse(any(round(temp$UTM_N) != UTM1), 1, 0))

}
sum(prob.houses>0)
prob.houses <- trap.names[which(prob.houses>0)]
#write.csv(prob.houses, "probhouses.csv")
rat$TrapName <- as.character(rat$Trap.ID)		#paste(rat$Line, rat$Trap)
#Any trap names in capture file not in effort file?
rat.houses <- unique(rat$TrapName)
discrep <- rat.houses %in% trap.names
rat.houses[which(discrep==FALSE)]

#Fix rat houses that don't match
rat$TrapName[rat$TrapName=="TF-17W"] <- "TF-17w"
rat$TrapName[rat$TrapName=="MZ-32D"] <- "MZ-32"
rat$TrapName[rat$TrapName=="SC-16A"] <- "SC-16"
rat$TrapName[rat$TrapName=="SC-156B"] <- "SC-156"


#Are we hitting every zone in every year?
effort$month <- substr(effort$DATE, 1, 7)
coords <- coordinates(effort[,c(14,13)])
spdf <- SpatialPointsDataFrame(coords=coords, data=effort)
#spdf <- spdf[spdf@data$Year==2009,]
#plot(spdf, pch=20, col=gray(0.5))
#spdf <- SpatialPointsDataFrame(coords=coords, data=effort)
#spdf <- spdf[spdf@data$season=="LATE",]
#plot(spdf, pch=8, add=T)

#Within the heavily sampled central area, they sampled every year, only did the extremes in 2008 and 2010
#Within that central area, they were not sampling every area every month
#But they were sampling every area within the central area in both early/late seasons

#Need to remove captures and effort in those extremity areas (at least for now)
effort.coords <- coordinates(effort[,c(14,13)])
effort.spdf <- SpatialPointsDataFrame(coords=effort.coords, data=effort)
#plot(effort.spdf, pch=20, col=gray(0.5), axes=T)
e <- extent(698500, 700500, 3960000, 3963000)
effort.clip <- crop(effort.spdf, e)
effort.c <- effort.clip@data
#plot(effort.clip, add=T)

rat.coords <- coordinates(rat[,c(38,37)])
rat.spdf <- SpatialPointsDataFrame(coords=rat.coords, data=rat)
rat.clip <- crop(rat.spdf, e)
rat.c <- rat.clip@data

	#nrow(rat.c) #13629 observations in central area
	#length(unique(rat.c$PJM.ID)) #1924 unique individuals trapped



####################################
##########  DEFINE PERIODS
setwd(DataDir)

filename <- "PeriodDefinition.csv"
sink(filename)
cat("Period, Start, End, Median
EARLY2008, 2008-05-21, 2009-03-13, 2008-06-01
EARLY2009, 2009-03-14, 2009-06-30, 2009-06-01
LATE2009, 2009-07-01, 2010-04-20, 2009-07-30
EARLY2010, 2010-04-21, 2010-06-30, 2010-06-01
LATE2010, 2010-07-01, 2011-04-12, 2010-07-30
EARLY2011, 2011-04-13, 2011-06-30, 2011-06-01
LATE2011, 2011-07-01, 2012-04-04, 2011-07-30
EARLY2012, 2012-04-05, 2012-06-30, 2012-06-01
LATE2012, 2012-07-01, 2013-05-02, 2012-07-30
EARLY2013, 2013-05-03, 2013-09-22, 2013-06-01
LATE2013, 2013-09-23, 2013-09-26, 2013-09-24")
sink()

setwd(DataDir)
PeriodDef <- read.csv(filename, header=TRUE)

PeriodDef$Start <- as.Date(PeriodDef$Start, format="%Y-%m-%d")
PeriodDef$End <- as.Date(PeriodDef$End, format="%Y-%m-%d")
#PeriodDef$StartTrap <- as.Date(PeriodDef$StartTrap, format="%Y-%m-%d")
#PeriodDef$EndTrap <- as.Date(PeriodDef$EndTrap, format="%Y-%m-%d")
PeriodDef$Median <- as.Date(PeriodDef$Median, format="%Y-%m-%d")

nPeriods <- nrow(PeriodDef)

#Determine period duration
duration <- numeric(nPeriods - 1)
for (i in 2:nPeriods){
	duration[i-1] = as.numeric(difftime(PeriodDef$Median[i], PeriodDef$Median[i-1], units="days"))/365
}

PeriodDef$Duration <- c(duration, NA)

#Period intervals
PeriodDef$Intervals <- interval(PeriodDef$Start, PeriodDef$End)


#######PRECIPITATION INTERVALS
#Classify winter rain for each year as Oct. of previous year through April of current year
#Winter rain is the only thing that affects each period
#NEW: Two precip periods: winter= Oct-Feb rain, spring= Mar-May rain
setwd(DataDir)
w <- read.csv("WeatherData.csv")
w$DATE <- ymd(paste(w$DATE, "-15", sep=""))

setwd(DataDir)
filename1 <- "PeriodDefinition_PrecipWinter.csv"
sink(filename1)
cat("Period,Start,End
EARLY2008, 10/1/2007, 2/28/2008
EARLY2009, 10/1/2008, 2/28/2009
LATE2009, 10/1/2008, 2/28/2009
EARLY2010, 10/1/2009, 2/28/2010
LATE2010, 10/1/2009, 2/28/2010
EARLY2011, 10/1/2010, 2/28/2011
LATE2011, 10/1/2010, 2/28/2011
EARLY2012, 10/1/2011, 2/28/2012
LATE2012, 10/1/2011, 2/28/2012
EARLY2013, 10/1/2012, 2/28/2013
LATE2013, 10/1/2012, 2/28/2013")
sink()

setwd(DataDir)
filename2 <- "PeriodDefinition_PrecipSpring.csv"
sink(filename2)
cat("Period,Start,End
EARLY2008, 3/1/2008, 5/31/2008
EARLY2009, 3/1/2009, 5/31/2009
LATE2009, 3/1/2009, 5/31/2009
EARLY2010, 3/1/2010, 5/31/2010
LATE2010, 3/1/2010, 5/31/2010
EARLY2011, 3/1/2011, 5/31/2011
LATE2011, 3/1/2011, 5/31/2011
EARLY2012, 3/1/2012, 5/31/2012
LATE2012, 3/1/2012, 5/31/2012
EARLY2013, 3/1/2013, 5/31/2013
LATE2013, 3/1/2013, 5/31/2013")
sink()


setwd(DataDir)
PrecipWinterPeriodDef <- read.csv(filename1)
PrecipSpringPeriodDef <- read.csv(filename2)

PrecipWinterPeriodDef$Start <- mdy(PrecipWinterPeriodDef$Start)
PrecipWinterPeriodDef$End <- mdy(PrecipWinterPeriodDef$End)
PrecipWinterPeriodDef$Intervals <- interval(PrecipWinterPeriodDef$Start, PrecipWinterPeriodDef$End)

PrecipSpringPeriodDef$Start <- mdy(PrecipSpringPeriodDef$Start)
PrecipSpringPeriodDef$End <- mdy(PrecipSpringPeriodDef$End)
PrecipSpringPeriodDef$Intervals <- interval(PrecipSpringPeriodDef$Start, PrecipSpringPeriodDef$End)


#Summarize precip by period
precipwinter <- NULL
for (i in 1:nrow(PrecipWinterPeriodDef)){
	temp <- w[w$DATE %within% PrecipWinterPeriodDef$Intervals[i],]
	precipwinter <- c(precipwinter, sum(temp$PRCP))
}

#Summarize precip by period
precipspring <- NULL
for (i in 1:nrow(PrecipSpringPeriodDef)){
	temp <- w[w$DATE %within% PrecipSpringPeriodDef$Intervals[i],]
	precipspring <- c(precipspring, sum(temp$PRCP))
}

#Standardize
precipwinter.std <- (precipwinter-mean(precipwinter))/sd(precipwinter)
precipspring.std <- (precipspring-mean(precipspring))/sd(precipspring)

######MAST COVARIATE
#Data from Walt Koenig, California Acorn Survey
setwd(DataDir)
mast <- read.csv("MastData.csv")
#Remove chrysolepis and kelloggii (uncommon spp)
mast <- mast[mast$SPECIES!="chrysolepi" & mast$SPECIES!="kelloggii",]
#Average across each species for each site/year
mast.av <- aggregate(mast$XN30, by=list(mast$LOCALITY, mast$YEAR), FUN=mean)
#plot(mast.av$x ~ mast.av$Group.2, col=mast.av$Group.1)
#Just run with Pozo data for now (closest to field site)
mast.av.pz <- mast.av[mast.av$Group.1=="Pozo",]
#Years 2007-2012 (looking at mast of 1 year previous)
mast.av.pz <- mast.av.pz[3:8, 2:3]
mast.std <- c(mast.av.pz[1,2], rep(mast.av.pz[2:6,2], each=2))
mast.std <- (mast.std-mean(mast.std)) / sd(mast.std)


####################################
##########  SUMMARY DATA (global: suffix "G")

nobsG <- nrow(rat.c)				#13619 observations

nindG <- length(unique(rat.c$PJM.ID))	#1914 individuals
#nindG <- 100 #Temporarily run fewer individuals

realindG <- sort(unique(rat.c$PJM.ID))	

nyearsG <- length(unique(rat.c$Year))	#6 years

realyearsG <- sort(unique(rat.c$Year))

nperiodsG <- nrow(PeriodDef)			#11 or 113 periods

realperiodsG <- PeriodDef$Period

idxperiodsG <- 1:nperiodsG

perioddurG <- PeriodDef$Duration

nhousesG <- length(unique(effort.c$TrapName))	#1174 houses

realhousesG <- sort(unique(effort.c$TrapName))

realhousesnumG <- 1:nhousesG		#numeric house names

housematching <- as.data.frame(cbind(realhousesG, realhousesnumG))

nquad <- 12  	#Number of spatial quadrants

####################################
##########  SPATIAL CAPTURE RECAPTURE
#Eliminating robust design
#Each individual has the possibility of being captured in each capture period, probability is a function of distance from the traps

#First need to match house "field names" with numeric names
effort.c <- merge(effort.c, housematching, by.x="TrapName", by.y="realhousesG")
effort.c$realhousesnumG <- as.numeric(as.character(effort.c$realhousesnumG))
rat.c <- merge(rat.c, housematching, by.x="TrapName", by.y="realhousesG", sort=FALSE)
rat.c$realhousesnumG <- as.numeric(as.character(rat.c$realhousesnumG))

#Coordinates for houses 
coords <- NULL
for (i in 1:nhousesG){
	temp <- effort.c[effort.c$realhousesnumG==realhousesnumG[i],]
	temp <- temp[1,]
	temp.c <- cbind(temp$realhousesnumG, temp$UTM_E, temp$UTM_N)
	coords <- rbind(coords, temp.c)
}
coords <- as.data.frame(coords)
names(coords) <- c("realhousesnumG", "x", "y")
coords$realhousesnumG <- as.factor(coords$realhousesnumG)

#Create minimum convex polygon around points 
z <- chull(coords$x, coords$y)
mcp <- coords[c(z, z[1]), ]
mcp <- mcp[,2:3]
mcp <- SpatialPolygons(list(Polygons(list(Polygon(mcp)), ID=1)))
mcp.b <- buffer(mcp, width=100)


#Create equally spaced points along the line (using nquad), and lines perpendicular to center line (45 degree line)
st.area <- owin(xrange=mcp.b@bbox[1,], yrange=mcp.b@bbox[2,])
rlinegridNEW <- function (angle = 45, spacing = 0.1, win = owin()) 
{
    win <- as.owin(win)
    width <- diff(win$xrange)
    height <- diff(win$yrange)
    rmax <- sqrt(width^2 + height^2)/2
    xmid <- mean(win$xrange)
    ymid <- mean(win$yrange)
    u <- 1-rmax		#runif(1, min = 0, max = spacing) - rmax #EDITED FROM SPATSTAT version here (not random start point)
    p <- seq(from = u, to = rmax, by = spacing)
    q <- sqrt(rmax^2 - p^2)
    theta <- pi * ((angle - 90)/180)
    co <- cos(theta)
    si <- sin(theta)
    X <- psp(x0 = xmid + p * co + q * si, y0 = ymid + p * si - 
        q * co, x1 = xmid + p * co - q * si, y1 = ymid + p * 
        si + q * co, window = owin(xmid + c(-1, 1) * rmax, ymid + 
        c(-1, 1) * rmax), check = FALSE)
    X <- X[win]
    return(X)
}
x <- rlinegridNEW(win=st.area, angle=325, spacing=2100/nquad) 
lns <- x$ends
#Add in first line to fill in last part of study area (top right bounding box point)
#lns <- rbind(rep(c(mcp.b@bbox[1,2], mcp.b@bbox[2,2]), times=2), lns) 
names(lns) <- c("x", "y", "x", "y")
poly <- NULL
for (i in 1:nquad){
	ply.coord <- Polygon(rbind(lns[i,1:2], lns[i,3:4], lns[i+1, 3:4], lns[i+1,1:2], lns[i,1:2]))
	poly <- c(poly, ply.coord)	
}
polys <- SpatialPolygons(list(Polygons(poly, ID=1)))
plot(mcp)
plot(polys, add=TRUE, col="yellow")
points(coords$x, coords$y)

#Which houses are in which quad?
effort.c$quad <- 0
coords$quad <- 0
for(i in 1:nquad){
	temp <- which(point.in.polygon(coords$x, coords$y, poly[[i]]@coords[,1], poly[[i]]@coords[,2]) > 0)
	for(j in temp){
		effort.c$quad[effort.c$realhousesnumG==j] <- i
		coords$quad[coords$realhousesnumG==j] <- i
}
}

#Number of trap-nights trapped at each QUAD for each period
eff.per <- array(0, dim=c(nquad, nperiodsG))
for (i in 1:nquad){
	temp <- effort.c[effort.c$quad==i,]
	if(nrow(temp)!=0){
	temp$counter <- 1
	for (j in 1:nperiodsG){
		temp2 <- temp[temp$DATE %within% PeriodDef$Intervals[j],]
		eff.per[i,j] <- sum(temp2$counter)
	}}
}

#Which QUADS are open in each period? And how many?
ch <- apply(eff.per, c(1,2), function(i) any(i!=0))
trapped <- apply(ch, 2, function(i) which(i))
trapped.num <- unlist(lapply(trapped, function(i) length(i)))
trapped <- do.call(rbind, lapply(trapped, function(x){
	length(x) <- max(trapped.num)
	x }))

#Distance between quads
#Average x and y coordinates for points within each quad
quad.ctr <- data.frame(x=rep(NA, times=nquad), y=rep(NA, times=nquad))
for (i in 1:nquad){
	temp <- coords[coords$quad==i,]
	if(nrow(temp)!=0){
	quad.ctr[i,1] <- mean(temp$x)
	quad.ctr[i,2] <- mean(temp$y)
	}
}
quad.dist <- as.matrix(dist(quad.ctr))

#Also assign traps to zones:
coords$zone <- ifelse(coords$y < 3960975, 3, ifelse(coords$y > 3961297, 1, 2))

##########################
#######    MAKE CAPTURE HISTORY
caphist <- array(0, dim=c(nindG, nperiodsG))
for (i in 1:nindG){
	temp <- rat.c[rat.c$PJM.ID==realindG[i],]
	temp.date <- temp$DATE
	for (j in 1:nperiodsG){
		if (any (temp.date %within% PeriodDef$Intervals[j]))
		{caphist[i,j] <- 1}
	}
}

#Make a vector of the first (numeric) period an animal was captured
firsts <- numeric(nindG)
for (i in 1:nindG){
	for (j in 1:nperiodsG){
		if(caphist[i,j]==1 & firsts[i]==0) {firsts[i] <- j}
	}
}

#Need to link the quadrants to the rat.c dataset
quadmatching <- coords[,c(1,4)]
rat.c <- merge(rat.c, quadmatching, by="realhousesnumG", sort=FALSE)


#Count of number of times caught in each trap in each period (for BINOMIAL MODEL)
ind.quad.per <- array(0, dim=c(nindG, nquad, nperiodsG))
for (i in 1:nindG){
	temp <- rat.c[rat.c$PJM.ID==realindG[i],]
	temp.date <- temp$DATE
	temp.quad <- temp$quad
	for (j in 1:nperiodsG){
		if (any (temp.date %within% PeriodDef$Intervals[j]))
		{idx <- which(temp.date %within% PeriodDef$Intervals[j])
		for (r in unique(temp.quad[idx])){	
			temp2 <- temp.quad[idx] 
			ind.quad.per[i,r,j] <- length(temp2[temp2==r])
		}
		}
	}
}

#Check to see if any ind.trap.per are greater than eff.per 
problems <- NULL
problems.i <- NULL
for (i in 1:nindG){
	temp <- ind.quad.per[i,,]
	if (any(temp > eff.per)){
	problems <- rbind(problems, which(temp > eff.per, arr.ind=TRUE))
	problems.i <- c(problems.i, i)
	}
}


#Need an array of the first quad each individual was caught at in each period (to get initial x y coordinates)
ind.quad <- array(0, dim=c(nindG, nperiodsG))
for (i in 1:nindG){
	temp <- rat.c[rat.c$PJM.ID==realindG[i],]
	temp.date <- temp$DATE
	temp.quad <- temp$quad
	for (j in 1:nperiodsG){
		if (any (temp.date %within% PeriodDef$Intervals[j]))
		{idx <- which(temp.date %within% PeriodDef$Intervals[j])[1] #First trapped w/i period
		ind.quad[i,j] <- temp.quad[idx]}
	}
}

#Any individuals that have moved between quads? (YES, ~160 individuals when there are 10 quads)
movers <- NULL
for (i in 1:nindG){
	temp <- ind.quad.per[i,,]
	temp.rws <- unique(which(temp!=0, arr.ind=T)[,1])
	if(length(temp.rws)>1)
		{movers <- c(movers, i)}
}


###########
#Create covariates: sex, age, genotypes
#Sex
isMale <- numeric(nindG)
for(i in 1:nindG){
  ndx <- which(rat.c$PJM.ID==realindG[i])
  sexvec <- rat.c$SexFI[ndx]
  sexvec <- sexvec[which(!is.na(sexvec))]
  sex <- as.character(sexvec[length(sexvec)])    # final determined sex is the true sex...
  isMale[i] <- ifelse(length(sex>0), ifelse(sex=="M",1,ifelse(sex=="F",0,NA)), NA)
}
#There are some individuals whose sex cannot be determined, so we'll make a list of which individuals those are
miss.sex <- which(is.na(isMale))


#Age: for now, assume 2 age classes (adult, juvenile)
#if a juvenile is juvenile in early period, have it still be a juvenile in the late period of the same year
#if a juvenile is juvenile in late period, it can become an adult in the following spring
isJuv  <- array(0,dim=c(nindG,nperiodsG))
for(i in 1:nindG){
	temp <- rat.c[rat.c$PJM.ID==realindG[i],]
	temp2 <- temp[temp$DATE %within% PeriodDef$Intervals[firsts[i]],]
      tempage <- as.character(temp2$MnAgeBYr)[1]
      isJuv[i,firsts[i]] <- ifelse(tempage=="J",1,0)
	if(tempage=="J" & firsts[i] %% 2==0) #"early" periods are all even, so this works out great!
		{isJuv[i,firsts[i]+1] <- 1}
}

#Genotype
#5 categories: full hybrid (40-60%), Fu (90-100%), Ma(0-10%)
isFu <- numeric(nindG)
isMa <- numeric(nindG)
isHyb <- numeric(nindG)
for(i in 1:nindG){
	temp <- rat.c[rat.c$PJM.ID==realindG[i],]
	temp.qFu <- temp$qFu[1]
	isFu[i] <- ifelse(temp.qFu > 0.9, 1, ifelse(is.na(temp.qFu), NA, 0)) 
	isMa[i] <- ifelse(temp.qFu < 0.1, 1, ifelse(is.na(temp.qFu), NA, 0)) 
	isHyb[i] <- ifelse(temp.qFu >= 0.1 & temp.qFu <= 0.9, 1, ifelse(is.na(temp.qFu), NA, 0)) 
}


###########
#Proportion of each genotype in each quadrant
temp <- ind.quad.per[which(isHyb==1),,]
temp <- ifelse(temp>0,1,0)
temp <- apply(temp, c(1,2), FUN="sum")
temp <- ifelse(temp>0,1,0)
hyb.sum <- apply(temp, 2, FUN="sum")

temp <- ind.quad.per[which(isFu==1),,]
temp <- ifelse(temp>0,1,0)
temp <- apply(temp, c(1,2), FUN="sum")
temp <- ifelse(temp>0,1,0)
Fu.sum <- apply(temp, 2, FUN="sum")

temp <- ind.quad.per[which(isMa==1),,]
temp <- ifelse(temp>0,1,0)
temp <- apply(temp, c(1,2), FUN="sum")
temp <- ifelse(temp>0,1,0)
Ma.sum <- apply(temp, 2, FUN="sum")

all.sums <- NULL
for (i in 1:nquad){
	all.sums[i] <- sum(hyb.sum[i], Ma.sum[i], Fu.sum[i])  
}

prop.hyb <- hyb.sum / all.sums
prop.Fu <- Fu.sum / all.sums
prop.Ma <- Ma.sum / all.sums
#prop.bcFu <- bcFu.sum / all.sums
#prop.bcMa <- bcMa.sum / all.sums


#Proportion of each genotype in each quadrant in each period
temp <- ind.quad.per[which(isHyb==1),,]
temp <- ifelse(temp>0,1,0)
hyb.per.sum <- apply(temp, c(2,3), FUN="sum")

temp <- ind.quad.per[which(isFu==1),,]
temp <- ifelse(temp>0,1,0)
Fu.per.sum <- apply(temp, c(2,3), FUN="sum")

temp <- ind.quad.per[which(isMa==1),,]
temp <- ifelse(temp>0,1,0)
Ma.per.sum <- apply(temp, c(2,3), FUN="sum")

all.per.sums <- hyb.per.sum + Fu.per.sum + Ma.per.sum

prop.per.hyb <- hyb.per.sum / all.per.sums
prop.per.Fu <- Fu.per.sum / all.per.sums
prop.per.Ma <- Ma.per.sum / all.per.sums

#Interpolate missing values based on trend of other values in the same quadrant 
library(zoo)
temp <- zoo(t(prop.per.hyb))
prop.per.hyb <- na.approx(temp, rule=2)

temp <- zoo(t(prop.per.Fu))
prop.per.Fu <- na.approx(temp, rule=2)

temp <- zoo(t(prop.per.Ma))
prop.per.Ma <- na.approx(temp, rule=2)


setwd(BUGSDir)
#save.image(file="RatSurvivalModelData_20170424.RData")






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
      quadprob[thisquad,toquad] <-  (1-isthisquad[thisquad,toquad]) * (1/(n.quad-1)) #* (beta.quadprob / D2[thisquad,toquad])   #quadprob is function of distance
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

# set initial quad as the one where it was first trapped... (for now)
firstQuad <- numeric(nindG)
i=1
for(i in 1:nindG){
  firstQuad[i] <- which(ind.quad.per[i,,firsts[i]]>0)[1]
}

D2 <- as.matrix(dist(quad.ctr))^2

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

ch = apply(ind.quad, c(1,2), function(i) any(i!=0))
first.last = apply(ch, 1, function(i) range(which(i)))  #First and last periods trapped

Z <- array(0, dim=c(nindG, nperiodsG))  # nindG
for (i in 1:nindG){    # nindG
	Z[i,first.last[1,i]:first.last[2,i]] = 1    #1 when known to be alive, 0 otherwise
}

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

params <- c("phi0", "sigma", "male.prob", "maleEff", "lambda0.base", "p.move", #"beta.quadprob",
	"juvEff", "JuvHybEff", "JuvMaEff", "MaleHybEff", "MaleMaEff", "sd.ind", "maleSigma",
	"HybEff", "MaEff", "precipWFu", "precipSFu", "precipWMa", "precipSMa", "precipWHyb", "precipSHyb",
	"mastEff", "MaNeighEff", "HybNeighEff", "FuNeighEff")  


	Mod <- jags(data=data, inits=inits, parameters.to.save=params, model.file=model.file,
		n.iter=20000, n.burnin=10000, n.thin=5, n.chains=3)		#, bugs.directory=WinBUGSDir, debug=F, codaPkg=TRUE)


cat(paste("SCRIPT FINISHED RUNNING SUCCESSFULLY, ",Sys.time(),", on ",Sys.Date(),"\n",sep=""))

setwd(BUGSDir)
save(Mod,file=sprintf("Woodrat_survival15_grid_precipinteract_sigma_%s.RData",Sys.Date()))
									   

#load(file="Woodrat_survival13_grid_precipinteract_2017-03-10.RData")

Mod.mcmc <- as.mcmc(Mod)
Mod.mcmc.list <- mcmc.list(Mod.mcmc)
heidel.diag(Mod.mcmc.list)
gelman.diag(Mod.mcmc.list)
plot(Mod.mcmc, ask=T)


BUGSlist <- as.data.frame(Mod$BUGSoutput$sims.list)
#