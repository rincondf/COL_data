# DEGREE DAY CALCULATION AND PHENOLOGY MODELING OF BEET LEAF HOPPER IN 
# WASHINGTON

# PRERAPRED BY DIEGO RINCON (diego.rincon@wsu.edu)


CLEAN_BLH <- CLEAN_BLH[-seq(7187, 7208), ] # coordinates are incorrect

# 1. Formalize dates and calculate julians 

CLEAN_BLH$SURVEY_DATE <- as.Date(CLEAN_BLH$SURVEY_DATE, format = "%m/%d/%Y")
CLEAN_BLH$JULIAN <- as.numeric(format(CLEAN_BLH$SURVEY_DATE, "%j"))
CLEAN_BLH$TRAP_COUNT_RAW <- as.numeric(CLEAN_BLH$TRAP_COUNT_RAW)

# 2. Find unique locations

locations_WA <- subset(CLEAN_BLH, 
                       duplicated(CLEAN_BLH[c("LATITUDE", 
                                              "LONGITUDE", "year")]) == 
                         FALSE)[, c(1, 2, 3, 7)]

write.csv(locations_WA, file = "locationsBLH.csv")

# 3. Extract max and min temp

library(daymetr)

temp_recs <- list()

for(i in 1: length(locations_WA$`FIELD ID`)) {
  temp_recs[[i]] <- download_daymet(
    site = "Daymet",
    lat = locations_WA$LATITUDE[i],
    lon = locations_WA$LONGITUDE[i],
    start = locations_WA$year[i],
    end = locations_WA$year[i],
    path = tempdir(),
    internal = TRUE,
    silent = FALSE,
    force = FALSE,
    simplify = FALSE
  )
}


tmaxWA <- list()

for(i in 1: length(locations_WA$`FIELD ID`)) {
  tmaxWA[[i]] <- as.numeric(temp_recs[[i]]$data[,7])
}

tminWA <- list()

for(i in 1: length(locations_WA$`FIELD ID`)) {
  tminWA[[i]] <- as.numeric(temp_recs[[i]]$data[,8])
}



# 4. Calculate cumulative degree days using retrieved max and min temps 
# and the DAS function for DDs calculation

upT_BLH <- 10 + 19.8
baseT_BLH <- 10

DDs_WA <- list()

for(i in 1:391) {
  DDs_WA[[i]] <- calc_dd_vec(tmax = tmaxWA[[i]], tmin = tminWA[[i]], 
                              lower_threshold = baseT_BLH, 
                              upper_threshold = upT_BLH, 
                              cutoff = "horizontal")
  DDs_WA[[i]] <- cumsum(DDs_WA[[i]])
}


# Matching data collection days

CLEAN_BLH$DDs <- rep(NA, length(CLEAN_BLH$SURVEY_DATE))


for(i in 1: length(locations_WA$`FIELD ID`)) {
  CLEAN_BLH$DDs[which(CLEAN_BLH$LATITUDE == locations_WA$LATITUDE[i] & 
                        CLEAN_BLH$LONGITUDE == locations_WA$LONGITUDE[i])] <-
    DDs_WA[[i]][CLEAN_BLH$JULIAN[which(CLEAN_BLH$LATITUDE == 
                                         locations_WA$LATITUDE[i] & 
                                         CLEAN_BLH$LONGITUDE == 
                                         locations_WA$LONGITUDE[i])]]
}


# aggregate mean captures by julian days

blh_juliansWA <- aggregate(CLEAN_BLH$TRAP_COUNT_RAW, 
                           list(CLEAN_BLH$SURVEY_DATE), mean, na.rm = TRUE)

names(blh_juliansWA) <- c("date", "blh_mean")
blh_juliansWA$Julian <- as.numeric(format(blh_juliansWA$date, "%j"))

plot(blh_juliansWA$Julian, blh_juliansWA$blh_mean, xlab = "Julian days", 
     ylab = "Mean captures")

# aggregate mean captures by degree days

blh_ddsWA <- aggregate(CLEAN_BLH$TRAP_COUNT_RAW, 
                       list(CLEAN_BLH$year, CLEAN_BLH$DDs), 
                       mean, na.rm = TRUE)

names(blh_ddsWA) <- c("year", "dds", "blh_mean")

blh_ddsWA <- blh_ddsWA[-which(blh_ddsWA$blh_mean == "NaN"), ]


plot(blh_ddsWA$dds, blh_ddsWA$blh_mean, xlab = "Degree days", 
     ylab = "Mean captures")

# aggregate mean captures by site

blh_siteWA <- aggregate(CLEAN_BLH[, -c(3, 6, 8)], 
                        list(CLEAN_BLH$LONGITUDE, CLEAN_BLH$LATITUDE, 
                             CLEAN_BLH$SURVEY_DATE), mean, na.rm = TRUE)

blh_siteWA <- blh_siteWA[-which(blh_siteWA$TRAP_COUNT_RAW == "NaN"), -c(1: 3)]

# 5. Calculation of proportions and cumulative proportions captured per site 
# per collection date

yearsWA <- seq(2007, 2023)

# By julian days

# Proportions

blh_juliansWA$props <- rep(NA, length(blh_juliansWA$blh_mean))

for(i in 1: length(yearsWA)) {
  blh_juliansWA$props[which(blh_juliansWA$date > 
                              paste0(yearsWA[i], "-01-01") & 
                              blh_juliansWA$date < 
                              paste0(yearsWA[i], "-12-31"))] <- 
    blh_juliansWA$blh_mean[which(blh_juliansWA$date > 
                                   paste0(yearsWA[i], "-01-01") & 
                                   blh_juliansWA$date < 
                                   paste0(yearsWA[i], "-12-31"))] / 
    sum(blh_juliansWA$blh_mean[which(blh_juliansWA$date > 
                                       paste0(yearsWA[i], "-01-01") & 
                                       blh_juliansWA$date < 
                                       paste0(yearsWA[i], "-12-31"))])
}


plot(blh_juliansWA$Julian, blh_juliansWA$props, xlab = "Julian days", 
     ylab = "Proportion of captures")

points(blh_julians$Julian, blh_julians$props, lwd = 2, col = "red")

# Cumulative proportions


blh_juliansWA$Cum_props <- rep(NA, length(blh_juliansWA$blh_mean))

for(i in 1: length(yearsWA)) {
  blh_juliansWA$Cum_props[which(blh_juliansWA$date > 
                                  paste0(yearsWA[i], "-01-01") & 
                                  blh_juliansWA$date < 
                                  paste0(yearsWA[i], "-12-31"))] <- 
    cumsum(blh_juliansWA$blh_mean[which(blh_juliansWA$date > 
                                          paste0(yearsWA[i], "-01-01") & 
                                          blh_juliansWA$date < 
                                          paste0(yearsWA[i], "-12-31"))]) / 
    sum(blh_juliansWA$blh_mean[which(blh_juliansWA$date > 
                                       paste0(yearsWA[i], "-01-01") & 
                                       blh_juliansWA$date < 
                                       paste0(yearsWA[i], "-12-31"))])
}


plot(blh_juliansWA$Julian, blh_juliansWA$Cum_props, xlab = "Julian days", 
     ylab = "Proportion of captures")

points(blh_julians$Julian, blh_julians$Cum_props, lwd = 2, col = "red")


# By degree days

# Proportions

blh_ddsWA$props <- rep(NA, length(blh_ddsWA$blh_mean))

for(i in 1: length(yearsWA)) {
  blh_ddsWA$props[which(blh_ddsWA$year == yearsWA[i])] <- 
    blh_ddsWA$blh_mean[which(blh_ddsWA$year == yearsWA[i])] / 
    sum(blh_ddsWA$blh_mean[which(blh_ddsWA$year == yearsWA[i])],
        na.rm = TRUE)
}


plot(blh_ddsWA$dds, blh_ddsWA$props, xlab = "Degree days", 
     ylab = "Proportion of captures")
points(blh_dds$dds, blh_dds$props, lwd = 2, col = "red")

# Cumulative proportions


blh_ddsWA$Cum_props <- rep(NA, length(blh_ddsWA$blh_mean))

for(i in 1: length(yearsWA)) {
  blh_ddsWA$Cum_props[which(blh_ddsWA$year == yearsWA[i])] <- 
    cumsum(blh_ddsWA$blh_mean[which(blh_ddsWA$year == yearsWA[i])]) / 
    sum(blh_ddsWA$blh_mean[which(blh_ddsWA$year == yearsWA[i])],
        na.rm = TRUE)
}


plot(blh_ddsWA$dds, blh_ddsWA$Cum_props, xlab = "Degree days", 
     ylab = "Proportion of captures")

points(blh_dds$dds, blh_dds$Cum_props, lwd = 2, col = "red")

# By site

# Proportions

# Unique locations must be updated, because there are some locations 
# with only NAs

locs_WA_Cl <- subset(blh_siteWA, 
                     duplicated(blh_siteWA[c("LATITUDE", 
                                             "LONGITUDE", "year")]) == 
                       FALSE)[, c(1, 2, 3, 5)]

blh_siteWA$propsBLH <- rep(NA, length(blh_siteWA$TRAP_COUNT_RAW))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_siteWA$propsBLH[which(blh_siteWA$LATITUDE == 
                              locs_WA_Cl$LATITUDE[i] &
                              blh_siteWA$LONGITUDE ==
                              locs_WA_Cl$LONGITUDE[i] &
                              blh_siteWA$year == 
                              locs_WA_Cl$year[i])] <-
    blh_siteWA$TRAP_COUNT_RAW[which(blh_siteWA$LATITUDE ==
                                      locs_WA_Cl$LATITUDE[i] & 
                                      blh_siteWA$LONGITUDE ==
                                      locs_WA_Cl$LONGITUDE[i] &
                                      blh_siteWA$year == 
                                      locs_WA_Cl$year[i])] /
    sum(blh_siteWA$TRAP_COUNT_RAW[which(blh_siteWA$LATITUDE ==
                                          locs_WA_Cl$LATITUDE[i] &
                                          blh_siteWA$LONGITUDE ==
                                          locs_WA_Cl$LONGITUDE[i] &
                                          blh_siteWA$year == 
                                          locs_WA_Cl$year[i])])
}


plot(blh_siteWA$DDs, blh_siteWA$propsBLH)
plot(blh_siteWA$JULIAN, blh_siteWA$propsBLH)

# Cumulative proportions

# Checking if retrieved locations and years are in ascending order

# for julian days

a <- matrix(NA, length(locs_WA_Cl$LATITUDE), 35)

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  a[i, 1: length(which(blh_siteWA$LATITUDE == 
                         locs_WA_Cl$LATITUDE[i] &
                         blh_siteWA$LONGITUDE ==
                         locs_WA_Cl$LONGITUDE[i] &
                         blh_siteWA$year == 
                         locs_WA_Cl$year[i]))] <- 
    blh_siteWA$JULIAN[which(blh_siteWA$LATITUDE ==
                              locs_WA_Cl$LATITUDE[i] &
                              blh_siteWA$LONGITUDE ==
                              locs_WA_Cl$LONGITUDE[i] &
                              blh_siteWA$year == 
                              locs_WA_Cl$year[i])]
}

check1 <- rep(NA, length(locs_WA_Cl$`FIELD ID`))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  check1[i] <- any(diff(a[i, ]) < 0, na.rm = TRUE)
}

check1 # they are

# for degree days

b <- matrix(NA, length(locs_WA_Cl$LATITUDE), 35)

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  b[i, 1: length(which(blh_siteWA$LATITUDE == 
                         locs_WA_Cl$LATITUDE[i] & 
                         blh_siteWA$LONGITUDE == 
                         locs_WA_Cl$LONGITUDE[i] &
                         blh_siteWA$year == 
                         locs_WA_Cl$year[i]))] <- 
    blh_siteWA$DDs[which(blh_siteWA$LATITUDE == 
                           locs_WA_Cl$LATITUDE[i] & 
                           blh_siteWA$LONGITUDE == 
                           locs_WA_Cl$LONGITUDE[i] &
                           blh_siteWA$year == 
                           locs_WA_Cl$year[i])]
}

check2 <- rep(NA, length(locs_WA_Cl$LATITUDE))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  check2[i] <- any(diff(b[i, ]) < 0, na.rm = TRUE)
}

check2 # they are



blh_siteWA$CumProps <- rep(NA, length(blh_siteWA$propsBLH))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_siteWA$CumProps[which(blh_siteWA$LATITUDE == 
                              locs_WA_Cl$LATITUDE[i] &
                              blh_siteWA$LONGITUDE ==
                              locs_WA_Cl$LONGITUDE[i] &
                              blh_siteWA$year == 
                              locs_WA_Cl$year[i])] <-
    cumsum(blh_siteWA$TRAP_COUNT_RAW[which(blh_siteWA$LATITUDE ==
                                             locs_WA_Cl$LATITUDE[i] & 
                                             blh_siteWA$LONGITUDE ==
                                             locs_WA_Cl$LONGITUDE[i] &
                                             blh_siteWA$year == 
                                             locs_WA_Cl$year[i])]) /
    sum(blh_siteWA$TRAP_COUNT_RAW[which(blh_siteWA$LATITUDE ==
                                          locs_WA_Cl$LATITUDE[i] &
                                          blh_siteWA$LONGITUDE ==
                                          locs_WA_Cl$LONGITUDE[i] &
                                          blh_siteWA$year == 
                                          locs_WA_Cl$year[i])])
}


plot(blh_siteWA$DDs, blh_siteWA$CumProps)
plot(blh_siteWA$JULIAN, blh_siteWA$CumProps)


plot(blh_siteWA$DDs, blh_siteWA$CumProps)
points(blh_site$BLH_DDs, blh_site$CumPropsBLH2, col = "red", lwd = 2)

plot(blh_siteWA$JULIAN, blh_siteWA$CumProps)
points(blh_site$Julian, blh_site$CumPropsBLH2, col = "red", lwd = 2)


plot(blh_siteWA$DDs, blh_siteWA$propsBLH)
points(blh_site$BLH_DDs, blh_site$propsBLH2, col = "red", lwd = 2)

plot(blh_siteWA$JULIAN, blh_siteWA$propsBLH)
points(blh_site$Julian, blh_site$propsBLH2, col = "red", lwd = 2)



# Truncating WA dataset to 150 < DDs < 920

blh_trncWA <- subset(blh_siteWA, 
                     blh_siteWA$DDs < 920 & blh_siteWA$DDs > 180)[, -c(8, 9)]


blh_trncWA$propsBLH <- rep(NA, length(blh_trncWA$TRAP_COUNT_RAW))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_trncWA$propsBLH[which(blh_trncWA$LATITUDE == 
                              locs_WA_Cl$LATITUDE[i] &
                              blh_trncWA$LONGITUDE ==
                              locs_WA_Cl$LONGITUDE[i] &
                              blh_trncWA$year == 
                              locs_WA_Cl$year[i])] <-
    blh_trncWA$TRAP_COUNT_RAW[which(blh_trncWA$LATITUDE ==
                                      locs_WA_Cl$LATITUDE[i] & 
                                      blh_trncWA$LONGITUDE ==
                                      locs_WA_Cl$LONGITUDE[i] &
                                      blh_trncWA$year == 
                                      locs_WA_Cl$year[i])] /
    sum(blh_trncWA$TRAP_COUNT_RAW[which(blh_trncWA$LATITUDE ==
                                          locs_WA_Cl$LATITUDE[i] &
                                          blh_trncWA$LONGITUDE ==
                                          locs_WA_Cl$LONGITUDE[i] &
                                          blh_trncWA$year == 
                                          locs_WA_Cl$year[i])])
}


blh_trncWA$CumProps <- rep(NA, length(blh_trncWA$propsBLH))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_trncWA$CumProps[which(blh_trncWA$LATITUDE == 
                              locs_WA_Cl$LATITUDE[i] &
                              blh_trncWA$LONGITUDE ==
                              locs_WA_Cl$LONGITUDE[i] &
                              blh_trncWA$year == 
                              locs_WA_Cl$year[i])] <-
    cumsum(blh_trncWA$TRAP_COUNT_RAW[which(blh_trncWA$LATITUDE ==
                                             locs_WA_Cl$LATITUDE[i] & 
                                             blh_trncWA$LONGITUDE ==
                                             locs_WA_Cl$LONGITUDE[i] &
                                             blh_trncWA$year == 
                                             locs_WA_Cl$year[i])]) /
    sum(blh_trncWA$TRAP_COUNT_RAW[which(blh_trncWA$LATITUDE ==
                                          locs_WA_Cl$LATITUDE[i] &
                                          blh_trncWA$LONGITUDE ==
                                          locs_WA_Cl$LONGITUDE[i] &
                                          blh_trncWA$year == 
                                          locs_WA_Cl$year[i])])
}


plot(blh_trncWA$DDs, blh_trncWA$CumProps)
points(blh_site$BLH_DDs, blh_site$CumPropsBLH2, col = "red", lwd = 2)

plot(blh_trncWA$JULIAN, blh_trncWA$CumProps, xlim = c(80, 300))
points(blh_site$Julian, blh_site$CumPropsBLH2, col = "red", lwd = 2)



plot(blh_trncWA$DDs, blh_trncWA$propsBLH)
points(blh_site$BLH_DDs, blh_site$propsBLH2, col = "red", lwd = 2)

plot(blh_trncWA$JULIAN, blh_trncWA$propsBLH, xlim = c(80, 300))
points(blh_site$Julian, blh_site$propsBLH2, col = "red", lwd = 2)


# PARAMETER ESTIMATION


frqsWA <- rep(blh_siteWA$DDs, round(blh_siteWA$propsBLH * 100))
hist(frqsWA)

frqsWA_JD <- rep(blh_siteWA$JULIAN, round(blh_siteWA$propsBLH * 100))
hist(frqsWA_JD)

frqsCO <- rep(blh_site$BLH_DDs[which(!is.na(blh_site$propsBLH2))], 
              round(blh_site$propsBLH2[which(!is.na(blh_site$propsBLH2))] 
                    * 1000))
hist(frqsCO)

frqsWAtr <- rep(blh_trncWA$DDs, round(blh_trncWA$propsBLH * 1000))
hist(frqsWAtr)


# WA

library(mixtools)
library(bbmle)
library(SuppDists)
library(ExtDist)

# Degree days

WAsplit <- gammamixEM(frqsWA, maxit = 5000)

hist(frqsWA, breaks = seq(0, 2000, 10), freq = FALSE)
lines(seq(0, 2000, 0.01), dgamma(seq(0, 2000, 0.01), 
                                 shape = WAsplit$gamma.pars[1, 1],
                                 scale = WAsplit$gamma.pars[2, 1]) / 2, lwd = 2)

lines(seq(0, 2000, 0.01), dgamma(seq(0, 2000, 0.01), 
                                 shape = WAsplit$gamma.pars[1, 2],
                                 scale = WAsplit$gamma.pars[2, 2]) / 2, lwd = 2)


InterFun <- function(x, sh1, sc1, sh2, sc2) {
  (dgamma(x, shape = sh1, scale = sc1) - dgamma(x, shape = sh2, scale = sc2))^2
}


IntsctWA <- optimize(InterFun, interval =  c(500, 1500), 
                     sh1 = WAsplit$gamma.pars[1, 1], 
                     sc1 = WAsplit$gamma.pars[2, 1],
                     sh2 = WAsplit$gamma.pars[1, 2],
                     sc2 = WAsplit$gamma.pars[2, 2])


abline(v = IntsctWA$minimum, lwd = 2, lty = 4)


# Julian days

WAsplit_JD <- gammamixEM(frqsWA_JD, maxit = 5000)

hist(frqsWA_JD, breaks = seq(60, 320, 5), freq = FALSE)
lines(seq(60, 320, 0.1), dgamma(seq(60, 320, 0.1), 
                                 shape = WAsplit_JD$gamma.pars[1, 1],
                                 scale = WAsplit_JD$gamma.pars[2, 1]) / 2, 
      lwd = 2)

lines(seq(60, 320, 0.1), dgamma(seq(60, 320, 0.1), 
                                 shape = WAsplit_JD$gamma.pars[1, 2],
                                 scale = WAsplit_JD$gamma.pars[2, 2]) / 2, 
      lwd = 2)


InterFun <- function(x, sh1, sc1, sh2, sc2) {
  (dgamma(x, shape = sh1, scale = sc1) - dgamma(x, shape = sh2, scale = sc2))^2
}


IntsctWA_JD <- optimize(InterFun, interval =  c(170, 230), 
                        sh1 = WAsplit_JD$gamma.pars[1, 1], 
                        sc1 = WAsplit_JD$gamma.pars[2, 1],
                        sh2 = WAsplit_JD$gamma.pars[1, 2],
                        sc2 = WAsplit_JD$gamma.pars[2, 2])


abline(v = IntsctWA_JD$minimum, lwd = 2, lty = 4)


# Parameter estimation


blh_WA1 <- subset(blh_siteWA, blh_siteWA$DDs < IntsctWA$minimum)[, 1:7]
blh_WA2 <- subset(blh_siteWA, blh_siteWA$DDs > IntsctWA$minimum)[, 1:7]


blh_WA1$props <- rep(NA, length(blh_WA1$TRAP_COUNT_RAW))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_WA1$props[which(blh_WA1$LATITUDE == 
                        locs_WA_Cl$LATITUDE[i] &
                        blh_WA1$LONGITUDE ==
                        locs_WA_Cl$LONGITUDE[i] &
                        blh_WA1$year == 
                        locs_WA_Cl$year[i])] <-
    blh_WA1$TRAP_COUNT_RAW[which(blh_WA1$LATITUDE ==
                                   locs_WA_Cl$LATITUDE[i] & 
                                   blh_WA1$LONGITUDE ==
                                   locs_WA_Cl$LONGITUDE[i] &
                                   blh_WA1$year == 
                                   locs_WA_Cl$year[i])] /
    sum(blh_WA1$TRAP_COUNT_RAW[which(blh_WA1$LATITUDE ==
                                       locs_WA_Cl$LATITUDE[i] &
                                       blh_WA1$LONGITUDE ==
                                       locs_WA_Cl$LONGITUDE[i] &
                                       blh_WA1$year == 
                                       locs_WA_Cl$year[i])])
}


blh_WA1$CumProps <- rep(NA, length(blh_WA1$props))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_WA1$CumProps[which(blh_WA1$LATITUDE == 
                              locs_WA_Cl$LATITUDE[i] &
                              blh_WA1$LONGITUDE ==
                              locs_WA_Cl$LONGITUDE[i] &
                              blh_WA1$year == 
                              locs_WA_Cl$year[i])] <-
    cumsum(blh_WA1$TRAP_COUNT_RAW[which(blh_WA1$LATITUDE ==
                                             locs_WA_Cl$LATITUDE[i] & 
                                             blh_WA1$LONGITUDE ==
                                             locs_WA_Cl$LONGITUDE[i] &
                                             blh_WA1$year == 
                                             locs_WA_Cl$year[i])]) /
    sum(blh_WA1$TRAP_COUNT_RAW[which(blh_WA1$LATITUDE ==
                                          locs_WA_Cl$LATITUDE[i] &
                                          blh_WA1$LONGITUDE ==
                                          locs_WA_Cl$LONGITUDE[i] &
                                          blh_WA1$year == 
                                          locs_WA_Cl$year[i])])
}




blh_WA2$props <- rep(NA, length(blh_WA2$TRAP_COUNT_RAW))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_WA2$props[which(blh_WA2$LATITUDE == 
                        locs_WA_Cl$LATITUDE[i] &
                        blh_WA2$LONGITUDE ==
                        locs_WA_Cl$LONGITUDE[i] &
                        blh_WA2$year == 
                        locs_WA_Cl$year[i])] <-
    blh_WA2$TRAP_COUNT_RAW[which(blh_WA2$LATITUDE ==
                                   locs_WA_Cl$LATITUDE[i] & 
                                   blh_WA2$LONGITUDE ==
                                   locs_WA_Cl$LONGITUDE[i] &
                                   blh_WA2$year == 
                                   locs_WA_Cl$year[i])] /
    sum(blh_WA2$TRAP_COUNT_RAW[which(blh_WA2$LATITUDE ==
                                       locs_WA_Cl$LATITUDE[i] &
                                       blh_WA2$LONGITUDE ==
                                       locs_WA_Cl$LONGITUDE[i] &
                                       blh_WA2$year == 
                                       locs_WA_Cl$year[i])])
}


blh_WA2$CumProps <- rep(NA, length(blh_WA2$props))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_WA2$CumProps[which(blh_WA2$LATITUDE == 
                           locs_WA_Cl$LATITUDE[i] &
                           blh_WA2$LONGITUDE ==
                           locs_WA_Cl$LONGITUDE[i] &
                           blh_WA2$year == 
                           locs_WA_Cl$year[i])] <-
    cumsum(blh_WA2$TRAP_COUNT_RAW[which(blh_WA2$LATITUDE ==
                                          locs_WA_Cl$LATITUDE[i] & 
                                          blh_WA2$LONGITUDE ==
                                          locs_WA_Cl$LONGITUDE[i] &
                                          blh_WA2$year == 
                                          locs_WA_Cl$year[i])]) /
    sum(blh_WA2$TRAP_COUNT_RAW[which(blh_WA2$LATITUDE ==
                                       locs_WA_Cl$LATITUDE[i] &
                                       blh_WA2$LONGITUDE ==
                                       locs_WA_Cl$LONGITUDE[i] &
                                       blh_WA2$year == 
                                       locs_WA_Cl$year[i])])
}


blh_WA2 <- blh_WA2[!is.nan(blh_WA2$props), ]


LL_WA_DDs <- function(gamma, delta, a, b) {
  -sum(pr * log(dJohnsonSB_ab(x = x, 
                              gamma = gamma, 
                              delta = delta, 
                              a = a, b = b) 
                + 0.000001))
}



jns_WA_DDs1 <- mle2(LL_WA_DDs, start = list(gamma = -0.5, 
                                           delta = 2, 
                                           a = min(blh_WA1$DDs) - 1, 
                                           b = max(blh_WA1$DDs) + 1), 
                   data = list(x = blh_WA1$DDs, 
                               pr = blh_WA1$props),
                   method = "Nelder-Mead")
summary(jns_WA_DDs1)

jns_WA1 <- function(x) {
  gamma = coef(jns_WA_DDs1)[1] 
  delta = coef(jns_WA_DDs1)[2] 
  xi = coef(jns_WA_DDs1)[3]
  lambda = coef(jns_WA_DDs1)[4] - coef(jns_WA_DDs1)[3]
  dJohnsonSB(x, gamma = gamma, delta = delta, xi = xi, lambda = lambda) * 100
}

jns_WA1A <- function(x) {
  gamma = coef(jns_WA_DDs1)[1] 
  delta = coef(jns_WA_DDs1)[2] 
  xi = coef(jns_WA_DDs1)[3]
  lambda = coef(jns_WA_DDs1)[4] - coef(jns_WA_DDs1)[3]
  pJohnsonSB(x, gamma = gamma, delta = delta, xi = xi, lambda = lambda)
}

plot(blh_WA1$DDs, blh_WA1$props)
lines(seq(1, 1800), jns_WA1(seq(1, 1800)), col = "blue", lwd = 2)

plot(blh_WA1$DDs, blh_WA1$CumProps)
lines(seq(1, 1800), jns_WA1A(seq(1, 1800)), col = "blue", lwd = 2)





jns_WA_DDs2 <- mle2(LL_WA_DDs, start = list(gamma = -0.5, 
                                            delta = 2, 
                                            a = IntsctWA$minimum - 1, 
                                            b = max(blh_WA2$DDs) + 1), 
                    data = list(x = blh_WA2$DDs, 
                                pr = blh_WA2$props),
                    method = "Nelder-Mead")
summary(jns_WA_DDs2)

jns_WA2 <- function(x) {
  gamma = coef(jns_WA_DDs2)[1] 
  delta = coef(jns_WA_DDs2)[2] 
  xi = coef(jns_WA_DDs2)[3]
  lambda = coef(jns_WA_DDs2)[4] - coef(jns_WA_DDs2)[3]
  dJohnsonSB(x, gamma = gamma, delta = delta, xi = xi, lambda = lambda) * 100
}

jns_WA2A <- function(x) {
  gamma = coef(jns_WA_DDs2)[1] 
  delta = coef(jns_WA_DDs2)[2] 
  xi = coef(jns_WA_DDs2)[3]
  lambda = coef(jns_WA_DDs2)[4] - coef(jns_WA_DDs2)[3]
  pJohnsonSB(x, gamma = gamma, delta = delta, xi = xi, lambda = lambda)
}

plot(blh_WA2$DDs, blh_WA2$props)
lines(seq(1, 1800), jns_WA2(seq(1, 1800)), col = "red", lwd = 2)

plot(blh_WA2$DDs, blh_WA2$CumProps)
lines(seq(1, 1800), jns_WA2A(seq(1, 1800)), col = "blue", lwd = 2)


jns_WA_jnt <- function(x) {
  resp <- rep(NA, length(x))
  for(i in 1: length(x)) {
    if(x[i] < IntsctWA$minimum) {
      resp[i] <- jns_WA1(x[i])
    } else {
      resp[i] <- jns_WA2(x[i])
    }
  }
  
  resp
  
}


jns_WA_jntA <- function(x) {
  resp <- rep(NA, length(x))
  for(i in 1: length(x)) {
    if(x[i] < IntsctWA$minimum) {
      resp[i] <- jns_WA1A(x[i])
    } else {
      resp[i] <- jns_WA2A(x[i])
    }
  }
  
  resp
  
}





jns_WA_jnt2 <- function(x) {
  resp <- rep(NA, length(x))
  for(i in 1: length(x)) {
    resp[i] <- (jns_WA1(x[i]) + jns_WA2(x[i])) / 2
  }
  
  resp
  
}


jns_WA_jnt2A <- function(x) {
  resp <- rep(NA, length(x))
  for(i in 1: length(x)) {
    resp[i] <- (jns_WA1A(x[i]) + jns_WA2A(x[i])) / 2
  }
  
  resp
  
}




plot(blh_WA1$DDs, blh_WA1$props, xlim = c(0, 2000))
points(blh_WA2$DDs, blh_WA2$props)
lines(seq(1, 1800), jns_WA_jnt(seq(1, 1800)), col = "blue", lwd = 2)

plot(blh_WA1$DDs, blh_WA1$CumProps, xlim = c(0, 2000))
points(blh_WA2$DDs, blh_WA2$CumProps)
lines(seq(1, 1800), jns_WA_jntA(seq(1, 1800)), col = "blue", lwd = 2)





plot(blh_siteWA$DDs, blh_siteWA$propsBLH, xlim = c(0, 2000))
lines(seq(1, 1800), jns_WA_jnt2(seq(1, 1800)), col = "blue", lwd = 2)

plot(blh_siteWA$DDs, blh_siteWA$CumProps, xlim = c(0, 2000))
lines(seq(1, 1800), jns_WA_jnt2A(seq(1, 1800)), col = "blue", lwd = 2)


# Julian days


blh_WA1JD <- subset(blh_siteWA, blh_siteWA$JULIAN < IntsctWA_JD$minimum)[, 1:7]
blh_WA2JD <- subset(blh_siteWA, blh_siteWA$JULIAN > IntsctWA_JD$minimum)[, 1:7]

blh_WA1JD$props <- rep(NA, length(blh_WA1JD$TRAP_COUNT_RAW))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_WA1JD$props[which(blh_WA1JD$LATITUDE == 
                        locs_WA_Cl$LATITUDE[i] &
                        blh_WA1JD$LONGITUDE ==
                        locs_WA_Cl$LONGITUDE[i] &
                        blh_WA1JD$year == 
                        locs_WA_Cl$year[i])] <-
    blh_WA1JD$TRAP_COUNT_RAW[which(blh_WA1JD$LATITUDE ==
                                   locs_WA_Cl$LATITUDE[i] & 
                                   blh_WA1JD$LONGITUDE ==
                                   locs_WA_Cl$LONGITUDE[i] &
                                   blh_WA1JD$year == 
                                   locs_WA_Cl$year[i])] /
    sum(blh_WA1JD$TRAP_COUNT_RAW[which(blh_WA1JD$LATITUDE ==
                                       locs_WA_Cl$LATITUDE[i] &
                                       blh_WA1JD$LONGITUDE ==
                                       locs_WA_Cl$LONGITUDE[i] &
                                       blh_WA1JD$year == 
                                       locs_WA_Cl$year[i])])
}



blh_WA1JD$CumProps <- rep(NA, length(blh_WA1JD$props))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_WA1JD$CumProps[which(blh_WA1JD$LATITUDE == 
                           locs_WA_Cl$LATITUDE[i] &
                           blh_WA1JD$LONGITUDE ==
                           locs_WA_Cl$LONGITUDE[i] &
                           blh_WA1JD$year == 
                           locs_WA_Cl$year[i])] <-
    cumsum(blh_WA1JD$TRAP_COUNT_RAW[which(blh_WA1JD$LATITUDE ==
                                          locs_WA_Cl$LATITUDE[i] & 
                                          blh_WA1JD$LONGITUDE ==
                                          locs_WA_Cl$LONGITUDE[i] &
                                          blh_WA1JD$year == 
                                          locs_WA_Cl$year[i])]) /
    sum(blh_WA1JD$TRAP_COUNT_RAW[which(blh_WA1JD$LATITUDE ==
                                       locs_WA_Cl$LATITUDE[i] &
                                       blh_WA1JD$LONGITUDE ==
                                       locs_WA_Cl$LONGITUDE[i] &
                                       blh_WA1JD$year == 
                                       locs_WA_Cl$year[i])])
}




blh_WA2JD$props <- rep(NA, length(blh_WA2JD$TRAP_COUNT_RAW))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_WA2JD$props[which(blh_WA2JD$LATITUDE == 
                        locs_WA_Cl$LATITUDE[i] &
                        blh_WA2JD$LONGITUDE ==
                        locs_WA_Cl$LONGITUDE[i] &
                        blh_WA2JD$year == 
                        locs_WA_Cl$year[i])] <-
    blh_WA2JD$TRAP_COUNT_RAW[which(blh_WA2JD$LATITUDE ==
                                   locs_WA_Cl$LATITUDE[i] & 
                                   blh_WA2JD$LONGITUDE ==
                                   locs_WA_Cl$LONGITUDE[i] &
                                   blh_WA2JD$year == 
                                   locs_WA_Cl$year[i])] /
    sum(blh_WA2JD$TRAP_COUNT_RAW[which(blh_WA2JD$LATITUDE ==
                                       locs_WA_Cl$LATITUDE[i] &
                                       blh_WA2JD$LONGITUDE ==
                                       locs_WA_Cl$LONGITUDE[i] &
                                       blh_WA2JD$year == 
                                       locs_WA_Cl$year[i])])
}


blh_WA2JD$CumProps <- rep(NA, length(blh_WA2JD$props))

for(i in 1: length(locs_WA_Cl$LATITUDE)) {
  blh_WA2JD$CumProps[which(blh_WA2JD$LATITUDE == 
                           locs_WA_Cl$LATITUDE[i] &
                           blh_WA2JD$LONGITUDE ==
                           locs_WA_Cl$LONGITUDE[i] &
                           blh_WA2JD$year == 
                           locs_WA_Cl$year[i])] <-
    cumsum(blh_WA2JD$TRAP_COUNT_RAW[which(blh_WA2JD$LATITUDE ==
                                          locs_WA_Cl$LATITUDE[i] & 
                                          blh_WA2JD$LONGITUDE ==
                                          locs_WA_Cl$LONGITUDE[i] &
                                          blh_WA2JD$year == 
                                          locs_WA_Cl$year[i])]) /
    sum(blh_WA2JD$TRAP_COUNT_RAW[which(blh_WA2JD$LATITUDE ==
                                       locs_WA_Cl$LATITUDE[i] &
                                       blh_WA2JD$LONGITUDE ==
                                       locs_WA_Cl$LONGITUDE[i] &
                                       blh_WA2JD$year == 
                                       locs_WA_Cl$year[i])])
}


blh_WA1JD <- blh_WA1JD[!is.nan(blh_WA1JD$props), ]
blh_WA2JD <- blh_WA2JD[!is.nan(blh_WA2JD$props), ]



LL_WA_JD <- function(gamma, delta, a, b) {
  -sum(pr * log(dJohnsonSB_ab(x = x, 
                              gamma = gamma, 
                              delta = delta, 
                              a = a, b = b) 
                + 0.000001))
}



jns_WA_JD1 <- mle2(LL_WA_JD, start = list(gamma = -0.5, 
                                            delta = 2, 
                                            a = min(blh_WA1JD$JULIAN) - 1, 
                                            b = max(blh_WA1JD$JULIAN) + 1), 
                    data = list(x = blh_WA1JD$JULIAN, 
                                pr = blh_WA1JD$props))
summary(jns_WA_JD1)

jns_WA1JD <- function(x) {
  gamma = coef(jns_WA_JD1)[1] 
  delta = coef(jns_WA_JD1)[2] 
  xi = coef(jns_WA_JD1)[3]
  lambda = coef(jns_WA_JD1)[4] - coef(jns_WA_JD1)[3]
  dJohnsonSB(x, gamma = gamma, delta = delta, xi = xi, lambda = lambda) * 10
}

jns_WA1AJD <- function(x) {
  gamma = coef(jns_WA_JD1)[1] 
  delta = coef(jns_WA_JD1)[2] 
  xi = coef(jns_WA_JD1)[3]
  lambda = coef(jns_WA_JD1)[4] - coef(jns_WA_JD1)[3]
  pJohnsonSB(x, gamma = gamma, delta = delta, xi = xi, lambda = lambda)
}

plot(blh_WA1JD$JULIAN, blh_WA1JD$props)
lines(seq(1, 250), jns_WA1JD(seq(1, 250)), col = "blue", lwd = 2)

plot(blh_WA1JD$JULIAN, blh_WA1JD$CumProps)
lines(seq(1, 250), jns_WA1AJD(seq(1, 250)), col = "blue", lwd = 2)




jns_WA_JD2 <- mle2(LL_WA_JD, start = list(gamma = -0.5, 
                                            delta = 2, 
                                            a = min(blh_WA2JD$JULIAN) - 1, 
                                            b = max(blh_WA2JD$JULIAN) + 1), 
                    data = list(x = blh_WA2$JULIAN, 
                                pr = blh_WA2$props),
                    method = "Nelder-Mead")
summary(jns_WA_JD2)

jns_WA2JD <- function(x) {
  gamma = coef(jns_WA_JD2)[1] 
  delta = coef(jns_WA_JD2)[2] 
  xi = coef(jns_WA_JD2)[3]
  lambda = coef(jns_WA_JD2)[4] - coef(jns_WA_JD2)[3]
  dJohnsonSB(x, gamma = gamma, delta = delta, xi = xi, lambda = lambda) * 10
}

jns_WA2AJD <- function(x) {
  gamma = coef(jns_WA_JD2)[1] 
  delta = coef(jns_WA_JD2)[2] 
  xi = coef(jns_WA_JD2)[3]
  lambda = coef(jns_WA_JD2)[4] - coef(jns_WA_JD2)[3]
  pJohnsonSB(x, gamma = gamma, delta = delta, xi = xi, lambda = lambda)
}

plot(blh_WA2JD$JULIAN, blh_WA2JD$props)
lines(seq(150, 310), jns_WA2JD(seq(150, 310)), col = "red", lwd = 2)

plot(blh_WA2JD$JULIAN, blh_WA2JD$CumProps)
lines(seq(150, 310), jns_WA2AJD(seq(150, 310)), col = "blue", lwd = 2)


jns_WA_jntJD <- function(x) {
  resp <- rep(NA, length(x))
  for(i in 1: length(x)) {
    if(x[i] < IntsctWA_JD$minimum) {
      resp[i] <- jns_WA1JD(x[i])
    } else {
      resp[i] <- jns_WA2JD(x[i])
    }
  }
  
  resp
  
}


jns_WA_jntAJD <- function(x) {
  resp <- rep(NA, length(x))
  for(i in 1: length(x)) {
    if(x[i] < IntsctWA_JD$minimum) {
      resp[i] <- jns_WA1AJD(x[i])
    } else {
      resp[i] <- jns_WA2AJD(x[i])
    }
  }
  
  resp
  
}





jns_WA_jnt2JD <- function(x) {
  resp <- rep(NA, length(x))
  for(i in 1: length(x)) {
    resp[i] <- (jns_WA1JD(x[i]) + jns_WA2JD(x[i])) / 2
  }
  
  resp
  
}


jns_WA_jnt2AJD <- function(x) {
  resp <- rep(NA, length(x))
  for(i in 1: length(x)) {
    resp[i] <- (jns_WA1AJD(x[i]) + jns_WA2AJD(x[i])) / 2
  }
  
  resp
  
}




plot(blh_WA1JD$JULIAN, blh_WA1JD$props, xlim = c(70, 310))
points(blh_WA2JD$JULIAN, blh_WA2JD$props)
lines(seq(1, 310), jns_WA_jntJD(seq(1, 310)), col = "blue", lwd = 2)

plot(blh_WA1JD$JULIAN, blh_WA1JD$CumProps, xlim = c(70, 310))
points(blh_WA2$JULIAN, blh_WA2$CumProps)
lines(seq(1, 310), jns_WA_jntAJD(seq(1, 310)), col = "blue", lwd = 2)




plot(blh_siteWA$JULIAN, blh_siteWA$propsBLH, xlim = c(70, 310))
lines(seq(1, 320), jns_WA_jnt2JD(seq(1, 320)), col = "blue", lwd = 2)

plot(blh_siteWA$JULIAN, blh_siteWA$CumProps, xlim = c(70, 310))
lines(seq(1, 320), jns_WA_jnt2AJD(seq(1, 320)), col = "blue", lwd = 2)









