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
                        list(CLEAN_BLH$`FIELD ID`, 
                             CLEAN_BLH$SURVEY_DATE), mean, na.rm = TRUE)

names(blh_siteWA)[1: 2] <- c("site", "date")
blh_siteWA <- blh_siteWA[-which(blh_siteWA$TRAP_COUNT_RAW == "NaN"), ]

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


# By site

# Proportions

blh_siteWA$propsBLH <- rep(NA, length(blh_siteWA$TRAP_COUNT_RAW))

for(i in 1: length(locations_WA$`FIELD ID`)) {
  blh_siteWA$propsBLH[which(blh_siteWA$site == 
                              locations_WA$`FIELD ID`[i] & 
                              blh_siteWA$year == 
                              locations_WA$year[i])] <-
    blh_siteWA$TRAP_COUNT_RAW[which(blh_siteWA$site == 
                                      locations_WA$`FIELD ID`[i] & 
                                      blh_siteWA$year == 
                                      locations_WA$year[i])] /
    sum(blh_siteWA$TRAP_COUNT_RAW[which(blh_siteWA$site == 
                                          locations_WA$`FIELD ID`[i] & 
                                          blh_siteWA$year == 
                                          locations_WA$year[i])])
}


plot(blh_siteWA$DDs, blh_siteWA$propsBLH)
plot(blh_siteWA$JULIAN, blh_siteWA$propsBLH)

# Cumulative proportions

# Checking if retrieved locations and years are in ascending order

# for julian days

a <- matrix(NA, length(locations_WA$`FIELD ID`), 30)

for(i in 1: length(locations_WA$`FIELD ID`)) {
  a[i, 1: length(which(blh_siteWA$site == locations_WA$`FIELD ID`[i] & 
                         blh_siteWA$year == locations_WA$year[i]))] <- 
    blh_siteWA$JULIAN[which(blh_siteWA$site == locations_WA$`FIELD ID`[i] & 
                              blh_siteWA$year == locations_WA$year[i])]
}

check1 <- rep(NA, length(locations_WA$`FIELD ID`))

for(i in 1: length(locations_WA$`FIELD ID`)) {
  check1[i] <- any(diff(a[i, ]) < 0, na.rm = TRUE)
}

check1 # they are

# for degree days

b <- matrix(NA, length(locations_COL$site), 15)

for(i in 1: length(locations_COL$site)) {
  b[i, 1: length(which(blh_site$site == locations_COL$site[i] & 
                         blh_site$Year == locations_COL$Year[i]))] <- 
    blh_site$BLH_DDs[which(blh_site$site == locations_COL$site[i] & 
                             blh_site$Year == locations_COL$Year[i])]
}

check2 <- rep(NA, length(locations_COL$site))

for(i in 1: length(locations_COL$site)) {
  check2[i] <- any(diff(b[i, ]) < 0, na.rm = TRUE)
}

check2 # they are



blh_site$CumPropsBLH <- rep(NA, length(blh_site$BLH))

for(i in 1: length(locations_COL$site)) {
  blh_site$CumPropsBLH[which(blh_site$site == locations_COL$site[i] & 
                               blh_site$Year == locations_COL$Year[i])] <-
    cumsum(blh_site$BLH[which(blh_site$site == locations_COL$site[i] & 
                                blh_site$Year == locations_COL$Year[i])]) /
    sum(blh_site$BLH[which(blh_site$site == locations_COL$site[i] & 
                             blh_site$Year == locations_COL$Year[i])])
}


plot(blh_site$BLH_DDs, blh_site$CumPropsBLH)
plot(blh_site$Julian, blh_site$CumPropsBLH)


totals_site <- rep(NA, length(locations_COL$site))


for(i in 1: length(locations_COL$site)) {
  totals_site[i] <- sum(blh_site$BLH[which(blh_site$site == 
                                             locations_COL$site[i] & 
                                             blh_site$Year == 
                                             locations_COL$Year[i])])
  
  
}


locations_COL2 <- locations_COL[which(totals_site > 30), ]


blh_site$CumPropsBLH2 <- rep(NA, length(blh_site$BLH))

for(i in 1: length(locations_COL2$site)) {
  blh_site$CumPropsBLH2[which(blh_site$site == locations_COL2$site[i] & 
                                blh_site$Year == locations_COL2$Year[i])] <-
    cumsum(blh_site$BLH[which(blh_site$site == locations_COL2$site[i] & 
                                blh_site$Year == locations_COL2$Year[i])]) /
    sum(blh_site$BLH[which(blh_site$site == locations_COL2$site[i] & 
                             blh_site$Year == locations_COL2$Year[i])])
}


plot(blh_site$BLH_DDs, blh_site$CumPropsBLH2)
plot(blh_site$Julian, blh_site$CumPropsBLH2)


blh_site$propsBLH2 <- rep(NA, length(blh_site$BLH))

for(i in 1: length(locations_COL2$site)) {
  blh_site$propsBLH2[which(blh_site$site == locations_COL2$site[i] & 
                             blh_site$Year == locations_COL2$Year[i])] <-
    blh_site$BLH[which(blh_site$site == locations_COL2$site[i] & 
                         blh_site$Year == locations_COL2$Year[i])] /
    sum(blh_site$BLH[which(blh_site$site == locations_COL2$site[i] & 
                             blh_site$Year == locations_COL2$Year[i])])
}


plot(blh_site$BLH_DDs, blh_site$propsBLH2)
plot(blh_site$Julian, blh_site$propsBLH2)
