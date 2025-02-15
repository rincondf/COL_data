# DEGREE DAY CALCULATION AND PHENOLOGY MODELING OF BEET LEAF HOPPER IN COLORADO

# PRERAPRED BY DIEGO RINCON (diego.rincon@wsu.edu)

# 1. Formalize dates and calculate julians 

COL_data$date <- as.Date(COL_data$date, format = "%m/%d/%Y")
COL_data$Julian <- as.numeric(format(COL_data$date, "%j"))

# 2. Find unique locations and years combinations

locations_COL <- subset(COL_data, 
                        duplicated(COL_data[c("latitude", 
                                              "longitude", 
                                              "Year")]) == 
                          FALSE)[, c(1, 2, 3, 4, 5)]

write.csv(locations_COL, file = "locationsCOL.csv")

# 3. Extract max and min temp per location and year

library(daymetr)

temp_recs <- list()

for(i in 1: length(locations_COL$site)) {
  temp_recs[[i]] <- download_daymet(
    site = "Daymet",
    lat = locations_COL$latitude[i],
    lon = locations_COL$longitude[i],
    start = locations_COL$Year[i],
    end = locations_COL$Year[i],
    path = tempdir(),
    internal = TRUE,
    silent = FALSE,
    force = FALSE,
    simplify = FALSE
  )
}


tmax <- list()

for(i in 1: length(locations_COL$site)) {
  tmax[[i]] <- as.numeric(temp_recs[[i]]$data[,7])
}

tmin <- list()

for(i in 1: length(locations_COL$site)) {
  tmin[[i]] <- as.numeric(temp_recs[[i]]$data[,8])
}


# 4. Calculate cumulative degree days for both BLH and GPaphid using retrieved
# max and min temps and the DAS function (single sin) for DDs calculation

# BLH

upT_BLH <- 10 + 19.8
baseT_BLH <- 10

DDs_BLH <- list()

for(i in 1: length(locations_COL$site)) {
  DDs_BLH[[i]] <- calc_dd_vec(tmax = tmax[[i]], tmin = tmin[[i]], 
                              lower_threshold = baseT_BLH, 
                              upper_threshold = upT_BLH, 
                              cutoff = "horizontal")
  DDs_BLH[[i]] <- cumsum(DDs_BLH[[i]])
}

# Matching data collection days

COL_data$BLH_DDs <- rep(NA, length(COL_data$date))


for(i in 1: length(locations_COL$site)) {
  COL_data$BLH_DDs[which(COL_data$latitude == locations_COL$latitude[i] & 
                           COL_data$longitude == locations_COL$longitude[i])] <-
    DDs_BLH[[i]][COL_data$Julian[which(COL_data$latitude == 
                                         locations_COL$latitude[i] &
                                         COL_data$longitude == 
                                         locations_COL$longitude[i])]]
}


# aggregate mean captures by julian days

blh_julians <- aggregate(COL_data$BLH, list(COL_data$date), mean, na.rm = TRUE)

names(blh_julians) <- c("date", "blh_mean")
blh_julians$Julian <- as.numeric(format(blh_julians$date, "%j"))

plot(blh_julians$Julian, blh_julians$blh_mean, xlab = "Julian days", 
     ylab = "Mean captures")

# aggregate mean captures by degree days

blh_dds <- aggregate(COL_data$BLH, list(COL_data$Year, COL_data$BLH_DDs), mean, 
                     na.rm = TRUE)

names(blh_dds) <- c("year", "dds", "blh_mean")

blh_dds <- blh_dds[-which(blh_dds$blh_mean == "NaN"), ]


plot(blh_dds$dds, blh_dds$blh_mean, xlab = "Degree days", 
     ylab = "Mean captures")

# aggregate mean captures by site

blh_site <- aggregate(COL_data[, -c(1, 2)], list(COL_data$site, COL_data$date),
                      mean, na.rm = TRUE)

names(blh_site)[1: 2] <- c("site", "date")
blh_site <- blh_site[-which(blh_site$BLH == "NaN"), ]

# 5. Calculation of proportions and cumulative proportions captured per site 
# per collection date

years <- seq(2015, 2024)

# By julian days

# Proportions

blh_julians$props <- rep(NA, length(blh_julians$blh_mean))

for(i in 1: length(years)) {
  blh_julians$props[which(blh_julians$date > 
                            paste0(years[i], "-01-01") & 
                            blh_julians$date < 
                            paste0(years[i], "-12-31"))] <- 
    blh_julians$blh_mean[which(blh_julians$date > 
                                 paste0(years[i], "-01-01") & 
                                 blh_julians$date < 
                                 paste0(years[i], "-12-31"))] / 
    sum(blh_julians$blh_mean[which(blh_julians$date > 
                                     paste0(years[i], "-01-01") & 
                                     blh_julians$date < 
                                     paste0(years[i], "-12-31"))])
}


plot(blh_julians$Julian, blh_julians$props, xlab = "Julian days", 
     ylab = "Proportion of captures")


# Cumulative proportions


blh_julians$Cum_props <- rep(NA, length(blh_julians$blh_mean))

for(i in 1: length(years)) {
  blh_julians$Cum_props[which(blh_julians$date > 
                                paste0(years[i], "-01-01") & 
                                blh_julians$date < 
                                paste0(years[i], "-12-31"))] <- 
    cumsum(blh_julians$blh_mean[which(blh_julians$date > 
                                        paste0(years[i], "-01-01") & 
                                        blh_julians$date < 
                                        paste0(years[i], "-12-31"))]) / 
    sum(blh_julians$blh_mean[which(blh_julians$date > 
                                     paste0(years[i], "-01-01") & 
                                     blh_julians$date < 
                                     paste0(years[i], "-12-31"))])
}


plot(blh_julians$Julian, blh_julians$Cum_props, xlab = "Julian days", 
     ylab = "Proportion of captures")


# By degree days

# Proportions

blh_dds$props <- rep(NA, length(blh_dds$blh_mean))

for(i in 1: length(years)) {
  blh_dds$props[which(blh_dds$year == years[i])] <- 
    blh_dds$blh_mean[which(blh_dds$year == years[i])] / 
    sum(blh_dds$blh_mean[which(blh_dds$year == years[i])],
        na.rm = TRUE)
}


plot(blh_dds$dds, blh_dds$props, xlab = "Degree days", 
     ylab = "Proportion of captures")


# Cumulative proportions


blh_dds$Cum_props <- rep(NA, length(blh_dds$blh_mean))

for(i in 1: length(years)) {
  blh_dds$Cum_props[which(blh_dds$year == years[i])] <- 
    cumsum(blh_dds$blh_mean[which(blh_dds$year == years[i])]) / 
    sum(blh_dds$blh_mean[which(blh_dds$year == years[i])],
        na.rm = TRUE)
}


plot(blh_dds$dds, blh_dds$Cum_props, xlab = "Degree days", 
     ylab = "Proportion of captures")


# By site

# Proportions

blh_site$propsBLH <- rep(NA, length(blh_site$BLH))

for(i in 1: length(locations_COL$site)) {
  blh_site$propsBLH[which(blh_site$site == locations_COL$site[i] & 
                            blh_site$Year == locations_COL$Year[i])] <-
    blh_site$BLH[which(blh_site$site == locations_COL$site[i] & 
                         blh_site$Year == locations_COL$Year[i])] /
    sum(blh_site$BLH[which(blh_site$site == locations_COL$site[i] & 
                             blh_site$Year == locations_COL$Year[i])])
}


plot(blh_site$BLH_DDs, blh_site$propsBLH)
plot(blh_site$Julian, blh_site$propsBLH)

# Cumulative proportions

# Checking if retrieved locations and years are in ascending order

# for julian days

a <- matrix(NA, length(locations_COL$site), 15)

for(i in 1: length(locations_COL$site)) {
  a[i, 1: length(which(blh_site$site == locations_COL$site[i] & 
                         blh_site$Year == locations_COL$Year[i]))] <- 
    blh_site$Julian[which(blh_site$site == locations_COL$site[i] & 
                            blh_site$Year == locations_COL$Year[i])]
}

check1 <- rep(NA, length(locations_COL$site))

for(i in 1: length(locations_COL$site)) {
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
