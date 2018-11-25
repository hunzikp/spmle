########################################################################
# APPLICATION FOR UNIVARIATE SPATIOTEMPORAL MODEL
# CIVIL WAR IN (SUB-SAHARAN) AFRICA
########################################################################

library(cshapes)
library(rgeos)
library(lubridate)

##################################################
# Generate country sample (back-date from 2015)
##################################################

## Get baseline sample
countries.spdf <- cshp(date = as.Date("2015-06-01"), useGW = TRUE)
countries.spdf <- countries.spdf[(countries.spdf$GWCODE >= 400 & countries.spdf$GWCODE <= 626) | countries.spdf$GWCODE == 651, ]
countries.spdf@data <- countries.spdf@data[,c("GWCODE", "CNTRY_NAME"),drop=FALSE]
names(countries.spdf) <- c("gwid", "cname")
countries.spdf$cname <- as.character(countries.spdf$cname)

## Mark SSA
countries.spdf$ssa <- !(countries.spdf$gwid %in% c(615, 402, 651, 600, 616, 620))

## Mark islands
intrs.mat <- gIntersects(countries.spdf, byid = TRUE)
diag(intrs.mat) <- FALSE
hasneighbor <- apply(intrs.mat, 1, any)
countries.spdf$island <- !hasneighbor

## Order
countries.spdf <- countries.spdf[order(countries.spdf$gwid),]


##################################################
# Compute spatial weights matrix
##################################################

B <- gIntersects(countries.spdf, byid = TRUE)
diag(B) <- 0
denominator <- ifelse(countries.spdf$island, 1, rowSums(B))
W <- B / denominator


##################################################
# Generate monthly and yearly panels, 1996-2017
##################################################

startdates.month <- seq(as.Date("1996-01-01"), as.Date("2018-09-01"), by = "month")
enddates.month <- seq(as.Date("1996-02-01"), as.Date("2018-10-01"), by = "month")-1
mpanel.df <- do.call("rbind", rep(list(countries.spdf@data), length(startdates.month)))
mpanel.df$startdate <- rep(startdates.month, each = nrow(countries.spdf))
mpanel.df$enddate <- rep(enddates.month, each = nrow(countries.spdf))
mpanel.df$year <- lubridate::year(mpanel.df$startdate)

startdates.year <- seq(as.Date("1996-01-01"), as.Date("2017-01-01"), by = "year")
enddates.year <- seq(as.Date("1997-01-01"), as.Date("2018-01-01"), by = "year")-1
ypanel.df <- do.call("rbind", rep(list(countries.spdf@data), length(startdates.year)))
ypanel.df$startdate <- rep(startdates.year, each = nrow(countries.spdf))
ypanel.df$enddate <- rep(enddates.year, each = nrow(countries.spdf))
ypanel.df$year <- lubridate::year(ypanel.df$startdate)

##################################################
# Add controls
#   - GDP per captia
#   - Population
#   - Democracy (Polity IV)
#   - Neighborhoods
##################################################

## Load, merge gdp / pop
gdppop.df <- read.csv("~/Projects/gdppop/gdp18.csv")
gdppop.df$lgdp <- log(gdppop.df$rgdppc05)
gdppop.df$lpop <- log(gdppop.df$pop)
mpanel.df <- merge(mpanel.df, gdppop.df, by = c("gwid", "year"), all.x = TRUE, all.y = FALSE, sort = FALSE)
mpanel.df <- mpanel.df[order(mpanel.df$gwid, mpanel.df$startdate),]
ypanel.df <- merge(ypanel.df, gdppop.df, by = c("gwid", "year"), all.x = TRUE, all.y = FALSE, sort = FALSE)
ypanel.df <- ypanel.df[order(ypanel.df$gwid, ypanel.df$startdate),]

## Load, merge polity
polity.df <- read.csv("~/Data/polity/xpolity17.csv")
polity.df <- polity.df[,c("year", "gwid", "polity2")]
polity.df$polity2_sqrd <- polity.df$polity2^2
mpanel.df <- merge(mpanel.df, polity.df, by = c("gwid", "year"), all.x = TRUE, all.y = FALSE, sort = FALSE)
mpanel.df <- mpanel.df[order(mpanel.df$gwid, mpanel.df$startdate),]
ypanel.df <- merge(ypanel.df, polity.df, by = c("gwid", "year"), all.x = TRUE, all.y = FALSE, sort = FALSE)
ypanel.df <- ypanel.df[order(ypanel.df$gwid, ypanel.df$startdate),]

## Extrapolate gdp / pop / polity
variables <- c("lpop", "lgdp", "polity2", "polity2_sqrd")
countries <- unique(mpanel.df$gwid)

# Function doing simple inter/extrapolation
#   If no past obs available -> Backdate time-series
#   If past obs available -> Extrapolate time-series
simplepolate <- function(x) {
  for (i in 1:length(x)) {
    if (is.na(x[i])) {
      last.observed <- NA
      if (i > 1) {
        x.prev <- x[1:(i-1)]
        if (!all(is.na(x.prev))) {
          last.observed <- x.prev[max(which(!is.na(x.prev)))]
        }
      }
      next.observed <- NA
      if (i < length(x)) {
        x.fut <- x[(i+1):length(x)]
        if (!all(is.na(x.fut))) {
          next.observed <- x.fut[min(which(!is.na(x.fut)))]
        }
      }
      if (!is.na(last.observed)) {
        x[i] <- last.observed
      } else {
        x[i] <- next.observed
      }
    }
  }
  return(x)
}

# Monthly
for (gw in countries) {
  for (vn in variables) {
    var <- mpanel.df[mpanel.df$gwid == gw, vn]
    mpanel.df[mpanel.df$gwid == gw, vn] <- simplepolate(var)
  }
}

# Yearly
for (gw in countries) {
  for (vn in variables) {
    var <- ypanel.df[ypanel.df$gwid == gw, vn]
    ypanel.df[ypanel.df$gwid == gw, vn] <- simplepolate(var)
  }
}

## Neighborhood vars
variables <- c("lpop", "lgdp", "polity2", "polity2_sqrd")

# Monthly
startdates <- mpanel.df$startdate
for (vn in variables) {
  var <- mpanel.df[,vn]
  spvar <- rep(NA, length(var))
  for (date in unique(startdates)) {
    var.t <- var[startdates == date]
    # Construct spatial-lag so that it includes non-missings only
    observed <- !is.na(var.t)
    B.t <- B[observed, observed]
    denominator <- ifelse(rowSums(B.t) == 0, 1, rowSums(B.t))
    W.t <- B.t/denominator
    spvar.t <- rep(NA, length(var.t))
    spvar.t[observed] <- var.t[observed]%*%W.t
    spvar[startdates == date] <- spvar.t
  }
  spvn <- paste0(vn, "_splag")
  mpanel.df$newvar <- spvar
  names(mpanel.df)[ncol(mpanel.df)] <- spvn
}

# Yearly
startdates <- ypanel.df$startdate
for (vn in variables) {
  var <- ypanel.df[,vn]
  spvar <- rep(NA, length(var))
  for (date in unique(startdates)) {
    var.t <- var[startdates == date]
    # Construct spatial-lag so that it includes non-missings only
    observed <- !is.na(var.t)
    B.t <- B[observed, observed]
    denominator <- ifelse(rowSums(B.t) == 0, 1, rowSums(B.t))
    W.t <- B.t/denominator
    spvar.t <- rep(NA, length(var.t))
    spvar.t[observed] <- var.t[observed]%*%W.t
    spvar[startdates == date] <- spvar.t
  }
  spvn <- paste0(vn, "_splag")
  ypanel.df$newvar <- spvar
  names(ypanel.df)[ncol(ypanel.df)] <- spvn
}

##################################################
# Add ACLED outcome variables
##################################################

## Load ACLED
acled18.df <- read.csv("/home/hunzikp/Data/acled/2018/Africa_1997-2018_Sep29.csv", stringsAsFactors = FALSE, fileEncoding="latin1")
keep <- c("ISO", "EVENT_ID_CNTY", "EVENT_DATE", "YEAR", "TIME_PRECISION",
          "EVENT_TYPE", "ACTOR1", "ACTOR2",
          "INTERACTION",
          "COUNTRY", "LOCATION", "LATITUDE",
          "LONGITUDE","FATALITIES", "NOTES")
acled18.df <- acled18.df[,keep]
acled.df <- acled18.df

## Correct dates
acled.df$EVENT_DATE <- as.Date(acled.df$EVENT_DATE, "%d-%B-%Y")
names(acled.df) <- tolower(names(acled.df))

## Add event type markers
acled.df$battle_event <- grepl(acled.df$event_type, pattern = 'attle')
acled.df$protest_event <- grepl(acled.df$event_type, pattern = 'rotest')
acled.df$civilian_event <- grepl(acled.df$event_type, pattern = 'ivilian') | grepl(acled.df$actor2, pattern = 'ivilian')

## Get rid of notes
acled.df <- acled.df[,-which(names(acled.df) == "notes")]

## Add new country names
acled.df$cname <- acled.df$country
countries.spdf$cname[!(countries.spdf$cname %in% unique(acled.df$cname))]
acled.df$cname[acled.df$cname == "Democratic Republic of Congo"] <- "Congo, DRC"
acled.df$cname[acled.df$cname == "Ivory Coast"] <- "Cote d'Ivoire"
acled.df$cname[acled.df$cname == "eSwatini"] <- "Swaziland"
acled.df$cname[acled.df$cname == "Republic of Congo"] <- "Congo"
acled.df$cname[acled.df$cname == "Gambia"] <- "The Gambia"
countries.spdf$cname[!(countries.spdf$cname %in% unique(acled.df$cname))]

## Map ACLED onto panel, monthly
variables <- c("battle_event", "protest_event", "civilian_event")
startdates <- mpanel.df$startdate
enddates <- mpanel.df$enddate
for (vn in variables) {
  newvar <- rep(NA, length(startdates))
  for (date in unique(startdates)) {
    start <- date
    end <- unique(enddates[startdates == start])
    eventvar <- acled.df[,vn]
    event.df <- acled.df[acled.df$event_date >= start & acled.df$event_date <= end & eventvar, "cname", drop = FALSE]
    if (nrow(event.df) == 0) {next}
    event.df$event <- 1
    event.df <- aggregate(event.df$event, by = list(event.df$cname), FUN = sum)
    names(event.df) <- c("cname", vn)
    cntr.df <- countries.spdf@data[,c("cname", "gwid")]
    event.df <- merge(cntr.df, event.df, by = "cname", all.x = TRUE, all.y = FALSE)
    event.df <- event.df[order(event.df$gwid),]
    event.df[is.na(event.df[,vn]),vn] <- 0
    newvar[startdates == start] <- event.df[,vn]
  }
  mpanel.df$newvar <- newvar
  names(mpanel.df)[ncol(mpanel.df)] <- vn
}

## Map ACLED onto panel, yearly
variables <- c("battle_event", "protest_event", "civilian_event")
startdates <- ypanel.df$startdate
enddates <- ypanel.df$enddate
for (vn in variables) {
  newvar <- rep(NA, length(startdates))
  for (date in unique(startdates)) {
    start <- date
    end <- unique(enddates[startdates == start])
    eventvar <- acled.df[,vn]
    event.df <- acled.df[acled.df$event_date >= start & acled.df$event_date <= end & eventvar, "cname", drop = FALSE]
    if (nrow(event.df) == 0) {next}
    event.df$event <- 1
    event.df <- aggregate(event.df$event, by = list(event.df$cname), FUN = sum)
    names(event.df) <- c("cname", vn)
    cntr.df <- countries.spdf@data[,c("cname", "gwid")]
    event.df <- merge(cntr.df, event.df, by = "cname", all.x = TRUE, all.y = FALSE)
    event.df <- event.df[order(event.df$gwid),]
    event.df[is.na(event.df[,vn]),vn] <- 0
    newvar[startdates == start] <- event.df[,vn]
  }
  ypanel.df$newvar <- newvar
  names(ypanel.df)[ncol(ypanel.df)] <- vn
}

##################################################
# Add UCDP outcome variable
##################################################

## Load
acd.df <- read.csv("~/Data/ucdp/acd/ucdp-prio-acd-181.csv", stringsAsFactors = FALSE)

## Subset to relevant conflicts / countries
acd.df <- acd.df[acd.df$type_of_conflict %in% c(3,4),]
acd.df$gwid <- as.numeric(acd.df$gwno_loc)
acd.df <- acd.df[acd.df$gwid %in% countries.spdf$gwid,]

## Get episodes table
last_date <- ""
acdep.df <- unique(acd.df[, c("conflict_id", "gwid", "start_date2", "ep_end_date")])
acdep.df$keep <- FALSE
for (i in 1:nrow(acdep.df)) {
  epd <- acdep.df$ep_end_date[i]
  if (epd != "") {
    acdep.df$keep[i] <- TRUE
  } else {
    # Identify ongoing episodes...
    # Check whether any other row with same conflictid-startdate has ep-end-date
    # If not, keep and assign last possible date as end date
    cid <- acdep.df$conflict_id[i]
    sd <- acdep.df$start_date2[i]
    equal.df <- acdep.df[acdep.df$conflict_id == cid & acdep.df$start_date2 == sd,]
    epds <- equal.df$ep_end_date
    if (all(epds == "")) {
      acdep.df$keep[i] <- TRUE
      acdep.df$ep_end_date[i] <- "2018-09-30"
    }
  }
}
acdep.df <- acdep.df[acdep.df$keep,]

## Aggregate to country year
yacd.df <- ypanel.df[,c("gwid", "year", "startdate", "enddate")]
yacd.df$incidence <- 0
for (i in 1:nrow(acdep.df)) {
  gwid <- acdep.df$gwid[i]
  sd <- acdep.df$start_date2[i]
  ed <- acdep.df$ep_end_date[i]
  yacd.df$incidence[yacd.df$gwid == gwid & yacd.df$startdate <= ed & yacd.df$enddate >= sd] <- 1
}
ypanel.df$incidence <- yacd.df$incidence

## Aggregate to country month
macd.df <- mpanel.df[,c("gwid", "year", "startdate", "enddate")]
macd.df$incidence <- 0
for (i in 1:nrow(acdep.df)) {
  gwid <- acdep.df$gwid[i]
  sd <- acdep.df$start_date2[i]
  ed <- acdep.df$ep_end_date[i]
  macd.df$incidence[macd.df$gwid == gwid & macd.df$startdate <= ed & macd.df$enddate >= sd] <- 1
}
mpanel.df$incidence <- macd.df$incidence

##################################################
# Temporal lags for predictors
##################################################

vars <- c("lpop", "lgdp", "polity2", "polity2_sqrd", "lgdp_splag", "polity2_splag", "polity2_sqrd_splag")

## Yearly panel
ypanel_lag.df <- ypanel.df[,c("gwid", "year", vars)]
names(ypanel_lag.df)[-c(1,2)] <- paste0(vars, "_lag")
ypanel_lag.df$year <- ypanel_lag.df$year + 1
ypanel.df <- merge(ypanel.df, ypanel_lag.df, by = c("gwid", "year"), sort = FALSE, all.x = TRUE, all.y = FALSE)
ypanel.df <- ypanel.df[order(ypanel.df$gwid, ypanel.df$startdate),]

## Monthly panel
mpanel_lag.df <- unique(mpanel.df[,c("gwid", "year", vars)])
names(mpanel_lag.df)[-c(1,2)] <- paste0(vars, "_lag")
mpanel_lag.df$year <- mpanel_lag.df$year + 1
mpanel.df <- merge(mpanel.df, mpanel_lag.df, by = c("gwid", "year"), sort = FALSE, all.x = TRUE, all.y = FALSE)
mpanel.df <- mpanel.df[order(mpanel.df$gwid, mpanel.df$startdate),]


##################################################
# Write out
##################################################

saveRDS(ypanel.df, file = "extdata/africa_ypanel.rds")
saveRDS(mpanel.df, file = "extdata/africa_mpanel.rds")
saveRDS(countries.spdf, file = "extdata/africa_countries.rds")
saveRDS(W, file = "extdata/africa_W.rds")

