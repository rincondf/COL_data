# DAS function to calculate degree days (single sin method)

calc_dd <- function(tmax, tmin, lower_threshold, upper_threshold, cutoff){
  
  # No degree day accumulation:
  if(tmin > tmax | tmax <= lower_threshold) {
    return(0)
  }
  
  sum_heat <- tmax + tmin
  diff_heat <- tmax - tmin
  
  # When vertical cutoff for upper threshold is used...
  if(cutoff == "vertical") {
    
    # No degree day accumulation: Options 1 & 2 from T5 in Zalom et al
    if(tmin >= upper_threshold | tmax <= lower_threshold) {
      return(0)
    }
    
    # Full temp range is within thresholds: Option 3 in T5 of Zalom et al
    if(lower_threshold <= tmin & upper_threshold >= tmax) {
      return((sum_heat / 2) - lower_threshold)
    }
    
    # Looks like "a" corresponds roughly to a portion of theta-1 in the Zalom diagrams
    a <- 2 * lower_threshold - sum_heat
    
    if(abs(diff_heat) > abs(a)) {
      b <- atan(a / sqrt(diff_heat * diff_heat - a * a))
    } else {
      b <- 0
    }
    
    # Looks like "c" corresponds roughly to a portion of theta-2 in the Zalom diagrams
    c <- 2 * upper_threshold - sum_heat
    
    if(abs(diff_heat) > abs(c)) {
      d <- atan(c / sqrt(diff_heat * diff_heat - c * c))
    } else {
      d <- 0
    }
    
    # Upper threshold limiting heat accumulation: Option 5 in T5 of Zalom et al
    if(lower_threshold <= tmin) {
      return((-diff_heat * cos(d) - a * (d + (0.5 * pi))) / (2 * pi))
    }
    
    # Both thresholds limiting heat accumulation: Option 6 in T5 of Zalom et al
    if(upper_threshold < tmax) {
      return((-diff_heat * (cos(d) - cos(b)) - a * (d - b)) / (2 * pi))
    }
    
    # Default to lower temp threshold limiting heat accumulation: Option 4
    # of T5 in Zalom et al
    return((diff_heat * cos(b) - a * ((0.5 * pi) - b)) / (2 * pi))
    
    
    # When a non-vertical cutoff method is used...
  } else {
    
    a <- sum_heat / 2 - lower_threshold
    
    # Full temp range is above both thresholds: Options 1 in T5 of Zalom et al
    if (tmin >= upper_threshold && tmax > upper_threshold) {
      return(upper_threshold - lower_threshold)}
    
    # Full temp range is within thresholds: Option 3 in T5 of Zalom et al
    if(lower_threshold <= tmin & upper_threshold >= tmax) {
      return(a)
    }
    
    # Part of theta-1
    b <- 2 * lower_threshold - sum_heat
    
    
    if(abs(diff_heat) > abs(b)) {
      c <- atan(b / sqrt(diff_heat * diff_heat - b * b))
      d <- (diff_heat * cos(c) - b * ((0.5 * pi) - c)) / (2 * pi)
      # Below both thresholds: Option 2 of T5 in Zalom et al
    } else {
      d <- 0
    }
    
    # Lower temp threshold limiting heat accumulation: Option 4 of T5 in Zalom
    # et al (or Option 2)
    if(upper_threshold >= tmax) {
      return(d)
    }
    
    # Part of theta-2
    e <- 2 * upper_threshold - sum_heat
    
    f <- atan(e / sqrt(diff_heat * diff_heat - e * e))
    
    g <- (diff_heat * cos(f) - e * ((0.5 * pi) - f)) / (2 * pi)
    
    
    # Both thresholds limiting heat accumulation: Option 6 in T5 of Zalom et al
    if(lower_threshold > tmin) {
      return((d - g))
    }
    
    # Default to upper threshold limiting heat accumulation: Option 5 in T5 of
    # Zalom et al
    return((a - g))
  }
  
}


# Allows sets of temperatures

calc_dd_vec <- function(tmax, tmin, lower_threshold, upper_threshold, cutoff) {
  ddss <- rep(NA, length(tmax))
  for(i in 1: length(tmax)) {
    ddss[i] <- calc_dd(tmax[i], tmin[i], lower_threshold, upper_threshold, cutoff)
  }
  ddss
}

