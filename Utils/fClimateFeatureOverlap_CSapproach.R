# CLIMATE-PRIORITY-AREA APPROACH (Buenafe et al. 2023)
# Modified to calculate the percentile cut-off for the distribution of each feature 

# ----- Description -----
# This function splits the feature data into climate-smart (CS) and non-climate-smart (NCS) areas depending on the percentile chosen by the user.
# Inputs:
# featuresDF: feature sf object which should have a column for cellID
# percentile: cut-off threshold for determining whether an area is a climate priority area or not. Note that the percentile here is the lower limit of the threshold.
# metricDF: climate metric data.frame with 'metric' as the column name of the metric values per planning unit. This should also have a column for the cellID
# direction: If direction = 1, metric values are from low (least climate-smart) to high (most climate-smart). If direction = -1, metric values are from high (least climate-smart) to low (most climate-smart).
# targetsDF: data.frame that has the list of features under "feature" column and their corresponding targets under "target" column
# refugiaTarget: target assigned to climate-smart areas



# Input: One metricDFs
# Input: Single metricDF
fClimateFeatureOverlap_CSapproach_Combined2_bin <- function(featuresDF,
                                                            percentile,
                                                            metricDF,
                                                            direction,
                                                            upper_percentile) {
  
  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  prct_half <- (percentile / 2) / 100
  
  if (any(is.na(metricDF$metric))) {
    stop("There are some NAs in the metric data. Please check.")
  }
  
  threshold_upper <- metricDF %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    summarise(qntl = quantile(metric, 1 - prct_half)[[1]]) %>%
    pull(qntl)
  
  threshold_lower <- metricDF %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    summarise(qntl = quantile(metric, prct_half)[[1]]) %>%
    pull(qntl)
  
  if (direction == 1) {
    message("Higher values refer to more climate-resilient areas.")
    metricCS <- metricDF %>%
      dplyr::mutate(V0 = ifelse(metric >= threshold_upper | metric <= threshold_lower, 1, 0))
  } else if (direction == -1) {
    message("Lower values refer to more climate-resilient areas.")
    metricCS <- metricDF %>%
      dplyr::mutate(V0 = ifelse(metric <= threshold_lower | metric >= threshold_upper, 1, 0))
  } else {
    stop("Please enter a valid direction: either 1 or -1.")
  }
  
  imptList <- list() # empty list to fill in with the important features
  for (i in seq_along(spp)) {
    df <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID")
    
    imptList[[i]] <- df %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS")) := V0 * .data[[spp[i]]]) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
  }
  
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    select(cellID, everything())
  
  repList <- list()
  for (i in seq_along(spp)) {
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID")
    
    df2 <- imptList %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
    
    repList[[i]] <- df1 %>%
      dplyr::left_join(df2, by = "cellID") %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_NCS")) := ifelse(!!sym(paste0(spp[i], "_CS")) == 1, 0, .data[[spp[i]]])) %>%
      dplyr::select(matches('_NCS|_CS'))
  }
  
  repList <- repList %>%
    do.call(bind_cols, .) %>%
    dplyr::bind_cols(featuresDF %>% dplyr::select(cellID, geometry)) %>%
    dplyr::select(cellID, everything()) %>%
    sf::st_as_sf(sf_column_name = "geometry")
  
  return(repList)
}





# Input: One metricDFs
fClimateFeatureOverlap_CSapproach_Combined2 <- function(featuresDF,
                                                       percentile,
                                                       metricDF,
                                                       direction,
                                                       upper_percentile, cut) {
  
  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  prct_half <- (percentile / 2) / 100
  
  if (any(is.na(metricDF$metric))) {
    stop("There are some NAs in the metric data. Please check.")
  }
  
  if(cut == "YES") {
    metricDF <- metricDF %>%
      dplyr::inner_join(featuresDF %>% dplyr::select(cellID), by = "cellID")
  } else {
    metricDF <- metricDF
  }
  
  threshold_upper <- metricDF %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    summarise(qntl = quantile(metric, 1 - prct_half)[[1]]) %>%
    pull(qntl)
  
  threshold_lower <- metricDF %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    summarise(qntl = quantile(metric, prct_half)[[1]]) %>%
    pull(qntl)
  
  if (direction == 1) {
    message("Higher values refer to more climate-resilient areas.")
    metricCS <- metricDF %>%
      dplyr::mutate(V0 = ifelse(metric >= threshold_upper | metric <= threshold_lower, 1, 0))
  } else if (direction == -1) {
    message("Lower values refer to more climate-resilient areas.")
    metricCS <- metricDF %>%
      dplyr::mutate(V0 = ifelse(metric <= threshold_lower | metric >= threshold_upper, 1, 0))
  } else {
    stop("Please enter a valid direction: either 1 or -1.")
  }
  
  imptList <- list() # empty list to fill in with the important features
  for (i in seq_along(spp)) {
    df <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID")
    
    imptList[[i]] <- df %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS")) := V0 * .data[[spp[i]]]) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
  }
  
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    select(cellID, everything())
  
  repList <- list()
  for (i in seq_along(spp)) {
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID")
    
    df2 <- imptList %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
    
    repList[[i]] <- df1 %>%
      dplyr::left_join(df2, by = "cellID") %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_NCS")) := ifelse(!!sym(paste0(spp[i], "_CS")) > 0, 0, .data[[spp[i]]])) %>%
      dplyr::select(matches('_NCS|_CS'))
  }
  
  repList <- repList %>%
    do.call(bind_cols, .) %>%
    dplyr::bind_cols(featuresDF %>% dplyr::select(cellID, geometry)) %>%
    dplyr::select(cellID, everything()) %>%
    sf::st_as_sf(sf_column_name = "geometry")
  
  return(repList)
}




# Input: Two metricDFs
fClimateFeatureOverlap_CSapproach_Combined_bin <- function(featuresDF,
                                                           percentile,
                                                           metricDF_AH,
                                                           metricDF_R,
                                                           direction,
                                                           upper_percentile) {
  
  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  prct_half <- (percentile / 2) / 100
  
  if (any(is.na(metricDF_AH$metric)) || any(is.na(metricDF_R$metric))) {
    stop("There are some NAs in the metric data. Please check.")
  }
  
  threshold_AH <- metricDF_AH %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    summarise(qntl = if (upper_percentile) quantile(metricDF_AH$metric, 1 - prct_half)[[1]] else quantile(metricDF_AH$metric, prct_half)[[1]]) %>%
    pull(qntl)
  
  threshold_R <- metricDF_R %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    summarise(qntl = if (upper_percentile) quantile(metricDF_R$metric, 1 - prct_half)[[1]] else quantile(metricDF_R$metric, prct_half)[[1]]) %>%
    pull(qntl)
  
  if (direction == 1) {
    message("Higher values refer to more climate-resilient areas.")
    metricCS_AH <- metricDF_AH %>%
      dplyr::mutate(V0_AH = ifelse(metricDF_AH$metric >= threshold_AH, 1, 0))
    metricCS_R <- metricDF_R %>%
      dplyr::mutate(V0_R = ifelse(metricDF_R$metric >= threshold_R, 1, 0))
  } else if (direction == -1) {
    message("Lower values refer to more climate-resilient areas.")
    metricCS_AH <- metricDF_AH %>%
      dplyr::mutate(V0_AH = ifelse(metricDF_AH$metric <= threshold_AH, 1, 0))
    metricCS_R <- metricDF_R %>%
      dplyr::mutate(V0_R = ifelse(metricDF_R$metric <= threshold_R, 1, 0))
  } else {
    stop("Please enter a valid direction: either 1 or -1.")
  }
  
  metricCS <- metricCS_AH %>%
    sf::st_drop_geometry() %>%
    dplyr::select(cellID, V0_AH) %>%
    dplyr::left_join(metricCS_R %>%
                       sf::st_drop_geometry() %>%
                       dplyr::select(cellID, V0_R), by = "cellID") %>%
    dplyr::mutate(V0 = ifelse(V0_AH == 1 | V0_R == 1, 1, 0)) %>%
    dplyr::select(cellID, V0) %>%
    sf::st_as_sf(., geometry = metricCS_AH$geometry)
  
  imptList <- list() # empty list to fill in with the important features
  for (i in seq_along(spp)) {
    df <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID") 
    
    imptList[[i]] <- df %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS")) := V0 * .data[[spp[i]]]) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
  }
  
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    select(cellID, everything())
  
  repList <- list()
  for (i in seq_along(spp)) {
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID")
    
    df2 <- imptList %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
    
    repList[[i]] <- df1 %>%
      dplyr::left_join(df2, by = "cellID") %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_NCS")) := ifelse(!!sym(paste0(spp[i], "_CS")) == 1, 0, .data[[spp[i]]])) %>%
      dplyr::select(matches('_NCS|_CS'))
  }
  
  repList <- repList %>%
    do.call(bind_cols, .) %>%
    dplyr::bind_cols(featuresDF %>% dplyr::select(cellID, geometry)) %>%
    dplyr::select(cellID, everything()) %>%
    sf::st_as_sf(sf_column_name = "geometry")
  
  return(repList)
}




# Input: Two metricDFs
fClimateFeatureOverlap_CSapproach_Combined <- function(featuresDF,
                                              percentile,
                                              metricDF_AH,
                                              metricDF_R,
                                              direction,
                                              upper_percentile) {

  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  prct_half <- (percentile / 2) /100
  
  if (any(is.na(metricDF_AH$metric)) || any(is.na(metricDF_R$metric))) {
    stop("There are some NAs in the metric data. Please check.")
  }
  
  threshold_AH <- metricDF_AH %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    summarise(qntl = if (upper_percentile) quantile(metricDF_AH$metric, 1 - prct_half)[[1]] else quantile(metricDF_AH$metric, prct_half)[[1]]) %>%
    pull(qntl)
  
  threshold_R <- metricDF_R %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    summarise(qntl = if (upper_percentile) quantile(metricDF_R$metric, 1 - prct_half)[[1]] else quantile(metricDF_R$metric, prct_half)[[1]]) %>%
    pull(qntl)
  
  if (direction == 1) {
    message("Higher values refer to more climate-resilient areas.")
    metricCS_AH <- metricDF_AH %>%
      dplyr::mutate(V0_AH = ifelse(metricDF_AH$metric >= threshold_AH, 1, 0))
    metricCS_R <- metricDF_R%>%
      dplyr::mutate(V0_R = ifelse(metricDF_R$metric >= threshold_R, 1, 0))
    
  } else if (direction == -1) {
    message("Lower values refer to more climate-resilient areas.")
    metricCS_AH <- metricDF_AH %>%
      dplyr::mutate(V0_AH = ifelse(metricDF_AH$metric <= threshold_AH, 1, 0))
    metricCS_R <- metricDF_R%>%
      dplyr::mutate(V0_R = ifelse(metricDF_R$metric <= threshold_R, 1, 0))
  } else {
    stop("Please enter a valid direction: either 1 or -1.")
  }
  
  metricCS <- metricCS_AH %>%
    sf::st_drop_geometry() %>%
    dplyr::select(cellID, V0_AH) %>%
    dplyr::left_join(metricCS_R %>%
                       sf::st_drop_geometry() %>%
                       dplyr::select(cellID, V0_R), by = "cellID") %>%
    dplyr::mutate(V0 = ifelse(V0_AH == 1 | V0_R == 1, 1, 0)) %>%
    dplyr::select(cellID, V0) %>%
    sf::st_as_sf(., geometry = metricCS_AH$geometry)
  
  
  imptList <- list() # empty list to fill in with the important features
  for (i in seq_along(spp)) {
    df <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID") 
    
    imptList[[i]] <- df %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS")) := V0 * .data[[spp[i]]]) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
  }
  
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    select(cellID, everything())
  
  repList <- list()
  for (i in seq_along(spp)) {
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID")
    
    df2 <- imptList %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
    
    repList[[i]] <- df1 %>%
      dplyr::left_join(df2, by = "cellID") %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_NCS")) := ifelse(!!sym(paste0(spp[i], "_CS")) > 0, 0, .data[[spp[i]]])) %>%
      dplyr::select(matches('_NCS|_CS'))
  }
  
  repList <- repList %>%
    do.call(bind_cols, .) %>%
    dplyr::bind_cols(featuresDF %>% dplyr::select(cellID, geometry)) %>%
    dplyr::select(cellID, everything()) %>%
    sf::st_as_sf(sf_column_name = "geometry")
  
  return(repList)
}




# Input: One metricDF but not combined
fClimateFeatureOverlap_CSapproach_cont<- function(featuresDF,
                                                  percentile,
                                                  metricDF,
                                                  direction,
                                                  upper_percentile,
                                                  cut) {
  
  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  prct <- percentile / 100
  
  if (any(is.na(metricDF$metric))) {
    stop("There are some NAs in the metric data. Please check.")
  }
  
  if(cut == "YES") {
    metricDF <- metricDF %>%
      dplyr::inner_join(featuresDF %>% dplyr::select(cellID), by = "cellID")
  } else {
    metricDF <- metricDF
  }
  
  threshold <- metricDF %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    summarise(qntl = if (upper_percentile) quantile(metricDF$metric, 1 - prct)[[1]] else quantile(metricDF$metric, prct)[[1]]) %>%
    pull(qntl)
  
  if (direction == 1) {
    message("Higher values refer to more climate-resilient areas.")
    metricCS <- metricDF %>%
      dplyr::mutate(V0 = ifelse(metricDF$metric >= threshold, 1, 0))
  } else if (direction == -1) {
    message("Lower values refer to more climate-resilient areas.")
    metricCS <- metricDF %>%
      dplyr::mutate(V0 = ifelse(metricDF$metric <= threshold, 1, 0))
  } else {
    stop("Please enter a valid direction: either 1 or -1.")
  }
  
  imptList <- list() # empty list to fill in with the important features
  for (i in seq_along(spp)) {
    df <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID") 
    
    imptList[[i]] <- df %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS")) := V0 * .data[[spp[i]]]) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
  }
  
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    select(cellID, everything())
  
  repList <- list()
  for (i in seq_along(spp)) {
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID")
    
    df2 <- imptList %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
    
    repList[[i]] <- df1 %>%
      dplyr::left_join(df2, by = "cellID") %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_NCS")) := ifelse(!!sym(paste0(spp[i], "_CS")) > 0, 0, .data[[spp[i]]])) %>%
      dplyr::select(matches('_NCS|_CS'))
  }
  
  repList <- repList %>%
    do.call(bind_cols, .) %>%
    dplyr::bind_cols(featuresDF %>% dplyr::select(cellID, geometry)) %>%
    dplyr::select(cellID, everything()) %>%
    sf::st_as_sf(sf_column_name = "geometry")
  
  return(repList)
}


# Input: One metricDF but not combined
fClimateFeatureOverlap_CSapproach_bin <- function(featuresDF,
                                              percentile,
                                              metricDF,
                                              direction,
                                              upper_percentile) {
  
  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  prct <- percentile / 100
  
  if (any(is.na(metricDF$metric))) {
    stop("There are some NAs in the metric data. Please check.")
  }
  
  threshold <- metricDF %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble() %>%
    summarise(qntl = if (upper_percentile) quantile(metricDF$metric, 1 - prct)[[1]] else quantile(metricDF$metric, prct)[[1]]) %>%
    pull(qntl)
  
  if (direction == 1) {
    message("Higher values refer to more climate-resilient areas.")
    metricCS <- metricDF %>%
      dplyr::mutate(V0 = ifelse(metricDF$metric >= threshold, 1, 0))
  } else if (direction == -1) {
    message("Lower values refer to more climate-resilient areas.")
    metricCS <- metricDF %>%
      dplyr::mutate(V0 = ifelse(metricDF$metric <= threshold, 1, 0))
  } else {
    stop("Please enter a valid direction: either 1 or -1.")
  }
  
  imptList <- list() # empty list to fill in with the important features
  for (i in seq_along(spp)) {
    df <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID") 
    
    imptList[[i]] <- df %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS")) := V0 * .data[[spp[i]]]) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
  }
  
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    select(cellID, everything())
  
  repList <- list()
  for (i in seq_along(spp)) {
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metricCS, by = "cellID")
    
    df2 <- imptList %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
    
    repList[[i]] <- df1 %>%
      dplyr::left_join(df2, by = "cellID") %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_NCS")) := ifelse(!!sym(paste0(spp[i], "_CS")) == 1, 0, .data[[spp[i]]])) %>%
      dplyr::select(matches('_NCS|_CS'))
  }
  
  repList <- repList %>%
    do.call(bind_cols, .) %>%
    dplyr::bind_cols(featuresDF %>% dplyr::select(cellID, geometry)) %>%
    dplyr::select(cellID, everything()) %>%
    sf::st_as_sf(sf_column_name = "geometry")
  
  return(repList)
}


# Assign targets for CPA approach
fAssignTargets_CPA <- function(climatesmartDF,
                               targetsDF,
                               refugiaTarget){
  
  spp <- targetsDF %>%
    dplyr::select(feature) %>%
    pull()
  
  featDF <- climatesmartDF %>%
    sf::st_drop_geometry() %>%
    dplyr::select(-cellID) %>%
    dplyr::summarize(across(everything(), sum)) %>%
    tidyr::pivot_longer(everything(), names_to = "feature", values_to = "planunit")
  
  finalList <- list() # empty list
  for(i in 1:length(spp)) {
    filteredTarget <- targetsDF %>% # Get the set target per feature
      dplyr::filter(feature == spp[i])
    
    trgt <- filteredTarget$target # Extracting the target set for each feature
    
    vars <- c(stringr::str_c(spp[i], "_CS"),
              stringr::str_c(spp[i], "_NCS"))
    
    suppressMessages({
      assignTarget <- featDF %>%
        dplyr::filter(.data$feature %in% vars) %>%
        dplyr::full_join(., filteredTarget)
    })
    sumUnits <- sum(assignTarget$planunit, na.rm = TRUE) # getting the total of the feature
    
    assignTarget1 <- assignTarget %>%
      dplyr::mutate(planunit = ifelse(str_ends(feature, paste0(spp[i])), yes = sumUnits, no = planunit)) %>%
      dplyr::mutate(proportion = planunit/sumUnits)
    
    reltargetCS <- assignTarget1[assignTarget1$feature == paste0(spp[i], "_CS"), "proportion"] %>% dplyr::pull()  # get the relative target for the climate-smart areas
    
    if(reltargetCS > assignTarget1[assignTarget1$feature == spp[i], "target"] %>% dplyr::pull()) { # Do this check; is the percentile greater than the assigned target for that feature?
      targetCS <- (assignTarget1[assignTarget1$feature == spp[i], "target"] %>% as.numeric())/(assignTarget1[assignTarget1$feature == paste0(spp[i], "_CS"), "proportion"] %>% as.numeric())
      
      targetNCS <- 0
    } else {
      targetCS <- refugiaTarget
      targetNCS <- ((assignTarget1[assignTarget1$feature == spp[i], "target"] %>% as.numeric()) - reltargetCS)/(assignTarget1[assignTarget1$feature == paste0(spp[i], "_NCS"), "proportion"] %>% as.numeric())
    }
    
    finalList[[i]] <- assignTarget1 %>%
      dplyr::mutate(target = case_when(str_ends(feature, "_CS") ~ targetCS,
                                       str_ends(feature, "_NCS") ~ targetNCS)) %>%
      dplyr::filter(feature != spp[i]) %>%
      dplyr::select(feature, target)
  }
  
  finalDF <- do.call(bind_rows, finalList)
  
  return(finalDF)
  
}