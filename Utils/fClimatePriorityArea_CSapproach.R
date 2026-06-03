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


fClimatePriorityArea_CSapproach_Combined2_bin <- function(featuresDF, 
                                                          percentile, 
                                                          metricDF, 
                                                          direction) {
  
  # Extract species names
  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  # Prepare the climate metric table
  metric <- metricDF %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble()
  
  imptList <- list() # Initialize empty list
  
  for(i in seq_along(spp)) {
    
    # Prepare dataframe for each species
    df <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric, by = "cellID") %>%
      na.omit()
    
    # Check for NA values
    if(any(is.na(df$metric))) {
      stop("There are some NAs in the metric data. Please check.")
    }
    
    # Filter for species presence (binary: presence = 1)
    filteredDF <- df %>%
      dplyr::filter(!!sym(spp[i]) == 1)
    
    # Calculate split percentiles
    prct_half <- (percentile / 2) / 100
    
    threshold_upper <- filteredDF %>%
      summarise(qntl = quantile(metric, 1 - prct_half)[[1]]) %>%
      pull(qntl)
    
    threshold_lower <- filteredDF %>%
      summarise(qntl = quantile(metric, prct_half)[[1]]) %>%
      pull(qntl)
    
    # Select cells based on climate thresholds
    if (direction == 1) {
      if (i == 1) message("Higher values refer to more climate-resilient areas.")
      df <- df %>%
        dplyr::mutate(V1 = ifelse(metric >= threshold_upper | metric <= threshold_lower, 1, 0),
                      V2 = ifelse(!!sym(spp[i]) == 1, 1, 0))
    } else if (direction == -1) {
      if (i == 1) message("Lower values refer to more climate-resilient areas.")
      df <- df %>%
        dplyr::mutate(V1 = ifelse(metric <= threshold_lower | metric >= threshold_upper, 1, 0),
                      V2 = ifelse(!!sym(spp[i]) == 1, 1, 0))
    } else {
      stop("Please enter a valid direction: either 1 or -1.")
    }
    
    # Create binary climate-resilient species layer
    imptList[[i]] <- df %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS")) := V1 * V2) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
  }
  
  # Combine all species layers
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    dplyr::select(cellID, everything())
  
  repList <- list()
  
  for(i in seq_along(spp)) {
    
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric, by = "cellID")
    
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
    sf::st_as_sf(sf_column_name = "geometry") %>%
    na.omit()
  
  return(repList)
}

fClimatePriorityArea_CSapproach_Combined_bin <- function(featuresDF,
                                                         percentile,
                                                         metricDF_R,
                                                         metricDF_AH,
                                                         direction) {
  
  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  metric_R <- metricDF_R %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble()
  
  metric_AH <- metricDF_AH %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble()
  
  imptList <- list()
  
  for(i in seq_along(spp)) {
    
    df_R <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric_R, by = "cellID") %>%
      na.omit()
    
    df_AH <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric_AH, by = "cellID") %>%
      na.omit()
    
    if(any(is.na(df_R$metric)) || any(is.na(df_AH$metric))) {
      stop("There are some NAs in the metric data. Please check.")
    }
    
    filteredDF_R <- df_R %>%
      dplyr::filter(!!sym(spp[i]) == 1)
    
    filteredDF_AH <- df_AH %>%
      dplyr::filter(!!sym(spp[i]) == 1)
    
    prct_half <- (percentile / 2) / 100
    
    qntl_R <- filteredDF_R %>%
      summarise(qntl = quantile(metric, 1 - prct_half)[[1]])
    
    qntl_AH <- filteredDF_AH %>%
      summarise(qntl = quantile(metric, 1 - prct_half)[[1]])
    
    if (direction == 1) {
      if (i == 1) message("Higher values refer to more climate-resilient areas.")
      
      df_R <- df_R %>%
        dplyr::mutate(V1_R = ifelse(metric >= qntl_R$qntl, 1, 0),
                      V2_R = ifelse(!!sym(spp[i]) == 1, 1, 0))
      
      df_AH <- df_AH %>%
        dplyr::mutate(V1_AH = ifelse(metric >= qntl_AH$qntl, 1, 0),
                      V2_AH = ifelse(!!sym(spp[i]) == 1, 1, 0))
      
    } else if (direction == -1) {
      if (i == 1) message("Lower values refer to more climate-resilient areas.")
      
      df_R <- df_R %>%
        dplyr::mutate(V1_R = ifelse(metric <= qntl_R$qntl, 1, 0),
                      V2_R = ifelse(!!sym(spp[i]) == 1, 1, 0))
      
      df_AH <- df_AH %>%
        dplyr::mutate(V1_AH = ifelse(metric <= qntl_AH$qntl, 1, 0),
                      V2_AH = ifelse(!!sym(spp[i]) == 1, 1, 0))
    } else {
      stop("Please enter a valid direction: either 1 or -1.")
    }
    
    imptList[[i]] <- df_R %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS_R")) := V1_R * V2_R) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS_R")), cellID) %>%
      dplyr::left_join(df_AH %>%
                         dplyr::mutate(!!sym(paste0(spp[i], "_CS_AH")) := V1_AH * V2_AH) %>%
                         dplyr::select(!!sym(paste0(spp[i], "_CS_AH")), cellID), by = "cellID")
  }
  
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    dplyr::select(cellID, everything())
  
  repList <- list()
  
  for(i in seq_along(spp)) {
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID)
    
    df2 <- imptList %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS_R")), !!sym(paste0(spp[i], "_CS_AH")), cellID)
    
    repList[[i]] <- df1 %>%
      dplyr::left_join(df2, by = "cellID") %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_NCS")) := ifelse(!!sym(paste0(spp[i], "_CS_R")) == 1 | !!sym(paste0(spp[i], "_CS_AH")) == 1, 0, .data[[spp[i]]])) %>%
      dplyr::select(matches('_NCS|_CS'))
  }
  
  repList <- repList %>%
    do.call(bind_cols, .) %>%
    dplyr::bind_cols(featuresDF %>% dplyr::select(cellID, geometry)) %>%
    dplyr::select(cellID, everything()) %>%
    sf::st_as_sf(sf_column_name = "geometry") %>%
    na.omit()
  
  return(repList)
}




# One metric DF input
# Updated Function 1
fClimatePriorityArea_CSapproach_Combined2 <- function(featuresDF, 
                                                      percentile, 
                                                      metricDF, 
                                                      direction) {
  
  # Convert featuresDF to tibble and extract species names
  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  # Drop geometry from metricDF and convert to tibble
  metric <- metricDF %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble()
  
  imptList <- list() # Initialize an empty list to store important features
  

  for(i in seq_along(spp)) {
    # Prepare dataframe for each species
    df <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric, by = "cellID") %>%
      na.omit()
    
    # Check for NAs in the metric column
    if(any(is.na(df$metric))) {
      stop("There are some NAs in the metric data. Please check.")
    }
    
    # Filter dataframe for positive species values
    filteredDF <- df %>%
      dplyr::filter(!!sym(spp[i]) > 0)
    
    # Calculate split percentiles
    prct_half <- (percentile / 2) / 100
    
    threshold_upper <- filteredDF %>%
      summarise(qntl = quantile(metric, 1 - prct_half)[[1]]) %>%
      pull(qntl)
    
    threshold_lower <- filteredDF %>%
      summarise(qntl = quantile(metric, prct_half)[[1]]) %>%
      pull(qntl)
    
    # Assign values based on the direction parameter
    if (direction == 1) {
      if (i == 1) {
        message("Higher values refer to more climate-resilient areas.")
      }
      df <- df %>%
        dplyr::mutate(V1 = ifelse(metric >= threshold_upper | metric <= threshold_lower, 1, 0))
    } else if (direction == -1) {
      if (i == 1) {
        message("Lower values refer to more climate-resilient areas.")
      }
      df <- df %>%
        dplyr::mutate(V1 = ifelse(metric <= threshold_lower | metric >= threshold_upper, 1, 0))
    } else {
      stop("Please enter a valid direction: either 1 or -1.")
    }
    
    # Create a new column for the climate resilience score
    imptList[[i]] <- df %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS")) := V1 * !!sym(spp[i])) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
  }
  
  
  # Combine all individual species data into one dataframe
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    select(cellID, everything())
  
  repList <- list()
  for(i in seq_along(spp)) {
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric, by = "cellID")
    
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
    sf::st_as_sf(sf_column_name = "geometry") %>%
    na.omit()
  
  return(repList)
}




# Two metric DF inputs
fClimatePriorityArea_CSapproach_Combined <- function(featuresDF,
                                            percentile,
                                            metricDF_R,
                                            metricDF_AH,
                                            direction) {

  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  metric_R <- metricDF_R %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble()
  
  metric_AH <- metricDF_AH %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble()
  
  imptList <- list() # empty list to fill in with the important features
  for(i in seq_along(spp)) {
    
    df_R <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric_R, by = "cellID") %>%
      na.omit()
    
    df_AH <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric_AH, by = "cellID") %>%
      na.omit()
    
    if(any(is.na(df_R$metric)) || any(is.na(df_AH$metric))) {
      stop("There are some NAs in the metric data. Please check.")
    }
    
    
    filteredDF_R <- df_R %>%
      dplyr::filter(!!sym(spp[i]) == 1)
    
    filteredDF_AH <- df_AH %>%
      dplyr::filter(!!sym(spp[i]) == 1)
    
    prct_half <- (percentile / 2) /100
    
    qntl_R <- filteredDF_R %>%
      summarise(qntl = quantile(metric, 1- prct_half)[[1]])
    
    qntl_AH <- filteredDF_AH %>%
      summarise(qntl = quantile(metric, 1 - prct_half)[[1]])
    
    if (direction == 1) {
      if (i == 1) {
        message("Higher values refer to more climate-resilient areas.")
      }
      df_R <- df_R %>%
        dplyr::mutate(V1_R = ifelse(metric >= qntl_R$qntl, 1, 0),
                      V2_R = ifelse(!!sym(spp[i]) == 1, 1, 0))
      
      df_AH <- df_AH %>%
        dplyr::mutate(V1_AH = ifelse(metric >= qntl_AH$qntl, 1, 0),
                      V2_AH = ifelse(!!sym(spp[i]) == 1, 1, 0))
    } else if (direction == -1) {
      if (i == 1) {
        message("Lower values refer to more climate-resilient areas.")
      }
      df_R <- df_R %>%
        dplyr::mutate(V1_R = ifelse(metric <= qntl_R$qntl, 1, 0),
                      V2_R = ifelse(!!sym(spp[i]) == 1, 1, 0))
      
      df_AH <- df_AH %>%
        dplyr::mutate(V1_AH = ifelse(metric <= qntl_AH$qntl, 1, 0),
                      V2_AH = ifelse(!!sym(spp[i]) == 1, 1, 0))
    } else {
      stop("Please enter a valid direction: either 1 or -1.")
    }
    
    imptList[[i]] <- df_R %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS_R")) := V1_R * V2_R) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS_R")), cellID) %>%
      dplyr::left_join(df_AH %>%
                         dplyr::mutate(!!sym(paste0(spp[i], "_CS_AH")) := V1_AH * V2_AH) %>%
                         dplyr::select(!!sym(paste0(spp[i], "_CS_AH")), cellID), by = "cellID")
  }

imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
  select(cellID, everything())

repList <- list()
for(i in seq_along(spp)) {
  
  df1 <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(!!sym(spp[i]), cellID) %>%
    dplyr::left_join(imptList %>%
                       dplyr::select(!!sym(paste0(spp[i], "_CS_R")), !!sym(paste0(spp[i], "_CS_AH")), cellID), by = "cellID")
  
  repList[[i]] <- df1 %>%
    dplyr::mutate(!!sym(paste0(spp[i], "_NCS")) := ifelse(!!sym(paste0(spp[i], "_CS_R")) == 1 | !!sym(paste0(spp[i], "_CS_AH")) == 1, 0, .data[[spp[i]]])) %>%
    dplyr::select(matches('_NCS|_CS'))
}

repList <- repList %>%
  do.call(bind_cols, .) %>%
  dplyr::bind_cols(featuresDF %>% dplyr::select(cellID, geometry)) %>%
  dplyr::select(cellID, everything()) %>%
  sf::st_as_sf(sf_column_name = "geometry") %>%
  na.omit()

return(repList)
}



fAssignTargets_CPA_Combined <- function(climatesmartDF,
                               targetsDF,
                               csTarget){
  
  
  targetsDF <- out_spp_C_bin %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names() %>%
    tibble::as_tibble() %>%
    dplyr::rename(feature = value) %>%
    dplyr::mutate(target = 0.3) # 30% protection for each species
  
  
  
  spp <- targetsDF %>%
    dplyr::select(feature) %>%
    pull()
  
  featDF <- climatesmartDF %>%
    sf::st_drop_geometry() %>%
    dplyr::select(-cellID) %>%
    dplyr::summarize(across(everything(), sum)) %>%
    tidyr::pivot_longer(everything(), names_to = "feature", values_to = "planunit")
  
  finalList <- list() # empty list
  for(i in seq_along(spp)) {
    filteredTarget <- targetsDF %>% # Get the set target per feature
      dplyr::filter(feature == spp[i])
    
    trgt <- filteredTarget$target # Extracting the target set for each feature
    
    vars <- c(stringr::str_c(spp[i], "_CS_R"),
              stringr::str_c(spp[i], "_CS_AH"),
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
    
    reltargetCS_R <- assignTarget1[assignTarget1$feature == paste0(spp[i], "_CS_R"), "proportion"] %>% dplyr::pull()
    reltargetCS_AH <- assignTarget1[assignTarget1$feature == paste0(spp[i], "_CS_AH"), "proportion"] %>% dplyr::pull()
    
    
    if(reltargetCS_R + reltargetCS_AH > assignTarget1[assignTarget1$feature == spp[i], "target"] %>% dplyr::pull()) {
      targetCS_R <- (assignTarget1[assignTarget1$feature == spp[i], "target"] %>% as.numeric()) / (reltargetCS_R + reltargetCS_AH) * reltargetCS_R
      targetCS_AH <- (assignTarget1[assignTarget1$feature == spp[i], "target"] %>% as.numeric()) / (reltargetCS_R + reltargetCS_AH) * reltargetCS_AH
      targetNCS <- 0
    } else {
      targetCS_R <- csTarget
      targetCS_AH <- csTarget
      targetNCS <- ((assignTarget1[assignTarget1$feature == spp[i], "target"] %>% 
                       as.numeric()) - (reltargetCS_R + reltargetCS_AH)) / assignTarget1[assignTarget1$feature == paste0(spp[i], "_NCS"), "proportion"] %>% as.numeric()
    }
    
    finalList[[i]] <- assignTarget1 %>%
      dplyr::mutate(target = case_when(
        str_ends(feature, "_CS_R") ~ targetCS_R,
        str_ends(feature, "_CS_AH") ~ targetCS_AH,
        str_ends(feature, "_NCS") ~ targetNCS
      )) %>%
      dplyr::filter(feature != spp[i]) %>%
      dplyr::select(feature, target)
  }
  
  finalDF <- do.call(bind_rows, finalList)
  return(finalDF)

}










########
fClimatePriorityArea_CSapproach_bin <- function(featuresDF,
                                            percentile,
                                            metricDF,
                                            direction
                                            
) {
  
  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  metric <- metricDF %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble()
  
  imptList <- list() # empty list to fill in with the important features
  for(i in seq_along(spp)) {
    
    df <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric, by = "cellID") %>%
      na.omit()
    
    if(any(is.na(df$metric))) {
      stop("There are some NAs in the metric data. Please check.")
    }
    
    filteredDF <- df %>%
      dplyr::filter(!!sym(spp[i]) == 1)
    
    prct <- percentile / 100
    
    qntl <- filteredDF %>%
      summarise(qntl = quantile(metric, 1 - prct)[[1]])
    
    if (direction == 1) {
      if (i == 1) {
        message("Higher values refer to more climate-resilient areas.")
      }
      df <- df %>%
        dplyr::mutate(V1 = ifelse(metric >= qntl$qntl, 1, 0),
                      V2 = ifelse(!!sym(spp[i]) == 1, 1, 0))
    } else if (direction == -1) {
      if (i == 1) {
        message("Lower values refer to more climate-resilient areas.")
      }
      df <- df %>%
        dplyr::mutate(V1 = ifelse(metric <= qntl$qntl, 1, 0),
                      V2 = ifelse(!!sym(spp[i]) == 1, 1, 0))
    } else {
      stop("Please enter a valid direction: either 1 or -1.")
    }
    
    imptList[[i]] <- df %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS")) := V1 * V2) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
  }
  
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    select(cellID, everything())
  
  
  repList <- list()
  for(i in seq_along(spp)) {
    
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric, by = "cellID")
    
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
    sf::st_as_sf(sf_column_name = "geometry") %>%
    na.omit()
  
  return(repList)
}

fClimatePriorityArea_CSapproach_cont <- function(featuresDF,
                                                 percentile,
                                                 metricDF,
                                                 direction) {
  # Convert featuresDF to tibble and extract species names
  spp <- featuresDF %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry, -cellID) %>%
    names()
  
  # Drop geometry from metricDF and convert to tibble
  metric <- metricDF %>%
    sf::st_drop_geometry() %>%
    tibble::as_tibble()
  
  imptList <- list() # Initialize an empty list to store important features
  
  for(i in seq_along(spp)) {
    # Prepare dataframe for each species
    df <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric, by = "cellID") %>%
      na.omit()
    
    # Check for NAs in the metric column
    if(any(is.na(df$metric))) {
      stop("There are some NAs in the metric data. Please check.")
    }
    
    # Filter dataframe for positive species values
    filteredDF <- df %>%
      dplyr::filter(!!sym(spp[i]) > 0)
    
    # Calculate the quantile threshold
    prct <- percentile / 100
    qntl <- filteredDF %>%
      summarise(qntl = quantile(metric, 1 - prct)[[1]])
    
    # Assign values based on the direction parameter
    if (direction == 1) {
      if (i == 1) {
        message("Higher values refer to more climate-resilient areas.")
      }
      df <- df %>%
        dplyr::mutate(V1 = ifelse(metric >= qntl$qntl, 1, 0))
    } else if (direction == -1) {
      if (i == 1) {
        message("Lower values refer to more climate-resilient areas.")
      }
      df <- df %>%
        dplyr::mutate(V1 = ifelse(metric <= qntl$qntl, 1, 0))
    } else {
      stop("Please enter a valid direction: either 1 or -1.")
    }
    
    # Create a new column for the climate resilience score
    imptList[[i]] <- df %>%
      dplyr::mutate(!!sym(paste0(spp[i], "_CS")) := V1 * !!sym(spp[i])) %>%
      dplyr::select(!!sym(paste0(spp[i], "_CS")), cellID)
  }
  
  # Combine all individual species data into one dataframe
  imptList <- purrr::reduce(imptList, function(x, y) left_join(x, y, by = "cellID")) %>%
    select(cellID, everything())
  
  
  repList <- list()
  for(i in seq_along(spp)) {
    
    df1 <- featuresDF %>%
      tibble::as_tibble() %>%
      dplyr::select(!!sym(spp[i]), cellID) %>%
      dplyr::left_join(metric, by = "cellID")
    
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
    sf::st_as_sf(sf_column_name = "geometry") %>%
    na.omit()
  
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