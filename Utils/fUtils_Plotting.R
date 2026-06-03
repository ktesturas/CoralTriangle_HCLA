# Utils for Plotting

# Get Corr Matrix
fGetCorrMatrix <- function(list_plans) {
  
  y = 1
  s_matrix <- list() # empty list
  for(i in 1:length(list_plans)){
    for(j in 1:length(list_plans)){ # do pairwise comparisons for all solutions
      kappa_temp <- irr::kappa2(bind_cols(list_plans[[i]], list_plans[[j]]))
      kappa_corrvalue <- kappa_temp$value # store coefficient
      kappa_pvalue <- kappa_temp$p.value # store pvalues
      s_matrix[[y]] <- cbind(colnames(list_plans[[i]]), # first plan
                             colnames(list_plans[[j]]), # second plan
                             kappa_corrvalue, # correlation value
                             kappa_pvalue) # p value
      y = y+1
    }
  }
  
  s_matrix_all <- do.call(rbind, s_matrix) %>% 
    tibble::as_tibble()
  colnames(s_matrix_all)[1:2] <- c('plan1','plan2')
  
  matrix <- s_matrix_all %>% 
    tibble::as_tibble() %>% 
    dplyr::select(-kappa_pvalue) %>% 
    tidyr::pivot_wider(names_from = plan2, values_from = kappa_corrvalue) %>% 
    as.matrix()
  
  return(matrix)
}

# Plot Cohen's Correlation Matrix
fPlot_CorrPlot <- function(matrix, num) {
  # creating corrplot
  rownames(matrix) <- matrix[,1]
  n <- num + 1 # num represents the number of inputted spatial plans
  matrix_f <- matrix[,2:n]
  class(matrix_f) <- "numeric"
  
  #col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  plot <- corrplot(matrix_f, method = "shade", cl.lim = c(-0.2,1), tl.col = "black", addCoef.col = "black",
                   col = COL2('BrBG', 200), tl.srt=45)
  return(plot)
}


# Create object for selection frequency across a set of solutions
fGetSelFrequency <- function(solution_list, 
                             col_names, 
                             PUs
) {
  
  df <- list() # empty list
  for (i in 1:length(col_names)) {
    df[[i]] <- solution_list[[i]] %>% 
      dplyr::select(solution_1, cellID) %>% # Select solution column
      dplyr::rename(!!sym(col_names[i]) := solution_1) %>%  # Rename solution column to the column name
      tibble::as_tibble()
    
  }
  
  tmp <- purrr::reduce(df, dplyr::left_join) %>% # join all the solutions
    dplyr::select(-geometry) %>%
    mutate(selection = rowSums(select(., -cellID))) # Calculate selection frequency for all PUs
  #dplyr::filter(selection != 0) # Exclude rows with 0 values in the selection column
  
  
  
  # Create the selection frequency sf object
  low_regret <- dplyr::full_join(tmp, PUs, 
                                 by = "cellID") %>% 
    sf::st_as_sf(sf_column_name = "geometry") %>% 
    dplyr::left_join(., 
                     UniformCost %>% 
                       tibble::as_tibble(), 
                     by = c("geometry", "cellID")) %>% 
    sf::st_as_sf(sf_column_name = "geometry") 
  
  return(low_regret)
}

# Create selection frequency plot
fPlot_SelFrequency <- function(data, land) {
  gg <- ggplot() + geom_sf(data = data, aes(fill = as.factor(selection)), color = NA, size = 0.01) +
    geom_sf(data = land, color = "grey20", fill = "grey20", alpha = 0.9, size = 0.1, show.legend = FALSE) +
    #   geom_sf(data = boundary, color = "black", fill = NA, size = 0.1, show.legend = FALSE) +
    coord_sf(xlim = st_bbox(data)$xlim, ylim = st_bbox(data)$ylim) +
    scale_fill_brewer(name = "Selection Frequency",
                      palette = "PuBu", aesthetics = "fill") +
    theme_bw() +
    theme(axis.ticks = element_line(color = "black", linewidth = 1),
          panel.border = element_rect(colour = "black", 
                                      fill=NA, 
                                      linewidth=1),
          axis.text = element_text(size = 10, color = "black"))
  return(gg)
}


# Mqke inset histogram
plot_inset <- function(sFreq) {
  
  temp <- sFreq %>% 
    tibble::as_tibble() %>% 
    dplyr::select(selection) %>% # Select the selection frequency column
    dplyr::group_by(selection) %>% 
    dplyr::summarize(total = n()) %>% # Summarize how many planning units have that selection frequency
    dplyr::mutate(proportion = total/nrow(PUs)) %>% 
    #arrange(., desc(selection)) %>%  # commenting this out because want to get areas with exact selection frequencies
    #dplyr::mutate(proportion = ifelse(selection == 0, yes = total/nrow(PUs), no = cumsum(total)/nrow(PUs))) %>% 
    dplyr::arrange(selection) # Ascending value of selection
  
  # Creating the inset plot
  inset <- ggplot(temp, 
                  aes(x = as.factor(selection), 
                      y = proportion, 
                      fill = as.factor(selection))) +
    ggplot2::scale_fill_brewer(name = "Selection Frequency",
                               palette = "PuBu", aesthetics = "fill") +
    ggplot2::geom_col(width = 1, show.legend = FALSE) +
    ggplot2::theme_bw() + 
    ggplot2::xlab(element_blank()) +
    ggplot2::ylab(element_blank()) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::labs(title = element_blank()) +
    ggplot2::theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())
  
  return(inset)
}


# Ceate kernel density plots
create_kernel_plot <- function(data, x, x_label, fill_label, filename, fill_option = "A") {
  plot <- ggplot() +
    geom_density_ridges_gradient(data = data %>% dplyr::filter(solution_1 == "1"), 
                                 aes_string(x = x, y = "approach", fill = "after_stat(x)"), 
                                 quantile_lines = TRUE, quantile_fun = mean,
                                 scale = 4, show.legend = TRUE) +
    scale_fill_viridis_c(name = fill_label, option = fill_option) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0))) +
    labs(x = x_label) +
    theme_classic() +
    theme(axis.ticks = element_line(color = "black", linewidth = 1),
          axis.line = element_line(colour = "black", linewidth = 1),
          axis.text = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_blank(),
          legend.key.height = unit(1, "inch"),
          legend.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 15, color = "black"))
  
  ggsave(plot, filename = here("figures", filename), width = 27, height = 20)
  return(plot)
}



# Ceate kernel density plots
create_kernel_plot2 <- function(data, x, x_label, fill_label, filename, fill_option = "A", legend_text_expr) {
  plot <- ggplot() +
    geom_density_ridges_gradient(data = data %>% dplyr::filter(solution_1 == "1"), 
                                 aes_string(x = x, y = "approach", fill = "after_stat(x)"), 
                                 quantile_lines = TRUE, quantile_fun = mean,
                                 scale = 1, show.legend = TRUE) +
    scale_fill_viridis_c(name = legend_text_expr, option = fill_option) +
    geom_density_ridges(data = data %>% dplyr::filter(solution_1 == "0"), 
                        aes_string(x = x, y = "approach"), 
                        quantile_lines = TRUE, quantile_fun = mean,
                        alpha = 0.2, color = "black",
                        linetype = "dashed", 
                        scale = 1, 
                        show.legend = TRUE)  +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0))) +
    labs(x = x_label) +
    theme_classic() +
    theme(axis.ticks.x = element_line(color = "black", linewidth = 1),
          axis.ticks.y = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 1),
          axis.text = element_text(color = "black", size = 30),
          axis.title.x = element_text(size = 35, face = "bold"),
          axis.title.y = element_blank(),
          legend.key.height = unit(1.5, "inch"),
          legend.text = element_text(size = 25, color = "black"),
          legend.title = element_text(size = 30, face = "bold", color = "black"))
  
  return(plot)
}


create_kernel_plot3 <- function(data, x, x_label, fill_label, filename, fill_option = "A", legend_text_expr) {
  plot <- ggplot() +
    geom_density_ridges_gradient(data = data %>% dplyr::filter(solution_1 == "1"), 
                                 aes_string(x = x, y = "approach", fill = "after_stat(x)"), 
                                 quantile_lines = TRUE, quantile_fun = mean,
                                 scale = 1, show.legend = TRUE) +
    scale_fill_viridis_c(name = legend_text_expr, option = fill_option) +
    geom_density_ridges(data = data %>% dplyr::filter(solution_1 == "0"), 
                        aes_string(x = x, y = "approach"), 
                        quantile_lines = TRUE, quantile_fun = mean,
                        alpha = 0.2, color = "black",
                        linetype = "dashed", 
                        scale = 1, 
                        show.legend = TRUE)  +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0))) +
    labs(x = x_label) +
    theme_classic() +
    theme(axis.ticks.x = element_line(color = "black", linewidth = 1),
          axis.ticks.y = element_blank(),  # Corrected this line
          axis.line = element_line(colour = "black", linewidth = 1),
          axis.text = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_blank(),
       #   axis.text.y = element_blank(),  # Hides the y-axis text labels
          legend.key.height = unit(1, "inch"),
          legend.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 20, face = "bold", color = "black"))
  
  ggsave(plot, filename = here("figures", filename), width = 27, height = 20)
  return(plot)
}




# Plot Targets Kernel
fPlot_RidgeTargetApproach <- function(datasets, approaches) {
  combined_df <- do.call(rbind, datasets) %>%
    dplyr::mutate(approach = rep(approaches, each = nrow(datasets[[1]]))) %>%
    dplyr::mutate(approach = fct_relevel(approach, 
                                         "CPA Overlap Difference", 
                                         "CPA Overlap Hotspot Refugia", 
                                         "CPA Overlap Refugia", 
                                         "CPA Overlap Hotspot", 
                                         "Climate-uninformed"))
  
  gg <- ggplot(data = combined_df) +
    geom_density_ridges(aes(x = PercentProtection, y = approach, group = approach, fill = approach), scale = 2) +
    scale_fill_manual(values = c(`CPA Overlap Difference` = "#FF7F0E", 
                                 `CPA Overlap Hotspot Refugia` = "#1F77B4", 
                                 `CPA Overlap Refugia` = "#2CA02C", 
                                 `CPA Overlap Hotspot` = "#9467BD", 
                                 `Climate-uninformed` = "#E377C2")) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "red", linewidth = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0))) +
    labs(x = "Protection (%)", y = "Selection") +
    theme_classic() +
    theme(axis.ticks = element_line(color = "black", linewidth = 1),
          axis.line = element_line(colour = "black", linewidth = 1),
          axis.text.x = element_text(color = "black", size = 20),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_blank(),
          legend.position = "none")
  
  return(gg)
}


create_kernel_plot_single <- function(data, x, x_label, fill_label, filename, fill_option = "A", legend_text_expr) {
  plot <- ggplot() +
    geom_density(data = data %>% dplyr::filter(solution_1 == "1"),
                 aes_string(x = x, fill = "after_stat(x)"), 
                 alpha = 0.8) +
    scale_fill_viridis_c(name = legend_text_expr, option = fill_option) +
    geom_density(data = data %>% dplyr::filter(solution_1 == "0"),
                 aes_string(x = x), 
                 linetype = "dashed", color = "black", alpha = 0.4) +
    labs(x = x_label) +
    theme_classic() +
    theme(axis.ticks.x = element_line(color = "black", linewidth = 1),
          axis.line = element_line(colour = "black", linewidth = 1),
          axis.text = element_text(color = "black", size = 45),
          axis.title.x = element_text(size = 50, face = "bold"),
          axis.title.y = element_blank(),
          legend.key.height = unit(1.5, "inch"),
          legend.text = element_text(size = 40, color = "black"),
          legend.title = element_text(size = 45, face = "bold", color = "black"))
  
  ggsave(plot, filename = here("figures", filename), width = 27, height = 15)
  return(plot)
}