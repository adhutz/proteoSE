# Antigen / patient-wise intensity plots
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.



#' Plot Antigen Intensities and Differences
#'
#' This function creates patient-wise plots of protein intensity differences for a given contrast between conditions in a `SummarizedExperiment` object. Optionally, intensities can be z-scaled. For z-score calculation, additional conditions can be included. The function also generates an overlap plot that depicts in how many patients the given proteins were above the thresholds.
#'
#' @param se A `SummarizedExperiment` object containing the experimental data.
#' @param contrast A string specifying the comparison between conditions in the format "test_vs_control".
#' @param additional_sets A character vector specifying additional condition sets to include in the analysis. Defaults to "none", meaning no additional sets are included. Set to "all" to include all conditions.
#' @param scale A logical value indicating whether to scale the data by z-transformation. Defaults to `FALSE`.
#' @param min_diff A numeric value specifying the minimum difference required for data points to be considered. Defaults to `1`.
#' @param min_intensity A numeric value specifying the minimum intensity for data points to be considered. Defaults to `10`. If scaled is set to TRUE, this needs to be changed to a lower value. 
#' @param max.overlaps A numeric value indicating the maximum number of label overlaps allowed in the plot. Defaults to `40`.
#' @param targets A data frame containing `name` and `target` columns, specifying target proteins or features of interest. Defaults to an empty data frame.
#' @param min.segment.length Minimum segment length for labels
#'
#' @return A combined plot displaying individual and overlapping data points for the specified contrast.
#' @details The function processes the `SummarizedExperiment` data to calculate the mean of control samples and the differences between test samples and control means. It generates individual plots for test samples and an overlap plot for the specified conditions.
#' 
#' If the mean for the control condition is not already present in the data, it is calculated using `proteoSE::add_stats()`.
#'
#' @examples
#' \dontrun{
#' # Assuming `se` is a valid SummarizedExperiment object and `contrast` is specified:
#' plot_antigen(se, contrast = "treatment_vs_control")
#'
#' }
#' @importFrom DEP2 get_df_wide
#' @importFrom dplyr select mutate across full_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
#' @export
plot_antigen <- function(se, contrast, additional_sets = "none", scale = FALSE, min_diff = 1, min_intensity = 10, max.overlaps = 40, targets = data.frame(name = character(0), target = character(0), min.segment.length = 0)){
  
  ctrl_condition <- gsub(".*_vs_(.*)$","\\1", contrast)
  ctrl_mean <- paste0(ctrl_condition, "_mean_intensity")
  
  test_condition <- gsub("(.*)_vs_.*","\\1", contrast)
  
  if(!ctrl_mean %in% colnames(rowData(se))){
    se <- add_stats(se, type = "mean")
    message(paste0("No mean for ctrl condition found. Mean for\"", ctrl_condition, "\" was calculated via proteoSE::add_stats(). This value is not retained in the se object."))
  }
  
  if(additional_sets == "all"){
    additional_sets <- se$condition
  } else if(additional_sets == "none"){
    additional_sets <- NULL
  }
  
  sets <- strsplit(contrast, "_vs_") %>% unlist() %>% c(., additional_sets)  %>% unique()
  
  df_diff <- se %>% DEP2::get_df_wide() 
  
  if(scale){
    temp_df <-  df_diff %>% select(matches(paste0(sets, "_[0-9]*$"))) %>% t() %>% scale() %>% t()%>% cbind(df_diff %>% select(name, contains("_diff"))) %>% mutate(across(everything() , ~ifelse(is.na(.), 0, .)))
  } else {
    temp_df <- df_diff %>% select(matches(paste0(sets, "_[0-9]*$"))) %>% cbind(df_diff %>% select(name, paste0(contrast, "_diff"))) %>% mutate(across(everything() , ~ifelse(is.na(.), 0, .)))
  }
  
  temp_df <- df_diff %>% select(name, matches(paste0(sets, "_[0-9]*$")), !!sym(ctrl_mean)) %>% 
    tidyr::pivot_longer(-c(name, !!sym(ctrl_mean)), names_to = "sample", values_to = "intensity") %>% 
    mutate(diff = intensity - !!sym(ctrl_mean)) %>% 
    tidyr::pivot_wider(names_from = "sample", values_from = c(intensity, diff)) %>% 
    select(name, starts_with("diff")) %>% 
    dplyr::full_join(temp_df, by = c("name" = "name"))
  
  
  test_samples <- grep(paste0("^", test_condition, "_[0-9]*$"), colnames(temp_df), value = TRUE)
  
  p_individuals <- lapply(test_samples, function(x) {
    plot_individuals(df = temp_df, x = paste0("diff_", x), y = x, 
                    cut_x = min_diff, cut_y = min_intensity, 
                    max.overlaps = max.overlaps, ctrl_condition = ctrl_condition, 
                    scaled = scale, min.segment.length = min.segment.length)
    }) 
  
  p_overlaps <- plot_overlap(df = temp_df, test_condition = test_condition, ctrl_condition = ctrl_condition, min_diff = min_diff, min_intensity = min_intensity, 
                             targets = targets, scaled = scale, min.segment.length = min.segment.length, max.overlaps = max.overlaps)
  
  p <- patchwork::wrap_plots(p_individuals, ncol = 2) / p_overlaps + 
    patchwork::plot_annotation(title = contrast) + patchwork::plot_layout(heights = c(1, 0.8))
}



#' Plot Antigen Missing Data Visualization
#'
#' This function creates a detailed plot to visualize missing data patterns, highlighting the percentage of missing data in test and control conditions. It also allows for the labeling of specific targets.
#'
#' @param se A `SummarizedExperiment` object containing the data to be analyzed.
#' @param test_condition A string naming the test condition to analyze.
#' @param ctrl_condition A string naming the control condition to compare against.
#' @param quantile A numeric value indicating the quantile cutoff for low intensity in the control condition. Defaults to `0.1`.
#' @param perc_low_ctrl A numeric threshold for filtering antigens based on the percentage of low-intensity measurements in the control condition. Defaults to `80`.
#' @param targets A data frame with `name` and `target` columns to highlight specific targets. Defaults to an empty data frame.
#'
#' @return A `ggplot2` object showing a plot with points representing antigens, color-coded based on their missing data percentages, and optionally labeled.
#' @details The function calculates missing and measured data percentages, applies quantile-based cutoffs for low intensity, and generates a scatter plot with facets showing different categories of missing data in the test and control conditions. The visualization helps identify patterns of data availability.
#'
#' @examples
#' \dontrun{
#' # Assuming `se` is a `SummarizedExperiment` object with relevant data:
#' plot_antigen_missing(
#'   se,
#'   test_condition = "test_cond",
#'   ctrl_condition = "control_cond",
#'   quantile       = 0.05,
#'   perc_low_ctrl  = 70
#' )
#'
#' }
#' @importFrom dplyr filter mutate group_by ungroup summarize_all summarize left_join
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_brewer facet_grid labs theme element_blank
#' @importFrom ggrepel geom_text_repel
#' @importFrom forcats fct_rev
#' @export
plot_antigen_missing <- function(se, test_condition, ctrl_condition, 
                                 quantile = 0.1, perc_low_ctrl = 80, targets = data.frame(name = character(0), target = character(0))){
  
  df_long <- DEP2::get_df_long(se[, se$condition %in% c(test_condition, ctrl_condition)]) %>% select(name, condition, intensity, label, replicate)
  n_ctrl <- sum(se$condition == ctrl_condition)
  n_test <- sum(se$condition == test_condition)
  
  df_pat <- df_long %>% 
    filter(!is.na(intensity), condition == test_condition) %>% 
    mutate(condition = "test_condition") %>%
    group_by(name) %>% 
    mutate(n = n()) %>% 
    ungroup() %>%  
    mutate(perc_miss = ((n_test-n)/n_test)*100, 
           perc_meas = ((n/n_test)*100)) %>% 
    group_by(condition, name) %>% 
    mutate(mean = mean(intensity, na.rm = TRUE)) %>% 
    ungroup() %>%   
    tidyr::pivot_wider(names_from = condition, 
                       values_from = c(intensity, perc_miss, perc_meas, mean, n))
  
  df_ctrl <- df_long %>% 
    filter(condition == ctrl_condition) %>% 
    mutate(condition = "ctrl_condition") %>%
    group_by(replicate, condition) %>% 
    mutate(cut_off = quantile(intensity, probs = quantile, na.rm = TRUE)) %>% 
    ungroup() %>%
    group_by(name, condition) %>% 
    summarize(n = sum(!is.na(intensity)), 
              n_low = sum(intensity < cut_off |is.na(intensity)), 
              n_miss = sum(is.na(intensity)),
              perc_miss = ((n_ctrl-sum(!is.na(intensity)))/n_ctrl)*100, 
              perc_low =  ((sum(intensity < cut_off |is.na(intensity)))/n_ctrl)*100, 
              mean = mean(intensity, na.rm = TRUE)) %>% 
    ungroup() %>% 
    tidyr::pivot_wider(names_from = condition, 
                       values_from = c(perc_miss, perc_low, mean, n_low, n))
  
  df_final <- inner_join(df_pat, df_ctrl, by = "name") %>% 
    #filter(perc_low_ctrl_condition >= perc_low_ctrl) %>% 
    group_by(name) %>% 
    summarize_all(mean, na.rm = TRUE) %>%
    mutate(perc_miss_ctrl_cut = cut(perc_miss_ctrl_condition, 
                                    breaks = c(0, 70, 80, 99, 100), 
                                    labels = c("0-70 %","70-80 %", "80-99 %", "100 %"), 
                                    include.lowest = TRUE)) %>% 
    mutate(x_coord = runif(nrow(.))) %>% 
    mutate(perc_low_ctrl_cut = cut(perc_low_ctrl_condition, 
                                   breaks = c(0, 80, 99, 100), 
                                   labels = c("0-80 %", "80-99 %", "100 %"), 
                                   include.lowest = TRUE, include.highest = TRUE)) %>%
    mutate(perc_miss_test_cut = cut(perc_miss_test_condition, 
                                    breaks = c(0, 10, 30, 50, 100), 
                                    labels = c("0-10 %", "10-30 %", "30-50 %", "50-100 %"), 
                                    include.lowest = TRUE))%>%
    mutate(perc_meas_test_cut = cut(perc_meas_test_condition, 
                                    breaks = c(0, 30, 50, 70, 90, 100), 
                                    labels = c("0-30 %", "30-50 %", "50-70 %", "70-90 %", "90-100 %"), 
                                    include.lowest = TRUE)) %>%
    mutate(perc_miss_ctrl_cut = forcats::fct_rev(factor(perc_miss_ctrl_cut, levels = c("0-70 %","70-80 %", "80-99 %", "100 %"))),
           perc_low_ctrl_cut = forcats::fct_rev(factor(perc_low_ctrl_cut, levels = c("0-80 %", "80-99 %", "100 %"))),
           perc_miss_test_cut = factor(perc_miss_test_cut, levels = c("0-10 %", "10-30 %", "30-50 %", "50-100 %")),
           perc_meas_test_cut = factor(perc_meas_test_cut, levels = c("0-30 %", "30-50 %", "50-70 %", "70-90 %", "90-100 %"))) %>% 
    dplyr::left_join(targets, by = "name") 
  
  set.seed(1)
  
  if(nrow(targets) == 0){
    p_final <- df_final  %>% 
      #filter(perc_meas_test_condition > 30) %>%
      ggplot(aes(x = x_coord, y = mean_test_condition, fill = perc_miss_ctrl_cut)) +
      geom_point(shape = 21, size = 3)+
      scale_fill_brewer(palette = "OrRd", direction = -1) +
      ggrepel::geom_text_repel(aes(label = name), min.segment.length = 0, max.overlaps = 10) +
      facet_grid(perc_low_ctrl_cut ~ perc_meas_test_cut, scales = "free") +
      labs(fill = paste0("Percent missing in ", ctrl_condition), y = paste0("Mean(Intensity) in ", test_condition), x = "") +
      theme_proteoSE() +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
  } else{
    p_final <- df_final  %>% 
      #filter(perc_meas_test_condition > 30) %>%
      ggplot(aes(x = x_coord, y = mean_test_condition, fill = target)) +
      geom_point(aes(shape = perc_miss_ctrl_cut), size = 3)+
      scale_fill_brewer(palette = "OrRd", direction = -1) +
      scale_shape_manual(values = c(21, 22, 24, 23, 25)) +
      ggrepel::geom_text_repel(aes(label = name), min.segment.length = 0, max.overlaps = 10) +
      facet_grid(perc_low_ctrl_cut ~ perc_meas_test_cut, scales = "free") +
      labs(fill = paste0("Percent missing in ", ctrl_condition), y = paste0("Mean(Intensity) in ", test_condition), x = "") +
      theme_proteoSE() +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  return(p_final)
}



#' Plot Individual Antigen Differences
#'
#' This function generates a scatter plot visualizing the difference between the intensity of a sample and the control mean for individual data points in a data frame. Points are colored based on whether they meet specified threshold criteria for both axes.
#'
#' @param df A data frame containing the data to be plotted.
#' @param x A string representing the column name for the x-axis variable (the difference between sample intensity and control mean).
#' @param y A string representing the column name for the y-axis variable (sample intensity).
#' @param cut_x A numeric value specifying the threshold for the x-axis variable, above which points are highlighted.
#' @param cut_y A numeric value specifying the threshold for the y-axis variable, above which points are highlighted.
#' @param max.overlaps A numeric value indicating the maximum number of overlapping labels for the text. Defaults to `10`.
#' @param ctrl_condition A string specifying the control condition name, used for labeling the plot. Defaults to `"ctrl"`.
#' @param scaled A logical value indicating whether the y-axis should be labeled as scaled. Defaults to `FALSE`.
#' @param min.segment.length Minimum length of segments for labels
#' 
#' @return A `ggplot2` object representing the scatter plot with points color-coded based on whether they meet the threshold criteria. Points above the threshold are labeled with their `name`.
#' @details The function uses `ggplot2` to create a scatter plot where points are highlighted and labeled if they meet or exceed the specified `cut_x` and `cut_y` values. Labels for points are added using `geom_text_repel` to reduce overlap.
#'
#' @examples
#' \dontrun{
#' # Assuming `df` is a valid data frame with relevant columns:
#' plot_individuals(df, x = "diff_sample1", y = "sample1", cut_x = 1, cut_y = 10)
#'
#' }
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs
#' @importFrom dplyr filter
#' @importFrom ggrepel geom_text_repel
plot_individuals <- function(df, x, y, cut_x, cut_y, max.overlaps = 10, ctrl_condition = "ctrl", scaled = FALSE, min.segment.length = 0){
  
  df <- df %>% filter(!is.na(!!sym(x)) & !is.na(!!sym(y)))
  p <- df %>% ggplot(aes(x = !!sym(x), y = !!sym(y), color = ifelse(!!sym(x) >= cut_x & !!sym(y) > cut_y, "TRUE", "FALSE"))) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("TRUE" = "darkorange", "FALSE" = "darkgrey")) +
    labs(color = "Above Thresholds", y = paste0(ifelse(scaled, "Scaled ", ""), "Intensity (", y, ")"), x = paste0(y, " - mean(", ctrl_condition, ")")) +
    geom_text_repel(data = filter(df, !!sym(x) >= cut_x & !!sym(y) >= cut_y), aes(label = name), max.overlaps = max.overlaps, min.segment.length = min.segment.length) +
    theme_proteoSE()
  return(p)
}



#' Plot Overlap of Antigen Intensity Differences
#'
#' This function creates a summary plot that visualizes the overlap of proteins based on the mean difference and mean intensity across multiple samples. It highlights the distribution of data points that meet certain criteria and optionally labels specific targets.
#'
#' @param df A data frame containing the processed data, with columns representing sample intensities and differences. Defaults to `temp_df`.
#' @param samples A character vector specifying which sample columns to include. If `NULL`, samples matching the condition pattern are selected automatically.
#' @param test_condition A string specifying the condition name used for filtering and plotting. Defaults to `"MGN_ag_neg"`.
#' @param ctrl_condition A string representing the control condition name, used for labeling the plot. Defaults to `"ctrl"`.
#' @param min_diff A numeric value specifying the minimum difference required for data points to be considered. Defaults to `0`.
#' @param min_intensity A numeric value specifying the minimum intensity for data points to be considered. Defaults to `1`.
#' @param targets A data frame with `name` and `target` columns for highlighting specific targets. Defaults to an empty data frame.
#' @param scaled A logical value indicating whether the y-axis should be labeled as scaled. Defaults to `FALSE`.
#' @param min.segment.length A numeric value specifying the minimum segment length for the text labels. Defaults to `0`.
#' @param max.overlaps A numeric value indicating the maximum number of overlapping labels for the text. Defaults to `10`.
#' 
#' @return A `ggplot2` object representing a scatter plot that shows antigen overlap based on mean differences and intensities, with options for customized labeling and color coding.
#' @details The function processes input data to compute summaries for each antigen and plots them based on the number of samples in which they are detected. Points can be color-coded by target if provided, and labeled to show individual antigens and patients.
#'
#' @examples
#' \dontrun{
#' # Assuming `df` is a valid data frame with appropriate columns:
#' plot_overlap(df, condition = "condition1", ctrl = "control1", targets = targets_df)
#'
#' }
#' @importFrom dplyr rename_with select filter group_by summarize left_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_discrete labs facet_wrap scale_fill_viridis_d
#' @importFrom ggrepel geom_text_repel
plot_overlap  <- function(df = temp_df, test_condition = "MGN_ag_neg", ctrl_condition = "ctrl", min_diff = 0, min_intensity = 1, targets = data.frame(name = character(0), target = character(0)), scaled = FALSE, samples = NULL, min.segment.length=0, max.overlaps = 10){ 
  if(is.null(samples)){
    samples <- grep(paste0("^", test_condition, "_[1-9]*$"), colnames(df), value = TRUE)
  }

  df_sum <- df %>% dplyr::rename_with(.cols = all_of(samples), ~ gsub("(.*)", "intensity_\\1", .x)) %>% 
    select(name, grep(paste0(samples, collapse = "|"), colnames(.), value = TRUE)) %>%
    tidyr::pivot_longer(-c(name), names_to = c("type", "id"), values_to = "intensity", names_pattern = "(.*?)_(.*)") %>% 
    tidyr::pivot_wider(names_from = "type", values_from = "intensity") %>%
    filter(diff >= min_diff & intensity >= min_intensity) %>% 
    group_by(name) %>% 
    summarize(n = factor(n(), levels = 1:length(samples)), 
              diff = mean(diff, na.rm = TRUE), 
              id = paste0(gsub(".*_(.*)", "\\1", id), collapse = ", "), 
              mean_intensity = mean(intensity, na.rm = TRUE)) %>% filter(n != 0)
  
  sample = paste0(".*", test_condition, "_[1-9]*$")
  
  df_sum <- df_sum %>% dplyr::left_join(targets, by = "name") 
  
  if(targets %>% nrow() > 0){
    p_summary <- df_sum %>% ggplot(aes(x = diff, y = mean_intensity, fill = target))+
      geom_point(shape = 21, size = 3, alpha = 0.8)+
      scale_fill_discrete(na.value = "grey90") +
      ggrepel::geom_text_repel(data = filter(df_sum, !n ==1 & !n == 2), aes(label = name), min.segment.length = min.segment.length, max.overlaps = max.overlaps) +
      ggrepel::geom_text_repel(data = filter(df_sum, n == 1 | n == 2), aes(label = paste0(name, " | ", id)), min.segment.length = min.segment.length, max.overlaps = max.overlaps) +
      facet_wrap(.~n, scales = "free") +
      labs(y = paste0("Mean(", ifelse(scaled, "Scaled ", ""), "Intensity)"), x = "Mean(Diff)", title = paste0(test_condition, "vs", ctrl_condition)) +
      theme_proteoSE() 
  } else {
    p_summary <- df_sum %>% ggplot(aes(x = diff, y = mean_intensity, fill = n))+
      geom_point(shape = 21, size = 3, alpha = 0.8)+
      scale_fill_viridis_d(option = "inferno", end = 0.8) +
      ggrepel::geom_text_repel(data = filter(df_sum, !n ==1 & !n == 2), aes(label = name), min.segment.length = min.segment.length, max.overlaps = max.overlaps) +
      ggrepel::geom_text_repel(data = filter(df_sum, n == 1 | n == 2), aes(label = paste0(name, " | ", id)), min.segment.length = min.segment.length, max.overlaps = max.overlaps) +
      facet_wrap(.~n, scales = "free") +
      labs(y = paste0("Mean(", ifelse(scaled, "Scaled ", ""), "Intensity)"), x = "Mean(Diff)", title = paste0(test_condition, "vs", ctrl_condition)) +
      theme_proteoSE()   
  }
  
  return(p_summary)
}

