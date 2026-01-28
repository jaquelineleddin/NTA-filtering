#' @importFrom stats sd
#' @importFrom utils read.csv write.table type.convert
NULL



#' Load feature list after preprocessing with MSDial
#'
#' Loads a csv or txt file as a data frame and drops the columns
#' MS1.isotopic.spectrum and MS.MS.spectrum.
#'
#' @param filepath path of a csv or txt file
#' @param skip number of rows to skip before the data starts. Default value is 4 (for data of MSDial).
#' @param sep separator of the columns. Default value is a tabstopp (for data of MSDial).
#' @param drop_last_cols number of last columns in the data frame to drop. Default value is 2.
#' @param drop_cols names of columns to delete
#' @param dec decimal separator
#'
#' @return a data frame
#' @export
#'
#' @examples

load_data <- function(filepath, skip = 4, sep="\t", drop_last_cols = 2,
                      drop_cols = c("MS1.isotopic.spectrum", "MS.MS.spectrum" ), dec = "."){

  df_data <- read.csv(filepath,
                      header = TRUE,
                      skip = skip,
                      sep = sep,
                      stringsAsFactors=FALSE,
                      na.strings = "null",
                      dec = dec)

  # Delete columns that are not needed
  if(length(drop_cols) != 0){
  df_data <- df_data[ , -which(names(df_data) %in% drop_cols)]
  }
  df <- df_data[1: (ncol(df_data) - drop_last_cols)]

  # replace True and False with logical
  df[df == "True"] <- TRUE
  df[df == "False"] <- FALSE

  df <- type.convert(df, as.is = TRUE)

  return(df)
}


#' Remove -1 from data frame with ion mobility
#'
#' @param df feature list exported from MS-Dial with CCS values
#'
#' @returns data frame without rows of negative CCS values
#' @export
#'
#' @examples
remove_neg_mobility <- function(df){

  if('Average.CCS' %in% colnames(df)){

    df_filtered <- df[which(df$Average.CCS > 0) , ]
    return(df_filtered)

  } else {
    print("Cannot find column of ion mobility. The column name must be Average CCS!")
    return(df)
    }
}

#' Get column names of all samples
#'
#' @param df data frame
#' @param first_sample_col name of the first column of samples
#' @param number_of_samples number of all sample incl. blanks and quality controls
#'
#' @return list of column names of samples incl. blanks and quality controls
#' @export
#'
#' @examples
get_all_sample_cols <- function(df, first_sample_col, number_of_samples){

  first_index <- which(colnames(df) == first_sample_col)

  samples <- colnames(df)[first_index:(first_index+number_of_samples-1)]

  return(samples)
}


#' Get column names of specific samples
#'
#' @param all_sample_cols list of all sample columns in data frame
#' @param target_word Pattern to get specific columns, e.g. "Blank" or "QC" for quality control.
#' Use NULL, to return an empty list
#'
#' @return list of column names
#' @export
#'
#' @examples
get_cols <- function(all_sample_cols, target_word = NULL){

  if(is.null(target_word)){
    res <- c()
  }
  else {

    res <- grep(target_word, all_sample_cols, value = TRUE)

    if(length(res) == 0){
      stop(target_word, " is in no part of a column name. Please check colnames in data frame."  )
    }
  }

  return(res)
}



#' replace negative intensities or zero with NA
#'
#' If there are negative intensities or zero in the samples, QC or blanks,
#' this function will replace them with NA (not available, so they are missing values).
#'
#' @param df data frame (e.g. the data frame,
#' that function load_data of the NTAWorkflow Package returns)
#' @param all_sample_cols list of column names of all samples
#'
#' @return data frame
#' @export
#'
#' @examples
replace_missing_values <- function(df, all_sample_cols){

  # replace 0 or negative values with NA
  df_new <- df
  df_new[all_sample_cols][df_new[all_sample_cols] <= 0] <- NA

  return(df_new)
}


#' build a vector of powers of 10
#'
#' function builds a sequence of different powers of 10 for the ticks
#' of logarithmic y-axis in the splots (ordered maximum intensity of samples)
#'
#' @param max_value maximum value in the sample columns in the data frame of raw data
#'
#' @return sequence of 10^i until the maximum value is reached.
#' @noRd
#'
#' @examples breaks <- breaks_for_splot(5000)
#' ## returns: c(10,100,1000,10000)
breaks_for_splot <- function(max_value){

  max_exp10 <- ceiling(log10(max_value))

  breaks <- rep(NA, max_exp10)

  for (i in 1:max_exp10) {
    breaks[i] = 10^i
  }

  return(breaks)
}


#' plot ordered maximum intensity of measured features
#'
#' function calculates the maximum intensity over all samples for every feature and plot
#' them with logarithmic y-axis and a subplot of the last 50% with line for 2, 5, 10% of the
#' non-zero values. These lines in the subplot define values (without NA) in the samples,
#' which will be lost by the threshold on the y-axis.
#'
#' @param df data frame with feature list after
#' @param sample_cols list of column names of samples without blanks and quality controls
#' @param fn_out optional export plot as png
#'
#'
#' @return list of data frame (df), plot and subplot
#' @export
#'
#' @examples
calc_max_intensity <- function(df, sample_cols, fn_out = NULL){

  df_samples <- df[sample_cols]


  # create array with maximum intensity for every feature
  max_intensity <- rep(NA, nrow(df_samples))

  for (i in 1:nrow(df_samples)){

    temp_row <- df_samples[i, ]

    if(sum(!is.na(temp_row)) > 0){
      max_intensity[i] <- max(temp_row, na.rm = TRUE)
    } else{
    max_intensity[i] <- NA
    }
  }

  # sort array with the maximum intensities
  max_sorted <- sort( max_intensity, decreasing = TRUE, na.last = NA)

  # calculate major ticks for logarithmic plot
  breaks <- breaks_for_splot(max(max_sorted))

  # plot
  splot <- ggplot()+
    geom_point(aes(x = 1:length(max_sorted), y = max_sorted), size = 0.5)+
    scale_y_log10( breaks = breaks) +
    theme_bw()+
    labs(x = "Features without QC and Blank",
         y = "max Intensity")+
    theme(axis.text = element_text(size = 11))


  # splot with 50% of the features and marker for the last 10 % and 5% and 2% that are not zero

  num_all_features <- length(max_sorted)
  subset_intensity <- max_sorted[ceiling(num_all_features / 2) : num_all_features]

  last_10_ind <- ceiling(num_all_features * 0.9)
  last_10_intens <- max_sorted[last_10_ind]

  last_5_ind <- ceiling(num_all_features * 0.95)
  last_5_intens <- max_sorted[last_5_ind]

  last_2_ind <- ceiling(num_all_features * 0.98)
  last_2_intens <- max_sorted[last_2_ind]

  num_zeros <- sum(is.na(max_intensity))

  splot_sub <- ggplot()+
    geom_point(aes(x = ceiling(num_all_features/2) : num_all_features, y = subset_intensity), size = 0.5)+
    scale_y_log10() +
    theme_bw()+
    geom_vline(xintercept = last_10_ind, color= "blue", size =0.5, linetype = "dashed")+
    geom_hline(yintercept = last_10_intens, color = "blue", size =0.5, linetype = "dashed")+
    geom_vline(xintercept = last_5_ind, color= "orange", size =0.5, linetype = "dashed")+
    geom_hline(yintercept = last_5_intens, color = "orange", size =0.5, linetype = "dashed")+
    geom_vline(xintercept = last_2_ind, color= "gray50", size =0.5, linetype = "dashed")+
    geom_hline(yintercept = last_2_intens, color = "gray50", size =0.5, linetype = "dashed")+
    annotate(geom="text", x = c(last_2_ind, last_5_ind, last_10_ind), y = rep(subset_intensity[1], 3) ,
              label = c("2%", "5%", "10%"), color = c("gray50", "orange", "blue"))+
    labs(x = "Features without QC and Blank",
         y = "max Intensity",
              subtitle = paste0("Missing Values: " , as.character(num_zeros) , " of " , as.character(length(max_intensity)), " (",
                  as.character(round(num_zeros/length(max_intensity) * 100)), "%)"))+
    theme(axis.text = element_text(size = 11), text = element_text(size = 11))


  # return array and plots as a list
  result <- list()
  result$max_intensity <- max_intensity
  result$plot <- splot
  result$subplot <- splot_sub

  # optional export of plot
  if(length(fn_out) != 0){
    ggsave(filename = fn_out, plot = splot, device = "png", width = 10, height = 7)
  }

  return(result)
}


#' Bar plot with the number of missing values in QC
#'
#' Function counts missing values in the quality control (QC) or other samples and
#' show the results in a bar plot, that can be exported as png.
#'
#' @param df data frame
#' @param list_of_columns List of column names of data frame for counting missing values
#' @param fn_out optional export of the plot as png
#'
#' @return list of array with count of missing values and bar plot
#' @export
#'
#' @examples
count_missing_values <- function(df, list_of_columns, fn_out = NULL){

  df_sub <- df[list_of_columns]

  # count missing values in quality control
  missing_values <- rep(NA, nrow(df_sub))

  for (i in 1:nrow(df_sub)) {
      missing_values[i] <- sum(is.na(df_sub[i, ]))
  }

  # return array and following plot as a list
  results <- list()
  results$missing_values <- missing_values

  # plot
  missing_values <- as.factor(missing_values)

  level_count <- as.numeric(table(missing_values))
  level_percent <- paste0(as.character(round(level_count/ sum(level_count) * 100, digits = 1)), "%")

  viridis_color <- scales::viridis_pal()(length(levels(missing_values)))

  plot <- ggplot()+
    geom_bar(aes(x = missing_values, y=after_stat(count), fill = missing_values), width = 0.8, show.legend = FALSE) +
    coord_flip() +
    labs(y = "Features", x = "Count of Missing Values")+
    ylim(0, max(level_count)*1.1)+
    scale_fill_manual(values = viridis_color)+
    annotate(geom = "text", x = levels(missing_values), y = level_count, label = level_percent, hjust = -0.1,
             color = "black", size = 11/.pt)+
    theme_bw()+
    theme(axis.text = element_text(size = 11))

  # optional export of bar plot
  if(length(fn_out) != 0){
    ggsave(filename = fn_out, plot = plot, device = "png", width = 10, height = 7)
  }


  results$plot <- plot

  return(results)
}


#' Bar plot with relative standard deviation of quality controls
#'
#' Calculates the relative standard deviation (RSD %) of quality controls and
#' creates a bar plot of the ranges 0 - 10 %, 10 - 20 % etc. The bars are
#' grouped be the number of missing values in the quality controls.
#'
#' @param df data frame
#' @param list_of_columns List of column names of data frame to calculate RSD
#' @param fn_out optional export of the plot as png
#'
#' @return list of rsd values, the plot and a subplot without the last bar ( > 100%)
#' @export
#'
#' @examples
calc_RSD <- function(df, list_of_columns, fn_out = NULL){

  df_qc <- df[list_of_columns]

  # calculate RSD of quality controls
  rsd <- rep(NA, nrow(df_qc))

  for(i in 1:nrow(df_qc)){

    temp_row <- as.numeric(df_qc[i, ])

    if(sum(!is.na(temp_row)) > 1){
      mw <- mean(temp_row, na.rm = TRUE)
      rsd[i] <- 100 / mw * sd(temp_row, na.rm = TRUE)
    } else {
        rsd[i] <- NA
    }
  }

  # count (again) missing values in quality controls
  missing_values <- rep(NA, nrow(df_qc))

  for (i in 1:nrow(df_qc)) {
    missing_values[i] <- sum(is.na(df_qc[i, ]))
  }

  df_rsd <- data.frame(QC.RSD = rsd, QC.miss_val = missing_values)
  df_rsd$RSD_bin <- NA

  # create bins for bar plot
  df_rsd$RSD_bin[df_rsd$QC.RSD < 10] <- "< 10"
  df_rsd$RSD_bin[df_rsd$QC.RSD >= 10 & df_rsd$QC.RSD < 20] <- "10 - 20"
  df_rsd$RSD_bin[df_rsd$QC.RSD >= 20 & df_rsd$QC.RSD < 30] <- "20 - 30"
  df_rsd$RSD_bin[df_rsd$QC.RSD >= 30 & df_rsd$QC.RSD < 40] <- "30 - 40"
  df_rsd$RSD_bin[df_rsd$QC.RSD >= 40 & df_rsd$QC.RSD < 50] <- "40 - 50"
  df_rsd$RSD_bin[df_rsd$QC.RSD >= 50 & df_rsd$QC.RSD < 75] <- "50 - 75"
  df_rsd$RSD_bin[df_rsd$QC.RSD >= 75 & df_rsd$QC.RSD < 100] <- "75 - 100"
  df_rsd$RSD_bin[df_rsd$QC.RSD >= 100 ] <- "> 100"
  df_rsd$RSD_bin[is.na(df_rsd$QC.RSD)] <- "no value"


  df_rsd$RSD_bin <- factor(df_rsd$RSD_bin,
                             levels = c("< 10", "10 - 20", "20 - 30", "30 - 40",
                                        "40 - 50", "50 - 75", "75 - 100", "> 100", "no value"))

  df_rsd$QC.miss_val <- as.factor(df_rsd$QC.miss_val)
  viridis_color <- scales::viridis_pal()(length(levels(df_rsd$QC.miss_val)))

  # count percentages of features for each bin
  level_count <- as.numeric(table(df_rsd$RSD_bin))
  level_percent <- paste0(as.character(round(level_count/ sum(level_count) * 100, digits = 1)), "%")

  # bar plot
  plot <- ggplot(df_rsd, aes(x= RSD_bin, fill = QC.miss_val))+
          geom_bar(width = 0.8)+
          ylim(0, max(level_count)*1.1)+
          scale_fill_manual(values = viridis_color)+
          labs(fill = "missing values in QC", x = "RSD %", y = "Counts")+
          annotate(geom = "text", x = levels(df_rsd$RSD_bin), y = level_count, label = level_percent, vjust = -0.5,
                   color = "black", size = 11/.pt)+
          theme_bw()+
          theme(axis.text = element_text(size = 11))


  subplot <-  ggplot(subset(df_rsd, RSD_bin != "> 100" & RSD_bin != "no value" ),aes(x= RSD_bin, fill = QC.miss_val))+
              geom_bar(width = 0.8, show.legend = FALSE)+
              scale_fill_manual(values = viridis_color)+
              labs(x = "RSD %", y = "Counts")+
              theme_bw()+
              theme(axis.text = element_text(size = 11))

  # return array and plots as a list
  result <- list()
  result$RSD <- rsd
  result$plot <- plot
  result$subplot <- subplot

  # optional export of bar plot
  if(length(fn_out) != 0){
    ggsave(filename = fn_out, plot = plot, device = "png", width = 10, height = 7)
  }

  return(result)
}









