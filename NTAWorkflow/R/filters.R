filter_rows <- function(df, colname, threshold, mode = "<="){

  if(mode == "<="){
    df_subset <- subset(df, df[[colname]] <= threshold)
  }
  else if(mode == ">="){
    df_subset <- subset(df, df[[colname]] >= threshold)
  }
  else{
    stop("choose for mode <= or >=")
  }

  return(df_subset)
}


#' Filter out rows in data frame
#'
#' The parameter colnames, thresholds and modes define the drop of rows by condition. E.g. the first
#' entry of the default values stands for the condition df$max.Intensity >= 10^5, so all rows with smaller
#' values than the threshold will be filtered out.
#'
#' @param df data frame
#' @param colnames array of column names for the filtering process
#' @param thresholds thresholds for every column define in colnames to filter out rows of the data frame
#' @param modes choose >= or <=
#'
#' @return data frame
#' @export
#'
#' @examples
filter_dataframe <- function(df, colnames = c("max.Intensity", "QC.missVal", "QC.RSD"),
                             thresholds = c(10^5, 0, 35), modes = c(">=", "<=", "<=")){

  df_subset <- df

  for(i in 1:length(colnames)){
    if(colnames[i] %in% colnames(df_subset)){
      df_subset <- filter_rows(df_subset, colnames[i], thresholds[i], modes[i])
    } else{
        stop(colnames[i], " is not a column name of the input data frame!")
    }

  }
 return(df_subset)
}


#' Logarithmic plot of maximum intensities
#'
#' @param intensity_array Array of maximum intensities calculated by the function calc_max_intensity
#'
#' @return Plot
#' @export
#'
#' @examples
plot_max_intensity <- function(intensity_array){

  # sort intensity array
  max_sorted <- sort(intensity_array, decreasing = TRUE, na.last = NA)

  # calculate major ticks for logarithmic plot
  breaks <- breaks_for_splot(max(max_sorted))

  # plot
  splot <- ggplot()+
           geom_point(aes(x = 1:length(max_sorted), y = max_sorted), size = 0.5)+
           scale_y_log10( breaks = breaks) +
           theme_bw()+
           labs(x = "Features of Samples",
                y = "Maximum of Intensities")+
           theme(axis.text = element_text(size = 11))

  return(splot)
}


#' Bar plot of missing values
#'
#' @param missing_values_array Array with missing values counted by the function count_missing_values
#'
#' @return bar plot
#' @export
#'
#' @examples
plot_missing_values <- function(missing_values_array){

  missing_values <- as.factor(missing_values_array)

  viridis_color <- scales::viridis_pal()(length(levels(missing_values)))

  barplot <- ggplot()+
             geom_bar(aes(x = missing_values, y=after_stat(count), fill = missing_values), width = 0.8, show.legend = FALSE) +
             coord_flip() +
             labs(y = "Features", x = "Count of Missing Values")+
             scale_fill_manual(values = viridis_color)+
             theme_bw()+
             theme(axis.text = element_text(size = 11))

  return(barplot)
}


#' Bar plot of relative standard deviation (rsd)
#'
#' @param rsd_array Array of rsd, that was calculated by the function calc_RSD
#'
#' @return bar plot
#' @export
#'
#' @examples
plot_rsd <- function(rsd_array){

  df_rsd <- data.frame(QC.RSD = rsd_array)
  df_rsd$RSD_bin <- NA

  # create bins for rsd
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

  level_count <- as.numeric(table(df_rsd$RSD_bin))

  barplot <- ggplot(df_rsd, aes(x= RSD_bin))+
             geom_bar(width = 0.8)+
             ylim(0, max(level_count)*1.1)+
             labs(fill = "missing values in QC", x = "RSD %", y = "Counts")+
             theme_bw()+
             theme(axis.text = element_text(size = 11))

  return(barplot)
}


#' Filter out rows of a data frame
#'
#' @param df data frame
#' @param list_of_colnames columns must be logical, thus contain False or True,
#' e.g. MS.MS.assigned, RT.matched, m.z.matched or MS.MS.matched
#'
#' @return data frame
#' @export
#'
#' @examples
filter_matched <- function(df, list_of_colnames){

  subdf <- df

  for(colname in list_of_colnames){
    subdf <- subdf[ which(subdf[ , colname]) , ]
  }

  return(subdf)

}


#' Adds column to a data frame for feature names
#'
#' Combines columns of the data frame for new feature names. Each item of the columns is
#' separated with _ . The new column is add at the second place of the data frame.
#'
#' @param df data frame with columns of MSDial
#' @param colnumbers_to_use Column number of the columns to combine.
#' E.g. Alignment.ID = 1, Average.Rt.min = 2, Average.Mz = 3, Metabolite.name = 4.
#' @param remove_wo If True and Metabolite.name is used,
#' then "w/o MS2:" will be removed in the combined name.
#'
#' @return data frame
#' @export
#'
#' @examples
add_feature_names <- function(df, colnumbers_to_use, remove_wo = TRUE){

  # convert first column to string
  feature_names <- as.character(df[ ,colnumbers_to_use[1]])

  # concatenate all given columns
  for(i in 2:length(colnumbers_to_use)){

    feature_names <- paste(feature_names, as.character(df[ ,colnumbers_to_use[i]]), sep = "_")
  }

  # insert new name column at the second place
  df_new <- cbind(df[1], data.frame(name = feature_names, stringsAsFactors = FALSE), df[2:ncol(df)])

  # remove w/o MS2:
  if(remove_wo){
    df_new$name <- sub("w/o MS2:", "", df_new$name)
  }

  return(df_new)

}


#' Exports data frame for further processing with MetaboAnalyst or SIMCA
#'
#' Function creates new name for the features, takes the columns of samples, transposes the data frame and adds
#' a column with class names to group the samples. At least the data frame can be saved as a text file.
#'
#' @param df data frame
#' @param colnum_of_name Number of the column in the data frame to use for feature names.
#' @param all_sample_cols list of column names of all samples
#' @param sort_out List of target words to sort out quality control (e.g. target word is "QC") or blanks
#' @param classes list of strings for grouping the samples
#' @param remove_X RStudio put an X before the name of a column, if it begins with a number.
#' After transpose of the matrix, X will be remove.
#'
#' @return data frame
#' @export
#'
#' @examples
transpose_df <- function(df, colnum_of_name, all_sample_cols, sort_out = NULL,
                           classes = NULL, remove_X = TRUE){

  # get sample columns
  df_export <- df[all_sample_cols]

  # sort out QC or Blanks
  if(!is.null(sort_out)){

    for (name in sort_out) {
      cols_out <- grep(name, colnames(df_export))

      if(length(cols_out) == 0){
        stop(name, " is in not part of a column name, so I can't sort out. Please check colnames in dataframe."  )
      }

    df_export <- subset(df_export, select = -cols_out)
    }
  }


  # transpose
  df_export <- as.data.frame(t(df_export))
  colnames(df_export) <- df[ , colnum_of_name]

  # add classes
  if(!is.null(classes)){

    df_export <- cbind(data.frame(class = rep(NA, nrow(df_export))), df_export)

    for(name in classes){
      rownum <- grep(name, row.names(df_export))

      if(length(rownum) == 0){
        stop(name, " is not part of a row name, so I can't associate a class. Please check row.names in dataframe."  )
      }

      df_export[rownum, "class"] <- name
    }
  }


  df_export <- cbind(data.frame("Name" = row.names(df_export)), df_export)
  row.names(df_export) <- 1:nrow(df_export)

  # is column name start with a number, R set an X before it. We will remove it now.
  if(remove_X){
      df_export$Name <- gsub("X", "", df_export$Name)
  }

  return(df_export)
}


#' Normalize values of sample columns
#'
#' Function divides values of sample columns, e.g. to normalize them by weight of samples.
#' For that, a data frame contains different target words for each sample and the factor for
#' division.
#'
#' @param df data frame, that contains columns with samples
#' @param df_norm_factors data frame, that contains in first column part of sample names, to
#' address the columns in data frame df. The second column contains the factors for division.
#' @param decimals numbers of decimals for rounding
#'
#' @return data frame
#' @export
#'
#' @examples
normalize_intensity <- function(df, df_norm_factors, decimals = 6){  # first col is target word in sample, second col is factor for division

  df_normalize <- df

  for(i in 1:nrow(df_norm_factors)){
    # go through dataframe of weights and find columns in the dataframe of samples
    col_index <- grep(df_norm_factors[i, 1], colnames(df_normalize))

    # columns must be unambiguous: every name in the dataframe of weights must reference one sample
    if(length(col_index) != 1){
      stop(df_norm_factors[i, 1], " is no part of sample names or belongs to more than one sample!")
    }

    df_normalize[ , col_index] <- round(df_normalize[ , col_index] / df_norm_factors[i, 2], digits = decimals)

  }
  return(df_normalize)
}
