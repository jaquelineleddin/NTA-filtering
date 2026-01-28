#' Converts lipid names
#'
#' If lipid name is of type 'class a:b_c:d' or 'class a:b;O/c:d',
#' function will sum up the numbers to 'class a+c:b+d' or 'class a+c:b+d;O'
#'
#' @param lipid_name name of lipid as string
#' @noRd
#'
#' @return corrected lipid name
sum_lipid_names <- function(lipid_name){

  # check class type
  # LP x:y
  if(grepl('^[A-Z, a-z]{1,6} [0-9]', lipid_name)){
    split_class <- unlist(strsplit(lipid_name, split = " ", fixed = TRUE))
    class_name <- paste0(split_class[1], ' ')
  }
  # LP O-x:y
  else if(grepl('^[A-Z, a-z]{1,3} [O,P]-', lipid_name)){
    split_class <- unlist(strsplit(lipid_name, split = "-", fixed = TRUE))
    class_name <- paste0(split_class[1], '-')
  }
  else{
    class_name <- 'None'
  }

  # check for standard
  if(class_name == 'None'){
    return(lipid_name)
  }else if(grepl('[(]d[0-9][)]', split_class[2])){
    standard_split <- unlist(strsplit(split_class[2], split = '(', fixed = TRUE))
    standard <- paste0('(', standard_split[2])
    number <- standard_split[1]
  } else{
    number <- split_class[2]
    standard <- NULL
  }

  # sum numbers
  # case: /
  if(grepl('O/', number)){
    number_split <- unlist(strsplit(number, split = '/', fixed = TRUE))

    o_split <- unlist(strsplit(number_split[1], split = ";", fixed = TRUE))

    vec1 <- as.numeric(unlist(strsplit(o_split[1], split = ':', fixed = TRUE)))
    vec2 <- as.numeric(unlist(strsplit(number_split[2], split = ':', fixed = TRUE)))

    sum <- as.character(vec1 + vec2)

    number <- paste0(sum[1], ":", sum[2], ';', o_split[2])
  }
  # case: _
  else if(grepl('_', number)){

    # case: _ and ;nO
    if(grepl(';[0-4]O', number)){
      split_o <- unlist(strsplit(number, split = ';', fixed = TRUE))
      number <- split_o[1]
      # now standard means the last position
      standard <- paste0(';', split_o[2])
    }

    under_split <- unlist(strsplit(number, split = '_', fixed = TRUE))

    sum <- c(0,0)

    for(str in under_split){
      sum <- sum +  as.numeric(unlist(strsplit(str, split = ':', fixed = TRUE)))
    }
    sum <- as.character(sum)
    number <- paste0(sum[1], ':', sum[2])
  }
  # case: x:y/0:0
  else if(grepl('/0:0', number)){
    number <- sub('/0:0', '', number)
  } else{
    return(lipid_name)
  }

  new_name <- paste0(class_name, number, standard)

  return(new_name)
}


#' Find index of lipid names, that must be corrected
#'
#' If lipid name is not like 'class x:y' or 'class x:y;O',
#' function will return index of lipid name in the array
#'
#' @param lipid_names_array array of lipid names
#' @noRd
#'
#' @return array with the indices
#'
find_lipid_index <- function(lipid_names_array){

  index <- c()

  for(i in 1:length(lipid_names_array)){
    if(!(
      grepl('^[A-Z, a-z]{1,6} [0-9]{1,2}:[0-9]{1,2}$', lipid_names_array[i]) |
      grepl('^[A-Z, a-z]{1,3} [0-9]{1,2}:[0-9]{1,2}[(]d[0-9][)]$', lipid_names_array[i]) |
      grepl('^[A-Z, a-z]{1,3} [O,P]-[0-9]{1,2}:[0-9]{1,2}$', lipid_names_array[i]) |
      grepl('^[A-Z, a-z]{1,6} [0-9]{1,2}:[0-9]{1,2};[0-9]O$', lipid_names_array[i])
    )){
      index <- c(index, i)
    }
  }
  return(index)
}

#' Correction of lipid names in data frame
#'
#' If the metabolites in a data frame are lipids, the names
#' will be corrected and insert as new column next to
#' metabolite.name. w/o MS2: and everything after | will be removed.
#' Summation of number of carbon atoms and double bonds is optional.
#'
#' @param lipid_df data frame with names in column 'Metabolite.name'
#' @param with_summation summation of number of carbon atoms and double bonds, e.g. x:y_a:b to x+a:y+b.
#' Default is False.
#'
#' @return data frame with new name column inserted
#' @export
lipid_names_correction <- function(lipid_df, with_summation = FALSE){
  # insert new name column
  names_corrected <- lipid_df$Metabolite.name

  # remove w/o
  names_corrected <- sub("w/o MS2:", "", names_corrected)

  # remove |...
  for(i in 1:length(names_corrected)){
    names_corrected[i] <- unlist(strsplit(names_corrected[i], split = "|", fixed = TRUE))[1]
  }

  if(with_summation){
    # find index for summation
    index <- find_lipid_index(names_corrected)

    # summation
    for(i in index){
      names_corrected[i] <- sum_lipid_names(names_corrected[i])
    }
  }

  # insert in df
  col_num <- which(colnames(lipid_df) == 'Metabolite.name')
  df <- cbind(lipid_df[1:col_num],
              data.frame("corrected names" = names_corrected), lipid_df[ (col_num+1):ncol(lipid_df)])

  return(df)
}
