#' amplicon finder
#'
#' @description
#' create a dataframe containing the coordinates of all potential
#' amplicons for single chromosome
#'
#'
#' @importFrom hash make.keys
#' @importFrom hash hash
#' @importFrom hash values
#' @importFrom stats ftable
#' @importFrom stats median
#' @param data A dataframe containing the location of each mutation.
#' @param chr The chromosome of interest
#' @param amp_len The length of amplicons in number of base pairs
#'
#' @keywords internal
#'
#' @return A dataframe containing the genomic coordinates of all potential amplicons on a single chromosome
#'
#' @examples
#'
#' data("mutation_data")
#' amplicon_finder(data, 1, 100)
#'
#' @noRd
#'

amplicon_finder <- function(data, chr, amp_len){
  pos <- data$pos
  keys <- make.keys(data$pos)
  dict <- hash(keys = keys, values = data$Freq) #create dictionary for mutation locations and mutation frequency
  make_bins <- sort(unique(pos))
  bins <- data.frame()
  diffs <- diff(make_bins)
  large_diffs <- c(0, which(diffs > amp_len))
  for (d in seq(1,(length(large_diffs)-1))){
    row <- data.frame("lowerbound" = make_bins[[(large_diffs[[d]] + 1)]], "upperbound" = make_bins[[large_diffs[[d + 1]]]], "chromosome" = chr)
    bins <- rbind(bins, row)
    gc()
  }
  mut_locations <- data.frame()
  possible_bins <- data.frame()
  for (i in seq_len(nrow(bins))){
    up_mut <- bins$upperbound[[i]]
    low_mut <- bins$lowerbound[[i]]
    vec <- up_mut:low_mut
    dif <- length(vec) - amp_len
    if (dif <= 0){
      up <- bins$upperbound[[i]]
      low <- bins$lowerbound[[i]]
      vec <- up:low
      mutations <- intersect(vec, pos) #find all mutations in total bin region
      mutations_keys <- make.keys(mutations)
      count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
      mut_all <- c()
      for (m in seq_along(mutations)){
        mut_all <- c(mut_all, rep(mutations[[m]], as.vector(values(dict, keys = mutations_keys))[[m]]))
      }
      weighted_mid <- round(mean(mut_all))
      weighted_mid_max <- (ceiling(mean(vec)) + (amp_len/2))
      weighted_mid_min <- (floor(mean(vec)) - (amp_len/2))
      if (weighted_mid < weighted_mid_min){final_mid <- weighted_mid_min}
      if (weighted_mid > weighted_mid_max){final_mid <- weighted_mid_max}
      if (weighted_mid > weighted_mid_min & weighted_mid < weighted_mid_max){final_mid <- weighted_mid}
      up <- final_mid + ceiling((amp_len - 1) / 2)
      low <- final_mid - floor((amp_len - 1) / 2)
      vec <- up:low
      mutations <- intersect(vec, pos) #find all mutations in total bin region
      mutations_keys <- make.keys(mutations)
      count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
      row <- data.frame("lowerbound" = low, "upperbound" = up,
                        "chromosome" = bins$chromosome[[i]], "count" = count,"id" = "x", "mut_lowerbound" = low_mut, "mut_upperbound" = up_mut)
      possible_bins <- rbind(possible_bins, row)
    }
    else{
      up_mut <- bins$upperbound[[i]]
      low_mut <- bins$lowerbound[[i]]
      vec <- up_mut:low_mut
      mut_list <- intersect(vec, pos)
      id = paste(as.character(bins$chromosome[[i]]), as.character(i), sep = "-")
      b_list <- list()
      for (mut in mut_list){
        bin1 <- (mut - amp_len):(mut - 1)
        b_list[[1]] <- bin1
        bin2 <- (mut - (amp_len - 1)):(mut)
        b_list[[2]] <- bin2
        bin3 <- (mut + amp_len):(mut + 1)
        b_list[[3]] <- bin3
        bin4 <- (mut + (amp_len - 1)):(mut)
        b_list[[4]] <- bin4
        bin5 <- (mut - ceiling(amp_len / 2)):(mut - floor(amp_len / 2))
        b_list[[5]] <- bin5

        for (b in b_list){
          low <- min(b)
          up <- max(b)
          mutations <- intersect(b, pos)
          if (length(mutations > 0)){
            mutations_keys <- make.keys(mutations)
            count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
            row <- data.frame("lowerbound" = low, "upperbound" = up,
                              "chromosome" = bins$chromosome[[i]], "count" = count,"id" = id, "mut_lowerbound" = low_mut, "mut_upperbound" = up_mut)
            possible_bins <- rbind(possible_bins, row)
          }
        }
      }

    }
  }
  return(possible_bins)       # all amplicons are added to dataframe
}

#' create amplicon pool
#'
#' @description
#' create a dataframe containing the coordinates of all potential
#' amplicons for hotspot testing
#'
#' @details
#' This algorithm searches the mutational dataset (input) for mutational
#' hotspot regions on each chromosome:
#'
#' 1.	Starting at the mutation with the lowest chromosomal position
#' (primary mutation), using a modified rank and recovery system,
#' the algorithm searches for the closest neighboring mutation.
#'
#' 2.	If the neighboring mutation is less than one amplicon, in distance,
#' away from the primary mutation, the neighboring mutation is included
#' within the hotspot region.
#'    a.	This rank and recovery system is repeated, integrating mutations
#'    into the hotspot region until the neighboring mutation is greater
#'    than or equal to the length of one amplicon in distance,
#'    from the primary mutation.
#'    b.	Once neighboring mutations equal or exceed one amplicon in distance
#'    from the primary mutation, incorporation into the hotspot region,
#'    halts incorporation.
#'
#' 3.	For hotspots within the one amplicon range, from the lowest to highest
#' mutation location, this area is covered by a single amplicon and added to
#' an amplicon pool, with a unique ID.
#'    a.	The center of these single amplicons is then defined by the weighted
#'    distribution of mutations.
#'
#' 4.	For all hotspots larger than one amplicon, the algorithm examines
#' 5 potential amplicons at each covered mutation in the hotspot:
#'    a.	one amplicon directly upstream of the primary mutation
#'    b.	one amplicon directly downstream of the primary mutation
#'    c.	one amplicon including the mutation at the end of the read and
#'    base pairs (amplicon length - 1) upstream
#'    d.	one amplicon including the mutation at the beginning of the read and
#'    base pairs (amplicon length - 1) downstream
#'    e.	one amplicon with the mutation directly in the center.
#'
#' 5.	All amplicons generated for each hotspot region of interest, are assigned a
#' unique ID and added to the amplicon pool.
#'
#' The mutation dataset should include two columns containing the chromosome and
#' genomic position, the columns should be names "chr" and "pos" respectively.
#' Optionally the gene names for each mutation may be included under a
#' column names "gene".
#'
#' @importFrom hash make.keys
#' @importFrom hash hash
#' @importFrom hash values
#' @importFrom stats ftable
#' @importFrom stats median
#' @param data A dataframe containing the location of each mutation.
#' @param amp The length of amplicons in number of base pairs
#'
#' @return A dataframe containing the genomic coordinates of all potential amplicons
#'
#' @examples
#'
#' data("mutation_data")
#' amp_pool(data, 100)
#'
#'
#' @export
#'



amp_pool <- function(data, amp){
  if (base::isFALSE(is.data.frame(data))){stop("Input data must be in the form of data.frame")}
  if (base::isFALSE("chr" %in% colnames(data))){stop("input data must contain column 'chr'")}
  if (base::isFALSE("pos" %in% colnames(data))){stop("input data must contain column 'pos'")}
  if (base::isFALSE(is.numeric(data$pos))){stop("column 'pos' must be in the form of numeric")}
  if (base::isFALSE(is.numeric(data$chr)) & base::isFALSE(is.character(data$chr)))
  {stop("column 'chr' must be in the form of either numeric or character")}
  if (amp < 1){stop("amplicon length must be greater than 1 base pair")}
  pos <- data$pos
  pos_freq <- data.frame(ftable(pos)) #make frequency table for each position
  data <- merge(data, pos_freq, by = "pos") #merge frequency table with original df
  data <- unique(data) #keep only unique positions
  chrom_list <- unique(as.list(data$chr)) # list of all chromosomes found in this dataset
  all_pos_bins <- data.frame()
  for (chrom in chrom_list){
    chromosome <- subset(data, chr == chrom)
    chromsort <- chromosome[order(chromosome$pos), ] #order df from lowest to highest position
    chrom_bins <- amplicon_finder(chromsort, chrom, amp) #run through above function
    all_pos_bins <- rbind(all_pos_bins, chrom_bins) #add all bins together from all chromosomes
    rm(list = ls()[! ls() %in% c("amplicon_finder", "data", "all_pos_bins", "chrom_list", "amp")])
    gc()
  }
  return(all_pos_bins)
}
