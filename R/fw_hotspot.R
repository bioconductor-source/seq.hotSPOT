#' forward amplicon ranking
#'
#' @description
#' create a targeted sequencing panel by finding which amplicons
#' will likely capture the most mutations
#'
#' @details
#' Forward Selection Sequencing Panel Identifier
#'
#' 1.	Amplicons covering hotspots less than or equal to one amplicon in length,
#' are added to the final sequencing panel dataset.
#'
#' 2.	For amplicons covering larger hotspot regions, the algorithm uses a
#' forward selection method to determine the optimal combination of amplicons
#' to use in the sequencing panel:
#'   a.	the algorithm first identifies the amplicon containing the highest
#'   number of mutations
#'   b.	the algorithm then identifies the next amplicon, which contains the
#'   highest number of new mutations.
#'   c.	this process continues until all mutations are covered by at least
#'   one amplicon
#'
#' 3.	Each of these amplicons are then added to the final sequencing panel,
#' with their own unique IDs.
#'
#' 4.	All amplicons in the final sequencing panel are ranked from highest to
#' lowest based on the number of mutations they cover.
#'
#' 5.	The algorithm then calculates the cumulative base-pair length and the
#' cumulative mutations covered by each amplicon.
#'
#' 6.	Dependent on the desired length of the targeted panel, a cutoff may be
#' applied to remove all amplicons which fall below a set cumulative length.

#'
#' @importFrom hash make.keys
#' @importFrom hash hash
#' @importFrom hash values
#' @importFrom stats ftable
#' @importFrom stats median
#' @param bins A dataframe containing all potential amplicons
#' @param data A dataframe containing the location of each mutation.
#' @param amp The length of amplicons in number of base pairs
#' @param len The total length of sequencing panel in number of base pairs
#' @param include_genes True or False based on whether dataset includes gene names
#'
#' @return A dataframe containing the genomic coordinates for targeted sequencing panel
#'
#'@examples
#'
#' data("mutation_data")
#' my_bins <- amp_pool(data, 100)
#'
#' fw_hotspot(my_bins, data, 100, 1000, TRUE)
#'
#' @export
#'


fw_hotspot <- function(bins, data, amp, len, include_genes){
  if (base::isFALSE(is.data.frame(data))){stop("Input data must be in the form of data.frame")}
  if (base::isFALSE("chr" %in% colnames(data))){stop("input data must contain column 'chr'")}
  if (base::isFALSE("pos" %in% colnames(data))){stop("input data must contain column 'pos'")}
  if (base::isFALSE(is.numeric(data$pos))){stop("column 'pos' must be in the form of numeric")}
  if (base::isFALSE(is.numeric(data$chr)) & base::isFALSE(is.character(data$chr)))
  {stop("column 'chr' must be in the form of either numeric or character")}
  if (base::isTRUE(include_genes) & base::isFALSE("gene" %in% colnames(data)))
  {stop("input data must contain column 'gene'")}
  if (base::isTRUE(include_genes) & base::isFALSE(is.character(data$gene))){
    stop("column 'gene' must be in the form of character")
  }
  if (amp < 1){stop("amplicon length must be greater than 1 base pair")}
  if (base::isFALSE(is.data.frame(bins))){stop("Bins must be in the form of data.frame")}
  if (len < amp){stop("Panel length must be equal to or greater than one amplicon")}
  if (len < nrow(data)){stop("Panel length cannot be greater than number of mutations in input data")}
  bins <- bins[,-c(6:7)]
  pos <- data$pos
  pos_freq <- data.frame(ftable(pos)) #make frequency table for each position
  data <- merge(data, pos_freq, by = "pos") #merge frequency table with original df
  data <- unique(data) #keep only unique positions
  panel_length <- len #set panel and amplicon length
  amp_length <- amp
  all_pos_bins <- bins
  ordered_bin <- all_pos_bins[order(-all_pos_bins$count), ]
  big_spots <- ordered_bin[!ordered_bin$id == "x", ]
  unique_spots <- unique(as.list(big_spots$id))      # we then need to find the optimal combination of amplicons to cover each large hotspot
  keys <- make.keys(data$pos)
  dict <- hash(keys = keys, values = data$Freq)      # create a new dictionary of mutation positions and frequency of mutations
  possible_bins <- data.frame()
  for (spot in unique_spots){
    spot_sub <- subset(big_spots, id == spot)
    big_bin_all <- list()
    while (nrow(spot_sub) > 0){
      row <- data.frame()                            # for each hotspot, subset a dataframe with only amplicons which cover that hotspot
      new_sub <- data.frame()
      big_spot <- spot_sub[1,]
      big_bin <- big_spot$upperbound:big_spot$lowerbound
      big_bin_all <- c(big_bin_all, big_bin)
      possible_bins <- rbind(possible_bins, big_spot)   # first find the amplicon which contains the most mutations and add to final list
      if (nrow(spot_sub) > 1){
        for (k in seq(2,nrow(spot_sub))){
          row <- data.frame()
          test_bin <- spot_sub$upperbound[k]:spot_sub$lowerbound[k]
          overlap <- intersect(big_bin_all, test_bin)    # then for each remaining bin, adjust mutation count so they are only showing unique mutations
          if (length(overlap > 0)){                      # continue process until there are no amplicons left which would capture additional unique mutations
            unique_bin <- setdiff(test_bin, overlap)
            mutations <- intersect(unique_bin, pos)
            if (length(mutations) > 0){
              mutations_keys <- make.keys(mutations)
              new_count <- sum(as.vector(values(dict, keys = mutations_keys)))
              row <- data.frame("lowerbound" = min(test_bin), "upperbound" = max(test_bin),
                                "chromosome" = spot_sub$chromosome[1], "count" = new_count,"id" = spot_sub$id[1])
            }
          }
          if (length(overlap) <= 0) {row <- spot_sub[k,]}
          new_sub <- rbind(new_sub, row)
        }
      }
      spot_sub <- new_sub
      if (nrow(spot_sub > 0)){spot_sub <- spot_sub[order(-spot_sub$count), ]}
    }
    gc()
  }
  final_bin <- rbind(subset(ordered_bin, id == "x"), possible_bins)
  final_bin <- final_bin[order(-final_bin$count), ]
  cutoff_point <- final_bin$count[[round(panel_length / amp_length)]]
  subset_final_bin <- subset(final_bin, count == cutoff_point)
  dif_list <- list()
  for (i in seq_len(nrow(subset_final_bin))){ ##
    whole_region <- subset_final_bin$lowerbound[[i]]:subset_final_bin$upperbound[[i]]
    bin_mutations <- intersect(pos, whole_region)
    mid_diff <- abs((max(whole_region) - max(bin_mutations)) - (min(bin_mutations) - min(whole_region)))
    dif_list <- c(dif_list, mid_diff)
  }
  subset_final_bin$dif <- unlist(dif_list)
  subset_final_bin <- subset_final_bin[order(-subset_final_bin$dif), ]
  subset_final_bin <- subset_final_bin[,-ncol(subset_final_bin)]
  df_range_min <- min(which(final_bin$count == cutoff_point))
  df_range_max <- max(which(final_bin$count == cutoff_point))
  final_bin[df_range_min:df_range_max,] <- subset_final_bin[seq_len(nrow(subset_final_bin)),]
  bin_len_list <- list()
  for (y in seq_len(nrow(final_bin))) {      #calculating bin length using upper and lower bound positions
    upper <- final_bin$upperbound[[y]]
    lower <- final_bin$lowerbound[[y]]
    bin_len <- upper - lower + 1# add 1 to include both upper and lowerbound plus all bp inbetween
    bin_len_list <- c(bin_len_list, bin_len) #make list of all bin lengths
  }
  final_bin$bin_length <- unlist(bin_len_list) # add to dataframe
  bin_len_list <- as.numeric(final_bin$bin_length) #new lists since dataframe has been reordered
  count_list <- as.numeric(final_bin$count)
  bin_len_cum <- list()
  cum_mut_list <- list()
  cum_mut_list[[1]] <- count_list[[1]]  #start list with first bin length and first mutation count
  bin_len_cum[[1]] <- bin_len_list[[1]]
  for (u in seq(2, length(count_list))) {     #repeat for all other bins
    cum_mut_list[[u]] <- cum_mut_list[[u-1]] + count_list[[u]]  #add all calculated bin lengths and mutation counts to value one row above
    bin_len_cum[[u]] <- bin_len_cum[[u-1]] + bin_len_list[[u]]
  }
  final_bin$Cummulative_Bin_Length <- unlist(bin_len_cum)
  final_bin$Cummulative_Mutations <- unlist(cum_mut_list)
  final_bin <- final_bin[final_bin$Cummulative_Bin_Length <= panel_length,] ##apply bp length cutoff if needed
  final_bin <- final_bin[,c(1:4,7:8)]
  colnames(final_bin) <- c("Lowerbound", "Upperbound", "Chromosome", "Mutation Count",
                           "Cumulative Panel Length", "Cumulative Mutations")
  if (include_genes == TRUE){
  gene_list <- c()
  for (i in seq_len(nrow(final_bin))){
    data.sub <- subset(data, chr == final_bin$Chromosome[[i]] & pos %in%
                         final_bin$Lowerbound[[i]]:final_bin$Upperbound[[i]])
    gene_list <- c(gene_list, paste(unique(data.sub$gene), collapse = ""))
  }
  final_bin$Gene <- unlist(gene_list)
  }

  return(final_bin)
}
