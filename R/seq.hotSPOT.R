

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
#' @import hash
#' @import rlist
#' @import R.utils
#' @param data A dataframe containing the location of each mutation.
#' @param amp The length of amplicons in number of base pairs
#' @param len The total length of sequencing panel in number of base pairs
#'
#' @return A dataframe containing the genomic coordinates of all potential amplicons
#'
#' @export
#'
#' @examples
#'
#' load(example.RData)
#' amp_pool(data, 100, 5000)
#'



amp_pool <- function(data, amp, len){
    panel_length = len #set panel and amplicon length
    amp_length = amp
    pos <- data$pos
    pos_freq <- data.frame(ftable(pos)) #make frequency table for each position
    data <- merge(data, pos_freq, by = "pos") #merge frequency table with original df
    data <- unique(data) #keep only unique positions
    amplicon_finder <- function(data, chr){
        pos <- data$pos
        keys <- make.keys(data$pos)
        dict <- hash(keys = keys, values = data$Freq) #create dictionary for mutation locations and mutation frequency
        make_bins <- sort(unique(pos))
        bins <- data.frame()
        while (length(make_bins) > 1){   ## groups all mutations by relative position
            temp_bins <- make_bins         ## starts with first mutation, checks for next closes mutation
            i <- 1                         ## mutations are added to the group if they are less than one amplicon length away from neighboring mutation
            start <- make_bins[1]          ## if mutations are further than one amplicon away, new group is started
            p_dif <- 1
            go = "yes"
            while (go == "yes"){
                p1 <- temp_bins[i]
                p_next <- temp_bins[i + 1]
                i <- i + 1
                if (i <= length(make_bins)){
                    p_dif <- p_next - p1
                }
                if (i > length(make_bins)){
                    p_dif = amp_length + 100
                    go = "no"
                }
                if (p_dif > amp_length){
                    p_next <- p1
                }
            }
            end <- p_next
            vec <- as.list(start:end)
            keep <- setdiff(make_bins, vec)
            make_bins <- keep
            chr = chr
            row <- data.frame("lowerbound" = start, "upperbound" = end, "chromosome" = chr)
            bins <- rbind(bins, row)
        }
        if (length(make_bins) == 1){
            row <- data.frame("lowerbound" = make_bins, "upperbound" = make_bins, "chromosome" = chr)  # all groups are added to dataframe of bins
            bins <- rbind(bins, row)
        }
        mut_locations <- data.frame()
        possible_bins <- data.frame()
        for (i in seq_len(nrow(bins))){             ## for each bin, finds amplicons which may optimally cover that region
            up_mut <- bins$upperbound[[i]]
            low_mut <- bins$lowerbound[[i]]   # find length of bin
            vec <- up_mut:low_mut
            dif <- length(vec) - amp_length
            if (dif <= 0){                    # if the bin is <= amplicon length, we only need to generate a single amplicon
                up <- bins$upperbound[[i]]
                low <- bins$lowerbound[[i]]
                vec <- up:low
                mutations <- intersect(vec, pos)    #find all mutations in total bin region
                mutations_keys <- make.keys(mutations)
                count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                mut_all <- c()
                for (m in seq_len(length(mutations))){
                    mut_all <- c(mut_all, rep(mutations[[m]], as.vector(values(dict, keys = mutations_keys))[[m]]))
                }
                weighted_mid <- round(mean(mut_all))
                weighted_mid_max <- (ceiling(mean(vec)) + (amp_length/2))        # the single amplicon we choose to cover this bin is weighted based on mutation distribution
                weighted_mid_min <- (floor(mean(vec)) - (amp_length/2))
                if (weighted_mid < weighted_mid_min){final_mid <- weighted_mid_min}
                if (weighted_mid > weighted_mid_max){final_mid <- weighted_mid_max}
                if (weighted_mid > weighted_mid_min & weighted_mid < weighted_mid_max){final_mid <- weighted_mid}
                up <- final_mid + ceiling((amp_length - 1) / 2)
                low <- final_mid - floor((amp_length - 1) / 2)
                vec <- up:low
                mutations <- intersect(vec, pos) #find all mutations in total bin region
                mutations_keys <- make.keys(mutations)
                count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                row <- data.frame("lowerbound" = low, "upperbound" = up,
                              "chromosome" = bins$chromosome[[i]], "count" = count,"id" = "x")
                possible_bins <- rbind(possible_bins, row)
            }
            else{                                    # for bins greater than one amplicon length, we generate 5 possible amplicons for each mutation within the bin
                up_mut <- bins$upperbound[[i]]
                low_mut <- bins$lowerbound[[i]]
                vec <- up_mut:low_mut
                mut_list <- intersect(vec, pos)
                id = paste(as.character(bins$chromosome[[i]]), as.character(i), sep = "-") # every time a new bin is created, it is assigned a unique id
                b_list <- list()
                for (mut in mut_list){                      # creating 5 amplicons per mutation
                    bin1 <- (mut - amp_length):(mut - 1)
                    b_list[[1]] <- bin1
                    bin2 <- (mut - (amp_length - 1)):(mut)
                    b_list[[2]] <- bin2
                    bin3 <- (mut + amp_length):(mut + 1)
                    b_list[[3]] <- bin3
                    bin4 <- (mut + (amp_length - 1)):(mut)
                    b_list[[4]] <- bin4
                    bin5 <- (mut - ceiling(amp_length / 2)):(mut - floor(amp_length / 2))
                    b_list[[5]] <- bin5
                    for (b in b_list){
                        low <- min(b)
                        up <- max(b)
                        mutations <- intersect(b, pos)
                        if (length(mutations > 0)){
                              mutations_keys <- make.keys(mutations)
                              count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                              row <- data.frame("lowerbound" = low, "upperbound" = up,
                                                "chromosome" = bins$chromosome[[i]], "count" = count,"id" = id)
                              possible_bins <- rbind(possible_bins, row)
                        }
                    }
                }
            }
        }
        return(possible_bins)       # all amplicons are added to dataframe
    }
    chrom_list <- unique(as.list(data$chr)) # list of all chromosomes found in this dataset
    all_pos_bins <- data.frame()
    for (chrom in chrom_list) {
        chromosome <- subset(data, chr == chrom)
        chromsort <- chromosome[order(chromosome$pos), ] #order df from lowest to highest position
        chrom_bins <- amplicon_finder(chromsort, chrom) #run through above function
        all_pos_bins <- rbind(all_pos_bins, chrom_bins) #add all bins together from all chromosomes
    }
    return(all_pos_bins)
}



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
#' @import hash
#' @import rlist
#' @import R.utils
#' @param bins A dataframe containing all potential amplicons
#' @param data A dataframe containing the location of each mutation.
#' @param amp The length of amplicons in number of base pairs
#' @param len The total length of sequencing panel in number of base pairs
#'
#' @return A dataframe containing the genomic coordinates for targeted sequencing panel
#'
#' @export
#'
#'
#'
#' @examples
#'
#' load(example.RData)
#' my_bins <- amp_pool(data, 100, 5000)
#'
#' fw_hotspot(my_bins, data, 100, 5000)


fw_hotspot <- function(bins, data, amp, len){
    panel_length = len #set panel and amplicon length
    amp_length = amp
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
                for (k in 2:nrow(spot_sub)){
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
        upper = final_bin$upperbound[[y]]
        lower = final_bin$lowerbound[[y]]
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
    for (u in 2:length(count_list)) {     #repeat for all other bins
        cum_mut_list[[u]] <- cum_mut_list[[u-1]] + count_list[[u]]  #add all calculated bin lengths and mutation counts to value one row above
        bin_len_cum[[u]] <- bin_len_cum[[u-1]] + bin_len_list[[u]]
    }
    final_bin$Cummulative_Bin_Length <- unlist(bin_len_cum)
    final_bin$Cummulative_Mutations <- unlist(cum_mut_list)
    final_bin <- final_bin[final_bin$Cummulative_Bin_Length <= panel_length,] ##apply bp length cutoff if needed
    final_bin <- final_bin[,c(1:4,7:8)]
    return(final_bin)
}



#' comprehensive amplicon ranking
#'
#' @description
#' create a targeted sequencing panel by finding which amplicons will likely
#' capture the most mutations using a pseudo-exhaustive selection method
#'
#' @details
#' Comprehensive Selection Sequencing Panel Identifier (Optimal mutation capture)
#'
#' 1.	To conserve computational power, the forward selection sequencing panel
#' identifier is run to determine the lowest number of
#' mutations per amplicon (mutation frequency) that need to be included in the
#' predetermined length sequencing panel.
#'   a.	any amplicon generated by the algorithm, which is less than this
#'   threshold value, will be removed.
#'
#' 2.	For the feasible exhaustive selection of amplicon combinations covering
#' hotspot areas larger than the predefined number of amplicons in length,
#' the algorithm breaks these large regions into multiple smaller regions.
#'   a.	The amplicons covering these regions are pulled from the amplicon pool,
#'   based on their unique IDs.
#'
#' 3.	The algorithm finds both the minimum number of amplicons overlap and
#' all positions with this value and identifies the region with the longest
#' continuous spot of minimum value.
#'   a.	The region is split at the center of this longest continuous minimum
#'   post values and continues the splitting process until all smaller regions
#'   are less than the “n” number amplicon length set by the user.
#'      i.	As this set number of amplicons decreases, the computation time
#'      required also often decreases.
#'
#' 4.	All amplicons contained in these bins are added back to the amplicon pool,
#' based on a new unique ID.
#'
#' 5.	Amplicons covering hotspots less than or equal to one amplicon length are
#' added to the final sequencing panel dataset.
#'
#' 6.	To determine the optimal combination of amplicons for each region, the number
#' of amplicons necessary for full coverage of the bin is calculated.
#'
#' 7.	A list is generated of every possible combination of n, number of amplicons,
#' needed. For each combination of amplicons:
#'   a.	 amplicons that would not meet the threshold of unique mutations are
#'   filtered out, and the number of all mutations captured by these amplicons
#'   is calculated.
#'   b.	the combination of amplicons that yields the highest number of mutations
#'   is added to the final sequencing panel.
#'
#' 8.	All amplicons in the final sequencing panel are ranked from highest to lowest
#' based on the number of mutations they cover.
#'
#' 9.	All amplicons capturing the number of mutations equal to the cutoff are
#' further ranked to favor amplicons that have mutations closer in location to
#' the center of the amplicon.
#'
#' 10.	Cumulative base-pair length and cumulative mutations covered by each
#' amplicon are calculated.
#'   a.	Depending on the desired length of the targeted panel, a cutoff may be
#'   applied to remove all amplicons which fall below a set cumulative length.
#'
#' @import hash
#' @import rlist
#' @import R.utils
#' @param bins A dataframe containing all potential amplicons
#' @param fw_panel a dataframe containing the sequencing panel designed by fw_hotspot
#' @param data A dataframe containing the location of each mutation.
#' @param amp The length of amplicons in number of base pairs
#' @param len The total length of sequencing panel in number of base pairs
#' @param size the threshold of hotspot size to split up in number of amplicons
#'
#' @return A dataframe containing the genomic coordinates for targeted sequencing panel
#'
#' @export
#'
#'
#' @examples
#'
#' load(example.RData)
#' my_bins <- amp_pool(data, 100, 5000)
#'
#' my_fw_panel <- fw_hotspot(my_bins, data, 100, 5000)
#'
#' com_hotspot(my_fw_panel, my_bins, data, 100, 5000, 3)
#'

com_hotspot <- function(fw_panel, bins, data, amp, len, size){
    panel_length = len #set panel and amplicon length
    amp_length = amp
    all_pos_bins <- bins
    num_amps <- size
    length <- nrow(fw_panel[fw_panel$Cummulative_Bin_Length <= len,])
    cutoff <- fw_panel$count[[length]]
    ordered_bin <- all_pos_bins[order(-all_pos_bins$count), ]
    big_spots <- ordered_bin[!ordered_bin$id == "x", ]
    unique_spots <- unique(as.list(big_spots$id))
    keys <- make.keys(data$pos)
    dict <- hash(keys = keys, values = data$Freq)
    possible_bins <- data.frame()
    maybe_possible_bins <- data.frame()
    y = 1
    new_pos <- pos
    new_big_spots <- data.frame()
    maybe_new_big_spots <- data.frame()
    too_big_split <- data.frame()
    new_big_split <- data.frame()
    for (spot in unique_spots){
        all_vec <- list()
        all_bin <- list()
        spot_sub <- subset(big_spots, id == spot) #subset df by spot id
        num_amplicons <- ceiling(length(spot_sub$mut_lowerbound[[1]]:spot_sub$mut_upperbound[[1]]) / amp_length)
        if (num_amplicons <= num_amps){
            new_big_spots <- rbind(new_big_spots, spot_sub)
        }
        if (num_amplicons > num_amps){
            spot_sub <- subset(spot_sub, count >= cutoff)         # to help split up very large hotspots, we filter out amplicons which to not meet the cutoff
            if (nrow(spot_sub) > 0){                              # mutation count established by forward binning method
                for (i in seq_len(nrow(spot_sub))){                        # hotspot is then split up if there are regions not covered by any amplicons
                    vec <- spot_sub$upperbound[[i]]:spot_sub$lowerbound[[i]]
                    all_vec <- c(all_vec, vec)
                    all_bin[[i]] <- vec
                }
                new_split <- as.data.frame(seqToIntervals(unique(all_vec)))
                too_big_split <- new_split
                for (j in seq_len(nrow(new_split))){
                    split_vec <- new_split$to[[j]]:new_split$from[[j]]
                    all_mutations <- intersect(split_vec, pos)
                    new_id <- paste(as.character(spot), as.character(j), sep = "-")
                    for (v in all_bin){
                        split_overlap <- intersect(v, split_vec)
                        if (length(split_overlap) == amp_length){
                            mutations <- intersect(v, pos) #find all mutations in total bin region
                            mutations_keys <- make.keys(mutations)
                            new_count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                            row <- data.frame("lowerbound" = min(v), "upperbound" = max(v), "count" = new_count,
                                      "chromosome" = spot_sub$chromosome[[1]], "id" = new_id, "mut_lowerbound" = min(all_mutations),
                                      "mut_upperbound" = max(all_mutations))
                            maybe_new_big_spots <- rbind(maybe_new_big_spots, row)
                        }
                    }
                }
                big_amp_ids <- unique(maybe_new_big_spots$id)
                for (p in big_amp_ids){
                    amp_list <- list()
                    amp_spot_sub <- subset(maybe_new_big_spots, id == p)
                    amp_test <- ceiling(length(amp_spot_sub$mut_lowerbound[[1]]:amp_spot_sub$mut_upperbound[[1]]) / 125)
                    if (amp_test <= num_amps){new_big_spots <- rbind(new_big_spots, amp_spot_sub)}
                    if (amp_test > num_amps){
                        too_big_split <- data.frame("to" = min(amp_spot_sub$mut_lowerbound), "from" = max(amp_spot_sub$mut_upperbound))
                        while (nrow(too_big_split) > 0){
                            new_amp_spot_sub <- data.frame()
                            too_big_list <- c()
                            for (i in seq_len(nrow(too_big_split))){
                                too_big_list <- c(too_big_list, too_big_split$from[[i]]:too_big_split$to[[i]])
                            }
                            for (v in seq_len(nrow(amp_spot_sub))){
                                vec <- amp_spot_sub$lowerbound[[v]]:amp_spot_sub$upperbound[[v]]
                                ov <- intersect(vec, too_big_list)
                                if (length(ov) == amp_length){new_amp_spot_sub <- rbind(new_amp_spot_sub, amp_spot_sub[v,])}
                            }
                            new_amp_spot_sub$mut_upperbound <- max(too_big_list)
                            new_amp_spot_sub$mut_lowerbound <- min(too_big_list)
                            all_big_vec <- list()
                            for (m in seq_len(nrow(new_amp_spot_sub))){
                                vec <- new_amp_spot_sub$upperbound[[m]]:new_amp_spot_sub$lowerbound[[m]]
                                all_big_vec <- c(all_big_vec, vec)
                                vec_freq <- data.frame(ftable(unlist(all_big_vec)))
                            }
                            low_point <- list()
                            for (x in seq_len(nrow(vec_freq))){
                                freq_list <- vec_freq$Freq
                                pos_list <- as.numeric(levels(factor(vec_freq$Var1)))
                                f <- freq_list[[x]]
                                minimum <- min(freq_list)
                                if (f == minimum){
                                    low_point <- c(low_point, pos_list[[x]])
                                }
                            }
                            low_point_split <- as.data.frame(seqToIntervals(low_point))
                            list_length <- list()
                            for (low in seq_len(nrow(low_point_split))){
                                length <- length(low_point_split$to[[low]]:low_point_split$from[[low]])
                                list_length <- c(list_length, length)
                            }
                            long_point <- which.max(list_length)
                            break_point <- ceiling(median(low_point_split$to[[long_point]]:low_point_split$from[[long_point]]))
                            all_big_vec <- vec_freq$Var1
                            all_big_vec <- as.numeric(levels(factor(all_big_vec)))
                            all_big_vec <- all_big_vec[all_big_vec != break_point]
                            test_big_split <- as.data.frame(seqToIntervals(unique(all_big_vec)))
                            too_big_split <- data.frame()
                            for (u in seq_len(nrow(test_big_split))){
                                vec <- test_big_split$to[[u]]:test_big_split$from[[u]]
                                len <- length(vec)
                                if (len <= (num_amps * amp_length)){
                                    new_big_split <- rbind(new_big_split, test_big_split[u,])
                                }
                                if (len > (num_amps * amp_length)){
                                    too_big_split <- rbind(too_big_split, test_big_split[u,])
                                }
                            }
                        }
                        for (j in seq_len(nrow(new_big_split))){
                            split_big_vec <- new_big_split$to[[j]]:new_big_split$from[[j]]
                            all_mutations <- intersect(split_big_vec, pos)
                            new_id <- paste(as.character(new_id), as.character(j), sep = "-")
                            for (v in all_bin){
                                split_overlap <- intersect(v, split_big_vec)
                                if (length(split_overlap) >= amp_length){
                                    mutations <- intersect(v, pos) #find all mutations in total bin region
                                    mutations_keys <- make.keys(mutations)
                                    new_count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                                    row <- data.frame("lowerbound" = min(v), "upperbound" = max(v), "count" = new_count,
                                                  "chromosome" = spot_sub$chromosome[[1]], "id" = new_id, "mut_lowerbound" = min(all_mutations),
                                                  "mut_upperbound" = max(all_mutations))
                                    new_big_spots <- rbind(new_big_spots, row)
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    all_best_combo <- list()
    fin_possible_bins <- data.frame()
    new_unique_spots <- unique(as.list(new_big_spots$id))
    for (spot in new_unique_spots){
        spot_sub <- subset(new_big_spots, id == spot) #subset df by spot id
        spot_sub <- subset(spot_sub, count >= cutoff)
        spot_sub <- unique(spot_sub)
        if (nrow(spot_sub) == 1){
            spot_sub <- spot_sub[,1:5]
            fin_possible_bins <- rbind(fin_possible_bins, spot_sub)
        }
        if (nrow(spot_sub) > 1){
            vec_list <- list()
            for (k in seq_len(nrow(spot_sub))){
                upper <- spot_sub$upperbound[[k]]
                lower <- spot_sub$lowerbound[[k]]
                vec <- upper:lower
                vec_list[[k]] <- vec
            }
            num_amplicons <- ceiling(length(spot_sub$mut_lowerbound[[1]]:spot_sub$mut_upperbound[[1]]) / amp_length)
            vec_list_num <- seq_len(length(vec_list))       # based on set length of large hotspots, generate list of every possible combination to amplicons
            if (num_amplicons == 1){
                pos_comb <- expand.grid(vec_list_num)
            }
            if (num_amplicons == 2){
                pos_comb <- expand.grid(vec_list_num, vec_list_num)
            }
            if (num_amplicons == 3){
                pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 4){
                pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 5){
                pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 6){
                pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 7){
                pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 8){
                pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 9){
                pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 10){
                pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            all_count <- list()
            fin_sum <- list()
            int_sum <- list()
            list_pos_bins <- list()
            for (c in seq_len(nrow(pos_comb))){           # use combination list to pull out amplicons from list
                com <- pos_comb[c,]                  # then find the total number of mutations captured by those amplicons
                un_com <- unique(unlist(com))        # pick combination which will capture the most mutations
                if (length(un_com) < length(com)){
                    pos_comb <- pos_comb[-c,]
                }
            }
            for (j in seq_len(nrow(pos_comb))){
                vec_comb <- list()
                temp_df <- data.frame()
                for (i in seq_len(ncol(pos_comb))){
                    vec <- unlist(vec_list[pos_comb[j,i]])
                    mutations <- intersect(vec, new_pos)
                    mutations_keys <- make.keys(mutations)
                    count <- sum(as.vector(values(dict, keys = mutations_keys)))
                    row <- data.frame("lowerbound" = min(vec), "upperbound" = max(vec),"count" = count)
                    temp_df <- rbind(temp_df, row)
                }
                maybe_bins_1 <- temp_df[order(-temp_df$count), ]
                maybe_bins_2 <- data.frame(maybe_bins_1) #
                maybe_bins_2$chromosome <- spot_sub$chromosome[[1]]
                maybe_bins_2$id <- spot_sub$id[[1]]
                list_pos_bins[[j]] <- maybe_bins_2
                for (x in seq_len(nrow(maybe_bins_2))){
                    sum_count <- list()
                    maybe_bins <- maybe_bins_2
                    maybe_bins[1,] <- maybe_bins_1[x,]
                    maybe_bins[x,] <- maybe_bins_1[1,]
                    big_bin_all <- list()
                    possible_bins <- data.frame()
                    while (nrow(maybe_bins) > 0){
                        row <- data.frame()
                        new_bin <- data.frame()
                        big_spot <- maybe_bins[1,]
                        big_bin <- big_spot$upperbound:big_spot$lowerbound
                        big_bin_all <- c(big_bin_all, big_bin)
                        possible_bins <- rbind(possible_bins, big_spot)
                        if (nrow(maybe_bins) > 1){
                            for (k in 2:nrow(maybe_bins)){
                                row <- data.frame()
                                test_bin <- maybe_bins$upperbound[k]:maybe_bins$lowerbound[k]
                                overlap <- intersect(big_bin_all, test_bin)
                                if (length(overlap) > 0){
                                    unique_bin <- setdiff(test_bin, overlap)
                                    mutations <- intersect(unique_bin, new_pos) #find all mutations in total bin region
                                    if (length(mutations) > 0){
                                        mutations_keys <- make.keys(mutations)
                                        new_count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                                        row <- data.frame("lowerbound" = min(test_bin), "upperbound" = max(test_bin), "count" = new_count,
                                                        "chromosome" = spot_sub$chromosome[[1]], "id" = spot_sub$id[[1]])
                                    }
                                }
                                if (length(overlap) <= 0) {row <- maybe_bins[k,]}
                                new_bin <- rbind(new_bin, row)
                            }
                        }
                        maybe_bins <- new_bin
                        if (nrow(maybe_bins > 0)){maybe_bins <- maybe_bins[order(-maybe_bins$count), ]}
                    }
                    testing_bins <- subset(possible_bins, count >= cutoff)
                    sum <- sum(unlist(testing_bins$count))
                    sum_count <- c(sum_count, sum)
                }
                fin_sum <- c(fin_sum, max(unlist(sum_count)))
                int_sum <- c(int_sum, which.max(sum_count))
            }
            best_combo <- which.max(fin_sum)
            best_one <- as.numeric(best_combo)
            bins <- list_pos_bins[[best_one]]
            best_int <- int_sum[[best_one]]
            maybe_bins <- bins
            maybe_bins[1,] <- bins[as.numeric(best_int),]               # additionally must favor combinations which contain the amplicons with unqiue mutations above cutoff value
            maybe_bins[as.numeric(best_int),] <- bins[1,]
            big_bin_all <- list()
            while (nrow(maybe_bins) > 0){
                row <- data.frame()
                new_bin <- data.frame()
                big_spot <- maybe_bins[1,]
                big_bin <- big_spot$upperbound:big_spot$lowerbound
                mutations <- intersect(big_bin, new_pos)
                mutations_keys <- make.keys(mutations)
                new_count <- sum(as.vector(values(dict, keys = mutations_keys)))
                big_bin_all <- c(big_bin_all, big_bin)
                fin_possible_bins <- rbind(fin_possible_bins, big_spot)
                if (nrow(maybe_bins) > 1){
                    for (k in 2:nrow(maybe_bins)){
                        row <- data.frame()
                        test_bin <- maybe_bins$upperbound[k]:maybe_bins$lowerbound[k]
                        overlap <- intersect(big_bin_all, test_bin)
                        if (length(overlap) > 0){
                            unique_bin <- setdiff(test_bin, overlap)
                            mutations <- intersect(unique_bin, new_pos)
                            if (length(mutations) > 0){
                                mutations_keys <- make.keys(mutations)
                                new_count <- sum(as.vector(values(dict, keys = mutations_keys)))
                                row <- data.frame("lowerbound" = min(test_bin), "upperbound" = max(test_bin), "count" = new_count,
                                            "chromosome" = spot_sub$chromosome[[1]], id = spot_sub$id[[1]])
                            }
                        }
                        if (length(overlap) <= 0) {row <- maybe_bins[k,]}
                        new_bin <- rbind(new_bin, row)
                    }
                }
                maybe_bins <- new_bin
                if (nrow(maybe_bins > 0)){maybe_bins <- maybe_bins[order(-maybe_bins$count), ]}
            }
        }
    }
    best_bin_all <- fin_possible_bins
    ordered_bin <- ordered_bin[-c(6:7)]
    final_bin <- rbind(subset(ordered_bin, id == "x"), best_bin_all)
    final_bin <- final_bin[order(-final_bin$count), ]
    cutoff_point <- final_bin$count[[round(panel_length / amp_length)]]
    subset_final_bin <- subset(final_bin, count == cutoff_point)
    dif_list <- list()
    for (i in nrow(subset_final_bin)){
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
    final_bin[df_range_min:df_range_max,] <- subset_final_bin[1:nrow(subset_final_bin),]
    bin_len_list <- list()
    for (y in seq_len(nrow(final_bin))) {      #calculating bin length using upper and lower bound positions
        upper = final_bin$upperbound[[y]]
        lower = final_bin$lowerbound[[y]]
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
    for (u in 2:length(count_list)) {     #repeat for all other bins
        cum_mut_list[[u]] <- cum_mut_list[[u-1]] + count_list[[u]]  #add all calculated bin lengths and mutation counts to value one row above
        bin_len_cum[[u]] <- bin_len_cum[[u-1]] + bin_len_list[[u]]
    }
    final_bin$Cummulative_Bin_Length <- unlist(bin_len_cum)
    final_bin$Cummulative_Mutations <- unlist(cum_mut_list)
    final_bin <- final_bin[final_bin$Cummulative_Bin_Length <= panel_length,] ##apply bp length cutoff if needed
    final_bin <- final_bin[,c(1:4,7:8)]
    return(final_bin)
}
