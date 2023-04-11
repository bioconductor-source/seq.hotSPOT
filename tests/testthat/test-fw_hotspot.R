data("mutation_data")
amps <- amp_pool(data = mutation_data, amp = 100)
fw_bins <- fw_hotspot(bins = amps, data = mutation_data, amp = 100, len = 1000, include_genes = TRUE)

test_that("forward hotspot panel is valid", {
  measured <- fw_bins[nrow(fw_bins),6]
  actual <- 0

  chr_list <- unique(fw_bins$Chromosome)
  for (c in chr_list){
    bins_sub <- subset(fw_bins, Chromosome == c)
    pos_list <- c()
    for (i in seq_len(nrow(bins_sub))){
      pos_list <- c(pos_list, bins_sub$Lowerbound[[i]]:bins_sub$Upperbound[[i]])
    }
    data_sub <- subset(mutation_data, chr == c & pos %in% pos_list)
    actual <- actual + nrow(data_sub)
  }
  expect_equal(measured, actual)
})
