data("mutation_data")
amps <- amp_pool(data = data, amp = 100)
fw_bins <- fw_hotspot(bins = amps, data = data, amp = 100, len = 1000, include_genes = TRUE)
com_bins <- com_hotspot(fw_panel = fw_bins, bins = amps, data = data,
                        amp = 100, len = 1000, size = 3, include_genes = TRUE)

test_that("comprehensive hotspot panel is valid", {
  measured <- com_bins[nrow(com_bins),6]
  actual <- 0

  chr_list <- unique(com_bins$Chromosome)
  for (c in chr_list){
    bins_sub <- subset(com_bins, Chromosome == c)
    pos_list <- c()
    for (i in seq_len(nrow(bins_sub))){
      pos_list <- c(pos_list, bins_sub$Lowerbound[[i]]:bins_sub$Upperbound[[i]])
    }
    data_sub <- subset(data, chr == c & pos %in% pos_list)
    actual <- actual + nrow(data_sub)
  }
  expect_equal(measured, actual)
})
