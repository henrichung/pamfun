
a <- data.table::fread("bin_cooc_sp.csv")
b <- mltools::sparsify(a)

saveRDS(b, "bin_cooc_sp.RDS")