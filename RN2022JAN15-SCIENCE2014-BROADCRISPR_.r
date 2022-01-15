#Rachel Yuan Nong Mikkelsen 2022-januari-15
#NBT2020-GSE132080__


sgRNA_barcode_table <- matrix(0, ncol = 20, nrow = 128)

for(i in 1:128){
        ki <- as.character(sg_barcodes[i])
                for(k in 1:20){
                        sgRNA_barcode_table[i,k] <- sapply(strsplit(ki, ""), function(x){x[k]})
                }
        }

mut_table_0_1 <- matrix(0, ncol = 20, nrow = 128)
mut_table_0 <- crispr_barcodes[crispr_barcodes[,6] == 1,]
mut_table_1 <- sgRNA_barcode_table[crispr_barcodes[,6] == 1,]
gene_table <- as.data.frame(table(crispr_barcodes[,3]))
k_0 <- 1
for(i in 1:25){
	ki <- as.numeric(gene_table[i,2])
	for(k in 1:ki){
		for(ii in 1:20){
			if(sgRNA_barcode_table[k_0, ii] != mut_table_1[i, ii]){
				mut_table_0_1[k_0, ii] <- -1
				print(paste0("(",k_0, ",", ii, ")"))
			}else{
				mut_table_0_1[k_0, ii] <- 0
			}
	k_0 <- k_0 + 1
	}
	}
	}

mut_table <- rowSums(mut_table_0_1)
table_mut_table <- as.data.frame(table(mut_table))
mut_table <- colSums(mut_table_0_1)
plot(mut_table)

mut_table_2 <- sgRNA_barcode_table[crispr_barcodes[,6] != 1,]
mut_nt_table_0_1 <- matrix(0, ncol = 20, nrow = 128)
k_0 <- 1
for(i in 1:dim(mut_table_2)[1]){
	ki <- as.numeric(gene_table[i,2])
	for(k in 1:ki){
		for(ii in 1:20){
			if(mut_table_0_1[k_0, ii] == -1){
				mut_nt_table_0_1[k_0, ii] <- sgRNA_barcode_table[k_0, ii]
			}	
		k_0 <- k_0 + 1
	}
	}
	}

mut_nt_table <- matrix(0, ncol = 20, nrow = 4)
ki_a <- c("T", "G", "C", "A")
for(i in 1:20){
	ki <- as.data.frame(table(mut_nt_table_0_1[,i]))
	ki <- ki[ki[,2] != 0,]
	for(k in 1:dim(ki)[1]){
		mut_nt_table[which(ki_a == ki[k,1]), i] <- ki[k,2]
		}
	}
rownames(mut_nt_table) <- ki_a
		
plot_table_1 <- matrix(0, ncol = 3, nrow = 20*4)
k_0 <- 1
for(i in 1:20){
        for(k in 1:4){
                plot_table_1[k_0, 1] <- i
                plot_table_1[k_0, 2] <- k
                plot_table_1[k_0, 3] <- mut_nt_table[k, i]
                k_0 <- k_0 + 1
                }
}

#mut_table_2 <- sgRNA_barcode_table[crispr_barcodes[,6] != 1,]
mut_table_1 <- sgRNA_barcode_table[crispr_barcodes[,6] == 1,]
mut_nt_table_0_1 <- matrix(0, ncol = 20, nrow = 128)
k_0 <- 1
for(i in 1:dim(mut_table_1)[1]){
        ki <- as.numeric(gene_table[i,2])
        for(k in 1:ki){
                for(ii in 1:20){
                        if(mut_table_0_1[k_0, ii] == 0){
                                mut_nt_table_0_1[k_0, ii] <- sgRNA_barcode_table[k_0, ii]
                        }
                k_0 <- k_0 + 1
        }
        }
        }

mut_nt_table <- matrix(0, ncol = 20, nrow = 4)
ki_a <- c("T", "G", "C", "A")
for(i in 1:20){
        ki <- as.data.frame(table(mut_nt_table_0_1[,i]))
        ki <- ki[ki[,2] != 0,]
        for(k in 1:dim(ki)[1]){
                mut_nt_table[which(ki_a == ki[k,1]), i] <- ki[k,2]
                }
        }
rownames(mut_nt_table) <- ki_a

plot_table_1 <- matrix(0, ncol = 3, nrow = 20*4)
k_0 <- 1
for(i in 1:20){
        for(k in 1:4){
                plot_table_1[k_0, 1] <- i
                plot_table_1[k_0, 2] <- k
                plot_table_1[k_0, 3] <- mut_nt_table[k, i]
                k_0 <- k_0 + 1
                }
}


matrix_in <- t(plot_table_1[,3])
rownames(matrix_in)[1] <- c("sgRNA_table")
source("/Users/rachelnong/Desktop/paperdna_codes/function_color_from_targets.R")
par(cex = 1)
plot(plot_table_1[,1], plot_table_1[,2], col = cols.by.raw.counts, pch = 19, ylim = c(0, 10), xlab = "position", ylab = "")
points(plot_table_1[,1], plot_table_1[,2] + 5, col = cols.by.raw.counts, cex = 1)
        par(cex = 0.5)
	for(i in 1:4){
                text(1, i, labels = paste0("to_", ki_a[i]), col = "blue")
        }
par(cex = 1)
text(6, 4.5, labels = "base was varied...", col = "purple")
	
	par(cex = 0.5)
        for(i in 1:4){
                text(1, i+5, labels = ki_a[i], col = "blue")
        }

par(cex = 1)
text(6, 10, labels = "at position from...", col = "purple")
