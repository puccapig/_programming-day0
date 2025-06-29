#Rachel Yuan Nong Mikkelsen 2025-JUN-28

print("GENE:")
gene_ <- readline()

print("aa_seq_:")
aa_seq_ <- readline()

print("clinvar_:")
clinvar_ <- readline()
clinvar_ <- read.delim(clinvar_)
clinvar_ <- clinvar_[grepl(gene_, clinvar_[,1]),]

a <- sapply(strsplit(clinvar_[,1], " "), function(x){x[2]})
a <- as.data.frame(a)
a[,2] <- sapply(strsplit(a[,1], ""), function(x){x[1]})
a <- a[which(a[,2] == "("),]

print("start cor fol:")
start_fol_ <- readline()
start_fol_ <- as.numeric(as.character(start_fol_))

print("feature length:")
feat_len_ <- readline()
feat_len_ <- as.numeric(as.character(feat_len_))

b <- 0
for(ii in 1:feat_len_){
	b[ii] <- start_fol_ + ii
}

ab <- 0
for(ii in 1:feat_len_){
        if(dim(as.data.frame(table(grepl(b[ii], a[,1]))))[1]==2){
        	ab[ii] <- a[grepl(b[ii], a[,1]),1][1]
        }
        if(dim(as.data.frame(table(grepl(b[ii], a[,1]))))[1]!=2){
        	ab[ii] <- -1
        }
}

print(ab)

ab <- as.data.frame(ab)
ab <- ab[which(ab[,1] != "-1"),]
ab <- as.data.frame(ab)
ab[,2] <- sapply(substr(ab[,1], 7, nchar(ab[,1])-1), function(x){x[1]})

for(aa_pos_ii in 1:dim(ab)[1]){
        ab_pos1_ <- ab[aa_pos_ii, 2]
        ab_pos2i_ <- 0
        for(ab_pos2_ii in 1:nchar(ab_pos1_)){
                ab_pos2_sub <- sapply(strsplit(ab_pos1_, ""), function(x){x[ab_pos2_ii]})
                if(ab_pos2_sub == 1 | ab_pos2_sub == 2 | ab_pos2_sub == 3 | ab_pos2_sub == 4 | ab_pos2_sub == 5 | ab_pos2_sub == 6 | ab_pos2_sub == 7 | ab_pos2_sub == 8 | ab_pos2_sub == 9 | ab_pos2_sub == 0){ 
                        ab_pos2i_ <- paste0(ab_pos2i_, ab_pos2_sub)
                }
                if(ab_pos2_sub == "_"){
                        break
                }
	}	
	ab_pos2i_ <- sapply(substr(ab_pos2i_, 2, nchar(ab_pos2i_)), function(x){x[1]})
	ab[aa_pos_ii,2] <- ab_pos2i_	
}
ab[,2] <- as.numeric(as.character(ab[,2]))

print(ab)

prop2_ <- matrix(ncol = 1, nrow = feat_len_)
for(ii in 1:feat_len_){
	prop2_[ii, 1] <- sapply(strsplit(aa_seq_, ""), function(x){x[ii]})
}

source("~col_to_aa_.r")

plot_ <- c(1:feat_len_)
plot_[1:feat_len_] <- 3

pdf("1.pdf", height = 3)

par(cex = 0.7)
plot(plot_, ylim = c(2, 6), pch = 19, col = col_to_aa_)

print(length(plot_))

par(cex = 1)
abline(v = ab[,2] - start_fol_, col = "green", lty = 3)

par(cex = 0.7)
title(sub = "clinvar_ |NCBI")
legend("topright", legend = gene_, text.col = "purple", bty = "n")
legend("topleft", legend = prop_aa_, text.col = aa_color_, bty = "n")

print("screen_aa:")
screen_aa_ <- readline()
screen_aa_ <- as.numeric(as.character(screen_aa_))

if(feat_len_ > 100){
	for(feat_ii in 1:100){
	#for(feat_ii in 1:1){
		part_start_ <- (feat_ii - 1)*screen_aa_ + 1
		part_end_ <- part_start_ + screen_aa_ - 1
		
		if(part_start_ > nchar(aa_seq_)){
			dev.off()
			break
		}

		if(part_end_ > nchar(aa_seq_)){
			part_end_ <- nchar(aa_seq_)
		}

		part_seq_ <- sapply(substr(aa_seq_, part_start_, part_end_), function(x){x[1]})
		part_fol_ <- 0
		part_len_ <- nchar(part_seq_)
		
		b <- 0
		for(part_ii in 1:part_len_){
        		b[part_ii] <- part_start_ + part_ii
		}

		ab <- 0
		for(part2_ii in 1:part_len_){
        		if(dim(as.data.frame(table(grepl(b[part2_ii], a[,1]))))[1]==2){
        			ab[part2_ii] <- a[grepl(b[part2_ii], a[,1]),1][1]
        		}
        		if(dim(as.data.frame(table(grepl(b[part2_ii], a[,1]))))[1]!=2){
        			ab[part2_ii] <- -1
        		}
		}	

		ab_i <- ab[grepl("-1", ab)]
		ab_i <- length(ab_i)
		if(ab_i != 100){
			ab <- as.data.frame(ab)
			ab <- ab[which(ab[,1] != "-1"),]
			ab <- as.data.frame(ab)
			ab[,2] <- sapply(substr(ab[,1], 7, nchar(ab[,1])-1), function(x){x[1]})

			for(aa_pos_ii in 1:dim(ab)[1]){
        			ab_pos1_ <- ab[aa_pos_ii, 2]
        			ab_pos2i_ <- 0
        			for(ab_pos2_ii in 1:nchar(ab_pos1_)){
                			ab_pos2_sub <- sapply(strsplit(ab_pos1_, ""), function(x){x[ab_pos2_ii]})
					if(ab_pos2_sub == 1 | ab_pos2_sub == 2 | ab_pos2_sub == 3 | ab_pos2_sub == 4 | ab_pos2_sub == 5 | ab_pos2_sub == 6 | ab_pos2_sub == 7 | ab_pos2_sub == 8 | ab_pos2_sub == 9 | ab_pos2_sub == 0){ 
                            ab_pos2i_ <- paste0(ab_pos2i_, ab_pos2_sub)
					}
					if(ab_pos2_sub == "_"){
                            break
                    }

				}
				ab_pos2i_ <- sapply(substr(ab_pos2i_, 2, nchar(ab_pos2i_)), function(x){x[1]})
				ab[aa_pos_ii,2] <- ab_pos2i_
			}
			ab[,2] <- as.numeric(as.character(ab[,2]))


			prop2_ <- matrix(ncol = 1, nrow = part_len_)
			for(part3_ii in 1:part_len_){
        			prop2_[part3_ii, 1] <- sapply(strsplit(part_seq_, ""), function(x){x[part3_ii]})
			}

			source("/Users/rachelnong/Desktop/paper_codes1/col_to_aa_.r")
	
			plot_ <- c(1:part_len_)
			plot_[1:part_len_] <- 3

			par(cex = 0.7)
			plot(c(part_start_:part_end_), plot_, ylim = c(2, 6), xlim = c(part_start_, part_end_), pch = 19, col = col_to_aa_, xlab = "aa pos")

			par(cex = 1)
			v_ <- ab[,2]
			abline(v = v_ - part_fol_, col = "green", lty = 3)
			#text(v_, 5, label = v_, col = "pink")

			par(cex = 0.7)
			sub_ <- paste0("clinvar_ |NCBI ", gene_, " | aa_pos: ", part_start_, " - ", part_end_, " of total ", nchar(aa_seq_), "aa") 
			title(sub = sub_)
	
			legend1_ <- paste0(gene_, " (", part_start_, "-",part_end_, ")")
			legend("topright", legend = legend1_, text.col = "purple", bty = "n")
			legend("topleft", legend = prop_aa_, text.col = aa_color_, bty = "n")
	
			legend2_ <- nchar(aa_seq_) - part_end_
        		legend2_ <- paste0("next ", legend2_, " aa")
       			print(legend2_)	
		}
	}
}

dev.off()
	



