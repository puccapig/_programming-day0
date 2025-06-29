#Rachel Yuan Nong Mikkelsen 2025JUN20

#col_to_aa_ <- c(1:dim(data_)[1])
print("col_to_aa_ [dim[1]]")
dim1_ <- readline()
col_to_aa_ <- as.numeric(as.character(dim1_))
col_to_aa_ <- c(1:dim1_)

prop_aa_ <- c("Aromatic", "Negative.charged", "Nonpolar", "Polar.uncharged", "Positive.charged", "*")
aa_color_ <- c("orange", "red", "darkblue", "lightblue", "brown", "purple")

for(iii in 1:dim1_){
	if(prop2_[iii] == "F" |prop2_[iii] == "W" |prop2_[iii] == "Y"){
		col_to_aa_[iii] <- aa_color_[1]
        }
	if(prop2_[iii] == "D" |prop2_[iii] == "E" |prop2_[iii] == "Q"){
		col_to_aa_[iii] <- aa_color_[2]
        }
	if(prop2_[iii] == "A" |prop2_[iii] == "G" |prop2_[iii] == "I" |prop2_[iii] == "L" |prop2_[iii] == "M" |prop2_[iii] == "P" |prop2_[iii] == "V"){
		col_to_aa_[iii] <- aa_color_[3]
	}
	if(prop2_[iii] == "C" |prop2_[iii] == "N" |prop2_[iii] == "S" |prop2_[iii] == "T"){
		col_to_aa_[iii] <- aa_color_[4]
       	}
       	if(prop2_[iii] == "H" |prop2_[iii] == "K" |prop2_[iii] == "R"){
		col_to_aa_[iii] <- aa_color_[5]
        }
	if(prop2_[iii] == "*"){
                col_to_aa_[iii] <- aa_color_[6]
	}
}

print("col_to_aa_")


