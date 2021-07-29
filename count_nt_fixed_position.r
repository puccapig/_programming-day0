#Rachel Yuan Nong Mikkelsen 2020-JUL-29 beijing time
#see how many same nuleotides are present in two sequences of the same length

#examples: 10x cell barcodes

print("assign barcodes: ie_bc")
ie_bc <- sapply(strsplit(ie_bc, "-"), function(x){x[1]})

score = 0
i = 0
k0 = 0
k1 = 0
nt_0 = 0
nt_1 = 0
nt_r = 0
nt_c = 0
seq_mat1 = 0
seq_r = 0
seq_c = 0

for( i in 1:nchar(ie_bc[1])){
	nt_0 <- sapply(strsplit(ie_bc[1], ""), function(x){x[i]})
	nt_1 <- sapply(strsplit(ie_bc[2], ""), function(x){x[i]})
	if(nt_0 == nt_1 & nt_0 != 0){ 
		score = score + 1
		}
	}

print(paste0("seq1:", ie_bc[1]))
print(paste0("seq2:", ie_bc[2]))
print(paste0("comparing two sequences of nt: ", nchar(ie_bc[1])))
print(paste0("the number of the same nucleotides of these two sequences are: ", score))
print(paste0("round ", round(score/nchar(ie_bc[1])*100,2), "% nucleotides are the same"))

seq_mat1 <- matrix(0, nrow = nchar(ie_bc[1]), ncol = nchar(ie_bc[2]))

	seq_r <- strsplit(ie_bc[1], "")[[1]]
	seq_c <- strsplit(ie_bc[2], "")[[1]]

for(k0 in 1:length(seq_r)){
	for(k1 in 1:length(seq_c)){
		nt_r <- seq_r[k0]
		nt_c <- seq_c[k1]
		if(nt_r == nt_c & nt_r != 0){
			seq_mat1[k0, k1] = 1
		}
		if(nt_r != nt_c){
			seq_mat1[k0, k1] = 0
		}
	}
}

print(seq_mat1)
	
