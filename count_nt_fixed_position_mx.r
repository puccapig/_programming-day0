#Rachel Yuan Nong Mikkelsen 2021-JUL-29 beijing time
#

print("assign mx_seq")
mx_seq <- sapply(strsplit(mx_seq, "-"), function(x){x[1]})

mx_seq_score <- matrix(0, ncol = length(mx_seq), nrow = length(mx_seq))
mx_seq0 = 0
mx_seq1 = 0
k_0 = 0
k_1 = 0

for(k_0 in 1:length(mx_seq)){
	for(k_1 in 1:length(mx_seq)){
	mx_seq0 <- mx_seq[k_0]
	mx_seq1 <- mx_seq[k_1]
	ie_bc <- c(mx_seq0, mx_seq1)
	#source("count_nt_fixed_position.r")
	source("count_nt_fixed_position1.r")
	mx_seq_score[k_0, k_1] <- score
	}
}
