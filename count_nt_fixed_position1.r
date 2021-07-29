#Rachel Yuan Nong Mikkelsen 2021-jul-29 beijing time
#examples: 10x cell barcodes

score = 0
i = 0
nt_0 = 0
nt_1 = 0


for( i in 1:nchar(ie_bc[1])){
        nt_0 <- sapply(strsplit(ie_bc[1], ""), function(x){x[i]})
        nt_1 <- sapply(strsplit(ie_bc[2], ""), function(x){x[i]})
        if(nt_0 == nt_1 & nt_0 != 0){
                score = score + 1
                }
        }

#print(paste0("round ", round(score/nchar(ie_bc[1])*100,2), "% nucleotides are the same"))


