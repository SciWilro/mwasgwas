full <- read.table('~/drive/research/prism/data/metadata/Anna-master-microbiome-2012-11-06-full-list.txt',sep='\t',head=T)

sampleids <- NULL; prismids <- NULL; diagnosis <- NULL

# add one "sample" per prism ID and per flora id
floraix <- grep('flora',colnames(full))
full <- apply(full, 2, as.character)
for(i in 1:nrow(full)){
    sampleids <- c(sampleids, full[i,'PRISM_Participant_ID'])
    prismids <- c(prismids, full[i,'PRISM_Participant_ID'])
    diagnosis <- c(diagnosis, full[i,'Diagnosis'])
    
    # iterate through flora IDs
    for(j in floraix){
        # if it's not blank or NA, add a sample
        if(!is.na(full[i,j]) && full[i,j] != ''){
            sampleids <- c(sampleids, full[i,j])
            prismids <- c(prismids, full[i,'PRISM_Participant_ID'])
            diagnosis <- c(diagnosis, full[i,'Diagnosis'])
        }
    }
}

newtable <- data.frame(SampleID=sampleids, PRISM_Participant_ID=prismids,
                        Diagnosis_0_CD_1_UC_2_IC_3_HC_4_DC=diagnosis)
sink('~/drive/research/prism/data/metadata/Anna-master-microbiome-2012-11-06-per-sample.txt')
write.table(newtable,sep='\t',quote=F,row.names=FALSE)
sink(NULL)
