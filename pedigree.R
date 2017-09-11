library("kinship2")
jpeg('Pedigree_Brazil_High_Mendel.jpg')
pedigree <- read.table('/home/brunettt/Brazil_nuclear_mendel_errors/high_mendel_error_famlies_Brazilian_LDpruned.fam', header=TRUE)
ped <- pedigree(pedigree$PATIENT, pedigree$FATHER,
                pedigree$MOTHER, pedigree$SEX, pedigree$PHENO)
plot(ped,
     id = pedigree$PATIENT)
dev.off()

ibs <- read.table('/home/brunettt/Brazil_nuclear_mendel_errors/high_mendel_error_famlies_Brazilian_LDpruned.genome', header=TRUE)
ibs['EZ-PI_HAT'] <- ibs$EZ - ibs$PI_HAT
erroneous_relationships <- subset(ibs, EZ-PI_HAT>0.1)
