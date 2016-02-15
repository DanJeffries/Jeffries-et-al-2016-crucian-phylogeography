#### Before and after histograms and stats for filtering hi Ho blacklisted loci ###

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/populations_r07_p17")

Pre_filter_covs <- read.delim("coverage.tsv", header = F)
length(Pre_filter_covs$V1)
summary(Pre_filter_covs$V1)


Post_filter_covs <- read.delim("../populations_r07_p17_Ho_filetered/coverage.tsv", header = F)
length(Post_filter_covs$V1)
summary(Post_filter_covs$V1)

par(mfrow = c(1,2))

hist(Pre_filter_covs$V1, breaks = max(Pre_filter_covs$V1), main = "Covereage per tag before Ho filtering", xlab = "Coverage")
textxy(50,1500, paste("mean coverage =", round(mean(Pre_filter_covs$V1), 3)), cex = 1, pos = 4)
textxy(50, 1400, paste("n Tags = ", length(Pre_filter_covs$V1)),cex = 1, pos = 4)

hist(Post_filter_covs$V1, breaks = max(Post_filter_covs$V1), main = "Covereage per tag after Ho filtering", xlab = "Coverage")
textxy(50,1400, paste("mean coverage =", round(mean(Post_filter_covs$V1), 3)), cex = 1, pos = 4)
textxy(50, 1300, paste("n Tags = ", length(Post_filter_covs$V1)),cex = 1, pos = 4)
#### Can also have a look to see how many snps are in each locus - look at the number for the kept loci, and look at the number for the blacklisted (homeolog) loci.
