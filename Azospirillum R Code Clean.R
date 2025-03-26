#read in data
geneData <- read.table("Genetics Data.txt", header = T)
pData <- read.table("Phosphate Solubilization Data.txt", header = T)

pDataNoNA <- pData[!is.na(pData$meanDiff), ]

#Load in libraries
library(dplyr)
library(car)


data_agg <- pData %>%
  group_by(strainNum) %>%
  summarise(
    meanPI = mean(meanDiff, na.rm = TRUE),
    SD = sd(meanDiff, na.rm = TRUE),
    gene = mean(sumGenes, na.rm = TRUE),
    colonyDiam = mean(colonyDiam1 + colonyDiam2 + colonyDiam3, na.rm = TRUE),
    cell_count = sqrt(mean(cellCount, na.rm = TRUE)),
    n = sum(!is.na(meanDiff))
  )

par(mfrow = c(1,3))

#colony diameter doesnt influence clearing size
plot(data_agg$meanPI, data_agg$colonyDiam,
     xlab = "Mean Colony Diameter",
     ylab = "Mean PI",
     pch = 19,
     main = "(A) Mean Colony Diameter Compared to Mean PI")

diamModel <- lm(meanPI ~ colonyDiam, data = data_agg)
abline(diamModel)
summary(diamModel)

#cell count doesnt influence clearing size
plot(sqrt(pData$cellCount),pData$meanDiff,
     xlab = "Cell Count",
     ylab = "Mean PI",
     pch = 19,
     main = "(B) Cell Count Compared to Mean PI")

cellCountModel <- lm(meanDiff ~ sqrt(cellCount), data = pData)
abline(cellCountModel)
summary(cellCountModel)

#total p-genes does not inflience clearing size
plot(pDataNoNA$sumGenes,pDataNoNA$meanDiff,
     xlab = "Sum of P-related Genes",
     ylab = "Mean PI",
     pch = 19,
     main = "(C) Sum of P-related Genes Compared to Mean PI")

geneModel <- lm(meanDiff ~ sumGenes, data = pDataNoNA)
abline(geneModel)
summary(geneModel)

#Variation
bp <- barplot(data_agg$meanPI,
              names.arg = data_agg$strainNum,
              col = "#A9A9A9",
              border = "white",
              main = "Mean Phosphate Solubilization for Each Strain",
              xlab = "Strain Number",
              ylab = "Mean Phosphate Solubilization Index (PI)",
              las = 2) 

arrows(x0 = bp, y0 = data_agg$meanPI - data_agg$SD, 
       x1 = bp, y1 = data_agg$meanPI + data_agg$SD,
       angle = 90, code = 3, length = 0.1, col = "black")



#No Significant Difference Between Sum of P-Related Genes
boxplot(geneData$sum_freq ~ geneData$clearing,
        ylab = "Sum of P-Related Genes",
        xlab = "Presence of Clearing",
        main = "A Sum of P-related Genes Compared Between Clearing and Non-Clearing Strains",
        col = c("#A9A9A9", "#E6E6E6"))

wilcox.test(geneData$sum_freq ~ geneData$clearing)

#Significant Difference Between Number of Unique P-Related Genes
boxplot(geneData$unique_level6 ~ geneData$clearing,
        ylab = "Number of Unique P-Related Genes",
        xlab = "Presence of Clearing",
        main = "Number of Unique P-Related Genes Compared Between Clearing and Non-Clearing Strains",
        col = c("#A9A9A9", "#E6E6E6"))

t.test(geneData$unique_level6 ~ geneData$clearing,var.equal=TRUE)


#Investigating Specific Genes
results <- data.frame(
  category = colnames(geneData)[2:(ncol(geneData) - 1)], # Adjust for phenotype column
  p_value = NA
)
for (i in 2:(ncol(geneData) - 1)) {
  category <- colnames(geneData)[i]
  results$p_value[i - 1] <- wilcox.test(geneData[[category]] ~ geneData$clearing)$p.value
}

print(results)


