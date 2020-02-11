#' ---
#' title: "Class 5: Data Visualization and graphs in R"
#' author: "Kendall Lin"
#' date: "2020-01-23"
#' ---

#Class 5
# Data Visualization

weight<- read.table("bimm143_05_rstats/weight_chart.txt", header = T)
plot(weight$Age, weight$Weight, type ="o", pch= 15, cex = 1.5, lwd = 2, ylim = c(2,10), xlab= "Age(months)", ylab = "Weight(kg)", main = "Baby weight with age", col = "blue")

#silly plot
plot(1:5, cex =1:5, pch =1:5)

#make a barplot
mouse<- read.table("bimm143_05_rstats/feature_counts.txt", header = T, sep= "\t")
par(mar= c(5,10,0,2))
barplot(mouse$Count, horiz = T, names.arg= mouse$Feature, las=1, col="light blue")

#hist
x <- c(rnorm(10000),rnorm(10000)+4)
hist(x, breaks=50)


#colors
color <- read.delim("bimm143_05_rstats/male_female_counts.txt", header = T)
barplot(color$Count, horiz = F,  names.arg = color$Sample, las= 2, col= rainbow(nrow(color)))
barplot(color$Count, horiz = F,  names.arg = color$Sample, las= 2, col= c(4,2))

#coloring by value
genes <- read.delim("bimm143_05_rstats/up_down_expression.txt", header = T)
table(genes$State)
palette(c("blue", "grey", "red"))
plot(genes$Condition1, genes$Condition2, col= genes$State, xlab = "Expression condition 1", ylab = "Expression condition 2")
levels(genes$State)

#Dynamic Use of Colors
meth <- read.delim("bimm143_05_rstats/expression_methylation.txt", header = T)
plot(meth$gene.meth, meth$expression, col = densCols(meth$gene.meth, meth$expression), pch=20)
zero <- meth$expression >0
plot(meth$gene.meth[zero], meth$expression[zero], col = densCols(meth$gene.meth[zero], meth$expression[zero]), pch=20)
colorchange <- colorRampPalette(c("blue", "green", "red", "yellow"))
plot(meth$gene.meth[zero], meth$expression[zero], col = densCols(meth$gene.meth[zero], meth$expression[zero], colramp = colorchange), pch=20)
     