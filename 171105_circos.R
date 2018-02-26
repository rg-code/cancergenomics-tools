skcm_maf_mutect <- GDCquery_Maf(tumor = "SKCM", pipelines = "mutect2", save.csv = TRUE)
skcm_maf_mutect <- skcm_maf_mutect[skcm_maf_mutect$Variant_Classification %in% c("Missense_Mutaion", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins"),]
mut.id <- paste0(skcm_maf_mutect$Chromosome, ":", skcm_maf_mutect$Start_Position, "-", skcm_maf_mutect$End_Position, "|", skcm_maf_mutect$Reference_Allele)
skcm_maf_mutect <- cbind(mut.id, skcm_maf_mutect)
numSamples <- table(mut.id)
s.mut <- names(which(numSamples>=2))
s.mut  <- skcm_maf_mutect[skcm_maf_mutect$mut.id %in% s.mut,]
s.mut <- s.mut[,c("Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Hugo_Symbol")]
s.mut <- unique(s.mut)
s.mut[,1] <-  as.character(s.mut[,1])
s.mut[,4] <-  as.character(s.mut[,4])
s.mut[,5] <-  as.character(s.mut[,5])
typeNames <- unique(s.mut[,4])
type <-  c(3:1)
names(type) <- typeNames[1:3]
Type <- type[s.mut[ ,4]]
s.mut <- cbind(s.mut,Type)
s.mut <- s.mut[,c(1:3,6,4,5)]

load("SKCM_CNV_results.rda")
s.cnv <- as.data.frame(RecCNV[RecCNV[, "q-value"]<=10^-2,c(1:4,6)])
s.cnv <- s.cnv[,c(1,3,4,2)]
xidx <- which(s.cnv$Chromosome==23)
yidx <- which(s.cnv$Chromosome==24)
s.cnv[xidx, "Chromosome"] <- "X"
s.cnv[yidx, "Chromosome"] <- "Y"
Chromosome <- sapply(s.cnv[,1],function(x) paste0("chr",x))
s.cnv <- cbind(Chromosome, s.cnv[,-1])
s.cnv[,1] <- as.character(s.cnv[,1])
s.cnv[,4] <- as.character(s.cnv[,4])
s.cnv <- cbind(s.cnv,CNV=1)
colnames(s.cnv) <- c("Chromosome","Start_position","End_position","Aberration_Kind","CNV")


library(circlize)
pdf("180130_circos_skcm.pdf", width=10, height=10)
par(mar=c(1,1,1,1), cex=1)
circos.initializeWithIdeogram()
colors <- c("forestgreen", "firebrick")
names(colors) <- c(0,1)
circos.genomicTrackPlotRegion(s.cnv, ylim = c(0,1.2),panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value,ytop.column = 2, ybottom = 0, col = colors[value[[1]]],
                     border="white")
  cell.xlim = get.cell.meta.data("cell.xlim")
  circos.lines(cell.xlim, c(0,0), lty = 2, col = "#00000040")
})
colors <- c("blue", "green", "red")
names(colors) <- typeNames[1:3]
circos.genomicTrackPlotRegion(s.mut, ylim = c(1.2,4.2), panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, cex=0.8, pch=16,  col = colors[value[[2]]], ...)
})
circos.clear()
legend(-0.2, 0.2, bty="n", y.intersp = 1, c("Amp","Del"), pch=15, col=c("firebrick", "forestgreen"), title="CNVs", text.font=3, cex=1.2, title.adj=0)
legend(-0.2, 0, bty="n", y.intersp = 1,names(colors), pch=16, col=colors, title="Mutations", text.font=3, cex=1.2, title.adj=0)
dev.off()