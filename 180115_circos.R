skcm_maf_mutect <- GDCquery_Maf(tumor = "SKCM", pipelines = "mutect2", save.csv = TRUE)
skcm_maf_mutect <- skcm_maf_mutect[skcm_maf_mutect$Variant_Classification %in% c("Missense_Mutaion", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins"),]
mut.id <- paste0(skcm_maf_mutect$Chromosome, ":", skcm_maf_mutect$Start_Position, "-", skcm_maf_mutect$End_Position, "|", skcm_maf_mutect$Reference_Allele)
skcm_maf_mutect <- cbind(mut.id, skcm_maf_mutect)
numSamples <- table(mut.id)
s.mut <- names(which(numSamples>=1))
s.mut  <- skcm_maf_mutect[skcm_maf_mutect$mut.id %in% s.mut,]
s.mut <- s.mut[,c("Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Hugo_Symbol")]
s.mut[,1] <-  as.character(s.mut[,1])
s.mut[,4] <-  as.character(s.mut[,4])
s.mut[,5] <-  as.character(s.mut[,5])
typeNames <- unique(s.mut[,4])
type <-  c(3:1)
names(type) <- typeNames[1:3]
Type <- type[s.mut[ ,4]]
s.mut <- cbind(s.mut,Type)
s.mut <- s.mut[,c(1:3,6,4,5)]

library(circlize)
pdf("180205_circos_skcm.pdf", width=10, height=10)
par(mar=c(1,1,1,1), cex=1)
circos.initializeWithIdeogram()
circos.genomicTrack(bed2, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "h", baseline = 0, col = ifelse(value > 0, "red", "blue")) 
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

pdf("180212_allDNMTs_hyper_hypo.pdf", width=10,height=10)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))
circos.genomicDensity(hyper, col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(hypo, col = c("#0000FF80"), track.height = 0.1)
legend(-0.2,0.2, bty="n", y.intersp = 1,c("Hyper", "Hypo"), pch=16, col=c("#FF000080", "#0000FF80"), title="Methylation Level", text.font=3, cex=1.2, title.adj=0)