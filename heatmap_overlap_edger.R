x<- as.matrix(overlap)
y <- melt(x)
library(ggplot2)
p <- ggplot(y, aes(y=X1, x=X2))
p + geom_tile(aes(fill=value), color="grey") + scale_fill_gradient(low="white", high="darkblue") + xlab("") + ylab("") + 
  theme(legend.position=c(0.9,0.2), legend.title = element_text(color="black", size=12,face="bold"), 
        legend.background=element_rect(fill="grey", size=0.8, linetype="solid"), 
        axis.text.x = element_text(angle=-90, face="bold"),
        axis.text.y = element_text(face="bold"),
        plot.title=element_text(size = 25, face = "bold", color = "black", vjust = -1)) + 
  ggtitle(expression(atop("Overlapped differentially expressed genes across cancer", atop(italic("Obtained by EdgeR"),""))))
save(p, file="overlap_genes_plot.RData")
load("overlap_genes_plot.RData")
ggsave("overlap_genes_plot.png", width=10, height=10)
