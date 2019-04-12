require(ggplot2)
require(reshape2)
args=commandArgs(T)
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

file<-read.csv(args[1],sep = '\t',header = T, row.names = 1)
file <- get_lower_tri(file)
table <- melt(as.matrix(file),na.rm = T)

ggplot(data = table, aes(Var1, Var2)) +
  geom_tile(aes(fill = value),colour = "white") +
  theme_classic() + 
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 90))+
  scale_fill_gradient(name="Jaccard Index", low = "white",high = "red")+
  labs(x = "", y = "")+
  theme(legend.title = element_text(size = 8))+
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 1.8)
ggsave(args[2])
print("Success!")

