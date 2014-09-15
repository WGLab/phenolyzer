.libPaths("/home/huiyang/R/x86_64-redhat-linux-gnu-library/3.0")
library("wordcloud");
args = commandArgs(TRUE);
disease_count = read.table(args[1]);
pal = brewer.pal(9,"Set1")
pal2 = brewer.pal(8,"Dark2")
pal = pal[-6]
pal = c(pal,pal2);
png(filename=paste(args[1],".png",sep=""), width = 2000, height = 1800);
wordcloud(disease_count[[1]],disease_count[[2]],c(20,3),1,600,F,T,0.10,pal,F,T);
dev.off();