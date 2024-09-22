library(openxlsx)
library(readxl)
act=read.xlsx("win_dic.xlsx",rowNames = TRUE)
act[which(is.na(act),arr.ind = TRUE)]=0
g=read_excel("Adj_Matrix_Weight0.xlsx")
g=g[,-1]
g=as.data.frame(g)
names_g=read.table("myfile.txt",sep = ":")
rownames(g)=names_g[,1]
colnames(g)=names_g[,1]

indices=vector()
for(i in 1:length(rownames(g))){
  if(length(which(g[i,]!=0))==0){
    indices=c(indices,i)
  }
}
g=g[-indices,-indices]
save(g,file="g1.RData")
rm(g)

g=read_excel("Adj_Matrix_Weight1.xlsx")
g=g[,-1]
g=as.data.frame(g)
rownames(g)=names_g[,1]
colnames(g)=names_g[,1]

indices=vector()
for(i in 1:length(rownames(g))){
  if(length(which(g[i,]!=0))==0){
    indices=c(indices,i)
  }
}
g=g[-indices,-indices]
save(g,file="g2.RData")
rm(g)

g=read_excel("Adj_Matrix_Weight2.xlsx")
g=g[,-1]
g=as.data.frame(g)
rownames(g)=names_g[,1]
colnames(g)=names_g[,1]

indices=vector()
for(i in 1:length(rownames(g))){
  if(length(which(g[i,]!=0))==0){
    indices=c(indices,i)
  }
}
g=g[-indices,-indices]
save(g,file="g3.RData")
rm(g)


g=read_excel("Adj_Matrix_Weight3.xlsx")
g=g[,-1]
g=as.data.frame(g)
rownames(g)=names_g[,1]
colnames(g)=names_g[,1]

indices=vector()
for(i in 1:length(rownames(g))){
  if(length(which(g[i,]!=0))==0){
    indices=c(indices,i)
  }
}
g=g[-indices,-indices]
save(g,file="g4.RData")
rm(g)

g=read_excel("Adj_Matrix_Weight4.xlsx")
g=g[,-1]
g=as.data.frame(g)
rownames(g)=names_g[,1]
colnames(g)=names_g[,1]

indices=vector()
for(i in 1:length(rownames(g))){
  if(length(which(g[i,]!=0))==0){
    indices=c(indices,i)
  }
}
g=g[-indices,-indices]
save(g,file="g5.RData")
rm(g)

g=read_excel("Adj_Matrix_Weight5.xlsx")
g=g[,-1]
g=as.data.frame(g)
rownames(g)=names_g[,1]
colnames(g)=names_g[,1]

indices=vector()
for(i in 1:length(rownames(g))){
  if(length(which(g[i,]!=0))==0){
    indices=c(indices,i)
  }
}
g=g[-indices,-indices]
save(g,file="g6.RData")
rm(g)

for(j in 1:6){
  load(paste("g",j,".RData",sep=""))
  print(paste("J------",j))
  for(i in 1:dim(g)[1]){
    print(paste("I------",i))
    name=colnames(g)[i]
    for(l in 1:dim(g)[1]){
      name2=rownames(g)[l]
      g[name2,name]=g[name2,name]*act[j,name]*act[j,name2]
    }
  }
  save(g,file=paste("g",j,".RData",sep=""))
}

#Saving the network and identity matrix for input to DGI
dir.create("Networks")
dir.create("Identity_Matrix")
dir.create("No_Feature")
library(igraph)
for(i in 1:6){
  load(paste("g",i,".RData",sep=""))
  a=graph.adjacency(as.matrix(g),weighted=TRUE,mode="undirected")
  z=cbind(get.edgelist(a),E(a)$weight)
  write.table(as.data.frame(z),file=paste("Networks/GN.",i,".csv",sep=""),sep=",",row.names=FALSE,quote=FALSE)
  x=diag(nrow(g))
  colnames(x)=colnames(g)
  rownames(x)=rownames(g)
  write.csv(as.data.frame(x),file=paste("Identity_Matrix/GI",i,".csv",sep=""),sep=",",row.names=TRUE,col.names=TRUE,quote=FALSE)
}