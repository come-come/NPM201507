vector<-read.table("G:/project2//NPM201507//code//1024vector.txt",header=TRUE,row.names = 1);
# rownames(vector)<-vector[,1]
# vector<-vector[,-1]

import_col<-vector[,1:3]
dis = dist(import_col, method="euclidean")
cluster3 = hclust(dis, method="average")
plot(cluster3,hang=-1)

# dis2 = dist(vector, method="euclidean")
# clust = hclust(dis2, method="average")
# plot(clust,hang=-1)

groups<-cutree(cluster3, k=8)
x<-cbind(groups,vector)
write.table (x, file ="G:/project2//NPM201507//code//R//col-3-all.txt", sep ="\t", row.names =TRUE, col.names =F, quote =F)
x1<- subset(x, groups==1)
write.table (x1, file ="G:/project2//NPM201507//code//R//col-3-1.txt", sep ="\t", row.names =TRUE, col.names =F, quote =F)
