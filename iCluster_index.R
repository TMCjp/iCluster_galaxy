args <- commandArgs(TRUE)

total <- read.table(args[1], sep = ",", header = TRUE)

#通过列名匹配，来区分不同组学数据
samplenames <- colnames(total)
meth <- agrep("meth", samplenames, max = list(sub = 0))
#samplenames[meth] <- "meth"
miRseq <- agrep("miRseq", samplenames, max = list(sub = 0))
#samplenames[miRseq] <- "miRseq"
mutation <- agrep("mutation", samplenames, max = list(sub = 0))
#samplenames[mutation] <- "mutation"
protein <- agrep("protein",samplenames, max = list(sub = 0))
#samplenames[protein] <- "protein"
snp <- agrep("snp",samplenames, max = list(sub = 0))
#samplenames[snp] <- "snp"
mRNAseq_raw_counts <- agrep("mRNAseq_raw_counts",samplenames, max = list(sub = 0))
#samplenames[mRNAseq_raw_counts] <- "mRNAseq_raw_counts"
mRNAseq_RSEM <- agrep("mRNAseq_RSEM",samplenames, max = list(sub = 0))
#samplenames[mRNAseq_RSEM] <- "mRNAseq_RSEM"

#samplenames <- unique(samplenames[-1])
#as.name(samplenames[1])
rm(samplenames)

#单独提取symbol，并与不同组学数据进行结合，拆分数据为不同组学数据表格形式
genesymbol <- total[,1]

meth <- total[,meth]
meth <- cbind(genesymbol,meth)
miRseq <- total[,miRseq]
miRseq <- cbind(genesymbol,miRseq)
mutation <- total[,mutation]
mutation <- cbind(genesymbol,mutation)
protein <- total[,protein]
protein <- cbind(genesymbol,protein)
snp <- total[,snp]
snp <- cbind(genesymbol,snp)
mRNAseq_raw_counts <- total[,mRNAseq_raw_counts]
mRNAseq_raw_counts <- cbind(genesymbol,mRNAseq_raw_counts)
mRNAseq_RSEM <- total[,mRNAseq_RSEM]
mRNAseq_RSEM <- cbind(genesymbol,mRNAseq_RSEM)
rm(genesymbol)

#去除无效的组学数据
# if(length(colnames(meth)) < 2){rm(meth)}
# if(length(colnames(miRseq)) < 2){rm(miRseq)}
# if(length(colnames(mutation)) < 2){rm(mutation)}
# if(length(colnames(protein)) < 2){rm(protein)}
# if(length(colnames(snp)) < 2){rm(snp)}
# if(length(colnames(mRNAseq_raw_counts)) < 2){rm(mRNAseq_raw_counts)}
# if(length(colnames(mRNAseq_RSEM)) < 2){rm(mRNAseq_RSEM)}


#通过循环迭代，将不确定数量的组学数据组合成list类型
a <- 1
datalist <- list()
if(length(colnames(meth)) > 1){
  datalist[[a]] <- meth
  b <- a+1
} else {
  b <- a
}

if(length(colnames(miRseq)) > 1){
  datalist[[b]] <- miRseq
  c <- b+1
} else {
  c <- b
}

if(length(colnames(mutation)) > 1){
  datalist[[c]] <- mutation
  d <- c+1
} else {
  d <- c
}

if(length(colnames(protein)) > 1){
  datalist[[d]] <- protein
  e <- d+1
} else {
  e <- d
}

if(length(colnames(snp)) > 1){
  datalist[[e]] <- snp
  f <- e+1
} else {
  f <- e
}

if(length(colnames(mRNAseq_raw_counts)) > 1){
  datalist[[f]] <- mRNAseq_raw_counts
  g <- f+1
} else {
  g <- f
}

if(length(colnames(mRNAseq_RSEM)) > 1){
  datalist[[g]] <- mRNAseq_RSEM
  h <- g+1
} else {
  h <- g
}

rm(a,b,c,d,e,f,g,h)
rm(total)
rm(meth,miRseq,mutation,protein,snp,mRNAseq_raw_counts,mRNAseq_RSEM)

#对list对象进行基本的去重、去缺失操作
for(i in 1:length(datalist)){
  datalist1 <- datalist[[i]]
  #去除缺失
  datalist1 <- na.omit(datalist1)
  #去除重复（基本不用）
  uniq <- unique(datalist1[,1])
  math <- match(uniq,datalist1[,1])
  datalist1 <- datalist1[math,]
  rownames(datalist1) <- as.factor(datalist1[,1])
  #转置并变为matrix格式
  datalist1 <- as.matrix(t(datalist1[,-1]))
  #将matrix对象的行名（样本名）进行统一
  samplenames <- rownames(datalist1)
  samplenames <- substr(samplenames, 1, 15) #因为样本名一定是前15个字符（根据条形码）
  rownames(datalist1) <- samplenames
  #回放入list对象中
  datalist[[i]] <- datalist1
}
rm(datalist1,i,math,uniq,samplenames)


#判断组学数据量，超过4个组学则停止运行
if(length(datalist)>4) {
  stop("You supplied more than 4 datasets ", call.=FALSE)
}


#确定n.lambda的值（这里只是取了相对较好的值，具体参见icluster的官方指导文档）
f <- length(datalist)
if (f == 1) 
{ 
  nlam = 100 
} else if(f == 2) {
  nlam = 89
} else if(f == 3) {
  nlam = 101
} else if(f == 4) {
  nlam = 307
}

#all(rownames(datas[[1]])==rownames(datas[[2]]))


#外部读入type值，并进行必要的处理
types <- args[3]
if(length(datalist)==1){
  types = types
} else {
  types <- strsplit(types,',')
  types <- as.vector(types[[1]])
}

#外部读入最大聚类数
maxk <- as.numeric(args[4])

#生成随机数文件夹
Prefix <- runif(1,0,1000000000)
Prefix <- round(Prefix)
dir.create(paste0("/var/www/iclusteroutput/",Prefix,"fit"))

save(datalist,file=paste0("/var/www/iclusteroutput/",Prefix,"fit/", "datalist.Rdata"))

suppressMessages(library(iClusterPlus))
#这一步会用很长时间
set.seed(123)

if (length(datalist)==1) {
  for(k in 1:maxk){
    cv.fit = tune.iClusterPlus(cpus=8,dt1=datalist[[1]],
                               type= types,
                               K=k,n.lambda= nlam, #此处type下的为用户输入的参数，对应输入数据的类型
                               scale.lambda=1,maxiter=20)                   #scale.lambda下为用户输入的参数（参数默认为1，个数和数据个数对应）
    save(cv.fit, file=paste0("/var/www/iclusteroutput/",Prefix,"fit/", "cv.fit.m",k,".Rdata"))
  }
  #datas <- list(get(data[1]), get(data[2]), get(data[3]), get(data[4]))
}

if (length(datalist)==2) {
  for(k in 1:maxk){
    cv.fit = tune.iClusterPlus(cpus=8,dt1=datalist[[1]],dt2=datalist[[2]],
                               type= types,
                               K=k,n.lambda= nlam, #此处type下的为用户输入的参数，对应输入数据的类型
                               scale.lambda=c(1,1),maxiter=20)                   #scale.lambda下为用户输入的参数（参数默认为1，个数和数据个数对应）
    save(cv.fit, file=paste0("/var/www/iclusteroutput/",Prefix,"fit/", "cv.fit.m",k,".Rdata"))
  }
  #datas <- list(get(data[1]), get(data[2]), get(data[3]), get(data[4]))
}

if (length(datalist)==3) {
  for(k in 1:maxk){
    cv.fit = tune.iClusterPlus(cpus=8,dt1=datalist[[1]],dt2=datalist[[2]],dt3 = datalist[[3]],
                               type= types,
                               K=k,n.lambda= nlam, #此处type下的为用户输入的参数，对应输入数据的类型
                               scale.lambda=c(1,1,1),maxiter=20)                   #scale.lambda下为用户输入的参数（参数默认为1，个数和数据个数对应）
    save(cv.fit, file=paste0("/var/www/iclusteroutput/",Prefix,"fit/", "cv.fit.m",k,".Rdata"))
  }
  #datas <- list(get(data[1]), get(data[2]), get(data[3]), get(data[4]))
}

if (length(datalist)==4) {
  for(k in 1:maxk){
    cv.fit = tune.iClusterPlus(cpus=8,dt1=datalist[[1]],dt2=datalist[[2]],dt3 = datalist[[3]], dt4 = datalist[[4]],
                               type= types,
                               K=k,n.lambda= nlam, #此处type下的为用户输入的参数，对应输入数据的类型
                               scale.lambda=c(1,1,1,1),maxiter=20)                   #scale.lambda下为用户输入的参数（参数默认为1，个数和数据个数对应）
    save(cv.fit, file=paste0("/var/www/iclusteroutput/",Prefix,"fit/", "cv.fit.m",k,".Rdata"))
  }
  #datas <- list(get(data[1]), get(data[2]), get(data[3]), get(data[4]))
}


# 加载上一步存储下来的Rdata数据，并进行BIC、devR等的计算
setwd(paste0("/var/www/iclusteroutput/",Prefix,"fit/"))  #将工作目录设定到上一步生成的随机文件夹中，为了取得上面存储的cv.fit数据
output=alist()
files=grep("cv.fit.m",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)


# Now we get the ID for the lambda vector at which the BIC is minimum
# Then we obtain the deviance ratio of the lambda vector at which the BIC is minimum
minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}

clusters=getClusters(output)
rownames(clusters)=rownames(datalist[[1]])
colnames(clusters)=paste("K=",2:(length(output)+1),sep="")

clusters <- cbind(sample=rownames(clusters),clusters)
write.table(clusters, args[5],sep=',',row.names=FALSE,quote=F)

png(args[6])
plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",             #此处作图，根据此图让使用者选取合适的K值
     ylab="%Explained Variation")
dev.off()

#Prefix
write.table(Prefix, args[7],sep=',',row.names=FALSE,col.names=FALSE)
