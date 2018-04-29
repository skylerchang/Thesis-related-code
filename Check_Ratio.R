files <- list.files(path="result", pattern="*.txt", full.names=T, recursive=FALSE)
details = file.info(files)
details = details[with(details, order(as.POSIXct(mtime))), ]
files = rownames(details)
datalist <- list()
out33a=NULL
for (k in 1:length(files)) {
  test<-read.table(files[k])
  tb<-as.data.frame(test,row.names = NULL)
  colnames(tb)<-c("ratio","count")
  bb<-rep(tb$ratio,tb$count)
  tc<-table(bb)
  tc<-as.data.frame(tc,row.names = NULL)
  colnames(tc)<-c("ratio","count")
  if (tail(tc$ratio,1)==1) {
  rt<-tc$count/sum(tc$count)
  ne<-as.data.frame(rt,row.names = NULL)
  d1<-tail(ne,1)
  kk<-data.frame(k,d1)
  colnames(kk)<-c("id","ratio")
  out33a=rbind(out33a,kk,row.names=NULL)
}
}
ind<-out33a[order(out33a$ratio),]
ID<-paste0(ind$id[1])
id<-cat(ID)
