#1. Missing value imputation
tempData <- mice(rawdata,m=5,maxit=50,meth='rf')

#2.Enrichment analysis
go <-enrichGO(gene=gene$ENTREZID,
               OrgDb='org.Hs.eg.db',
               ont="all",
               pAdjustMethod = "BH",
               readable = TRUE)

KEGG <- enrichKEGG(gene= gene$ENTREZID,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)
KEGG <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


#3. Subgroup analysis
#3.1 Correlation analysis
corr = cor(rawdata, method = c("spearman"))

corrplot.mixed(corr,lower='square',upper='pie',tl.col="black",
               lower.col=col1(200),upper.col=col1(200),tl.pos = 'lt')

#3.2 Mfuzz analysis
DEGs_exp_averp<-as.matrix(rawdata)
eset <- new("ExpressionSet",exprs = DEGs_exp_averp)
eset <- standardise(eset)
m <- mestimate(eset)
cl <- mfuzz(eset, c = 4, m = m) 
mfuzz.plot(eset,cl,mfrow=c(2,3),new.window= FALSE,time.labels=colnames(DEGs_exp_averp),colo = color.2)

#3.3 Visualization
yarrr::pirateplot(formula = RPS20.exon1.AP ~ Group+Age, 
                  data = subgroup,                  
                  main = "Age",                  
                  xlab = "Age",                  
                  ylab = "PSI",
                  gl.lwd = 0,
                  bean.f.o=1,
                  bean.f.col=c("#E65164","#66C2A5"),
                  inf.f.o=0.3,
                  point.o = 1,
                  point.col="azure4",
                  point.bg=0.05,
                  avg.line.lwd = 1,
                  theme=3)

#4. Classifiers
#4.1 Set division
set.seed(791)
training.rows <- sample(rownames(rawdata), dim(rawdata)[1]*0.8)
training <- rawdata[training.rows, ]

test.rows <- setdiff(rownames(test), training.rows)
test<- all[test.rows, ]


#4.2 Variable selection
variable_selectction=c()
for (i in 2:41) {
  glm1<-glm(formula = Group~training[,i],data=training,family= "binomial")
  glm2<-summary(glm1)
  OR<-round(exp(coef(glm1)[2]),3)
  SE<-glm2$coefficients[,2][2]
  CI5<-round(exp(coef(glm1)[2]-1.96*SE),3) #95%CI
  CI95<-round(exp(coef(glm1)[2]+1.96*SE),3)
  CI<-paste0(CI5,"-",CI95)
  P<-glm2$coefficients[,4][2]
  HR=paste0(colnames(train)[i],"#",OR,"#",CI,"#",P)
  variable_selectction=rbind(variable_selectction,HR)
}
write.csv(variable_selectction,"variable_selectction.csv")

#4.3 Classifier construction
#Adaboost
adaboost <- boosting(Group~., data = training)

#naive bayesian regression
NB<-naiveBayes(Group~., training)

#elastic net regression
lasso.cv = cv.glmnet(as.matrix(training[,2:5]), training[,1], nfolds = 10)

#gradient boosting machine
gbmFit1 <- gbm.fit(as.matrix(training[,2:5]), training[,1])

#K-nearest neighbor
knn <- knn3(Group~., training)

#logistic regression
fit3 <- glm(formula = Group~.,data = training)

#neural networks
network<-neuralnet(Group~.,training)

#support vector machine
fit.svm<- svm(Group~.,training)

#Xgboost
dtrain <- xgb.DMatrix(data = list(data=Matrix(data.matrix(training[,c(2:5)]),sparse=T),
                                  label=as.numeric(as.character(training[,1])))$data,
                      label = list(data=Matrix(data.matrix(training[,c(2:5)]),sparse=T),
                                   label=as.numeric(as.character(training[,1])))$label) 
xgb <- xgboost(dtrain, params=best_param, nrounds=nrounds, nthread=6)

#4.4 bootstrapping
rsq<- function(formula,data,indices){
  d<- data[indices,]
  adaboost<- boosting(formula,data = d,boos=TRUE, mfinal=200)
  adaboost.pred<-predict(adaboost,data)
  return (adaboost.pred[["prob"]][,2])
}

results<- boot(data=training,statistic=rsq,R=1000,formula=Group~.)

#4.5 Visualization
ROC1 <- roc(training$Group,as.numeric(train$prevalue))
plot(ROC1, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), 
     max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE)




