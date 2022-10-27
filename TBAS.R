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
training.samples <- rawdata$Group %>%
  createDataPartition(p = 0.8,list = FALSE)
model <- rawdata[training.samples,]
test <- rawdata[-training.samples,]

set.seed(234)
model.samples <- model$Group %>%
  createDataPartition(p = 0.75,list = FALSE)
train <- model[model.samples,]
validation <- model[-model.samples,]

#4.2 Variable selection
traindata1 <- data.matrix(train[,c(3:40)]) 
traindata2 <- Matrix(traindata1,sparse=T) 
traindata3 <- as.numeric(as.character(train[,2]))
traindata4 <- list(data=traindata2,label=traindata3) 
dtrain <- xgb.DMatrix(data = traindata4$data, label = traindata4$label) 

params <- list(max_depth = 3, 
               objective = "binary:logistic",
               silent = 0)

watchlist <- list(train = dtrain, eval = dtest)

cv_model <- xgb.cv(params = params,
                   data = drain, 
                   nrounds = 100, 
                   watchlist = watchlist,
                   nfold = 5,
                   verbose = FALSE,
                   prediction = TRUE,
                   eval_metric = "logloss")

cv_model$evaluation_log %>%
  filter(test_logloss_mean == min(test_logloss_mean))
min_logloss <- min(cv_model$evaluation_log[, test_logloss_mean])
min_logloss_index <- which.min(cv_model$evaluation_log[, test_logloss_mean])

nround <- min_logloss_index
best_param <- params
xgb <- xgb.train(data=dtest, params=best_param, nrounds=nround, nthread=6)

names <- dimnames(data.matrix(testv[,c(3:40)]))[[2]]
importance_matrix <- xgb.importance(names,model=xgb)
write.csv(importance_matrix,"importance_xgboost_trainv.csv")

#4.3 Classifier construction
#Adaboost
wine_adaboost <- boosting(Group~ ., 
                          data = train,boos=TRUE, mfinal=200)

#bayesian regression
stan_glm1 <- stan_glm(Group ~ ., 
                      data = train, family = poisson, 
                      prior = normal(0,2.5), prior_intercept = normal(0,5))

#elastic net regression
modelA <- train(
  as.factor(Group) ~., data = train, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)

modelA$bestTune
coef(modelA$finalModel, modelA$bestTune$lambda)

#gradient boosting machine
gbmFit1 <- gbm.fit(x,y,distribution = "bernoulli",shrinkage = 0.01,
                   interaction.depth = 3,n.minobsinnode = 10)

#K-nearest neighbor
golf.tkknn <- train.kknn(Group~.,data = train,
                         kernel = c("rectangular", "triangular", "epanechnikov", "optimal"),
                         distance=1,scale=T)

#logistic regression
fit3 <- lrm(formula = Group~.,
            data = train)

#neural networks
network<-neuralnet(Group~.,train,hidden = 3)


#support vector machine
svmfit <- svm (Group ~ ., data = train, kernel = "radial", 
               cost = 2.1, gamma=0.2) # radial svm, scaling turned OFF

#Xgboost
xgb <- xgb.train(data=dtest, params=best_param, nrounds=nround, nthread=6)

#4.4 Visualization
ROC1 <- roc(train$Treat,as.numeric(train$prevalue))
plot(ROC1, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='XGboostROC曲线')




