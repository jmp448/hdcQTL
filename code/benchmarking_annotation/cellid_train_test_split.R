## TRAIN TEST SPLIT
fca_full <- readRDS("/project2/gilad/jpopp/ebQTL/data/fca/counts.subsampled.mca_pca_fromhvgs.sce")
train_split <- sample(colnames(fca_full), 0.5*ncol(fca_full))
test_split <- setdiff(colnames(fca_full), train_split)

fca_train <- fca_full[,train_split]
fca_test <- fca_full[,test_split]

## Train a classifier on the training data only
fca_train.dec <- modelGeneVar(fca_train)
hvg.fca_train <- getTopHVGs(fca_train.dec, fdr.threshold=0.1)

fca_train.pca <- prcomp(t(logcounts(fca_train)[hvg.fca_train,]), rank=50)
reducedDims(fca_train) <- list(PCA=fca_train.pca$x)

fca_test.pca <- prcomp(t(logcounts(fca_test)[hvg.fca_train,]), rank=50)
reducedDims(fca_test) <- list(PCA=fca_test.pca$x)

### MCA
fca_train <- RunMCA(fca_train, nmcs=50, features=hvg.fca_train)
fca_test <- RunMCA(fca_test, nmcs=50, features=hvg.fca_train)

train_signatures <- GetGroupGeneSet(fca_train, group.by="celltype", n.features=100)

train_cellid_enrichments <- RunCellHGT(fca_train, pathways = train_signatures, features=hvg.fca_train, dims = 1:50, n.features=100)
test_cellid_enrichments <- RunCellHGT(fca_test, pathways = train_signatures, features=hvg.fca_train, dims = 1:50, n.features=100)

train_cellid_labels <- rownames(train_cellid_enrichments)[apply(train_cellid_enrichments, 2, which.max)]
train_cellid_labels <- ifelse(apply(train_cellid_enrichments, 2, max)>2, yes = train_cellid_labels, "unassigned") #3 unassigned cells here

test_cellid_labels <- rownames(test_cellid_enrichments)[apply(test_cellid_enrichments, 2, which.max)]
test_cellid_labels <- ifelse(apply(test_cellid_enrichments, 2, max)>2, yes = test_cellid_labels, "unassigned") #3 unassigned cells here

### View confusion matrix
celltypes <- c(unique(colData(fca_train)$celltype), "unassigned")
confusion_train <- confusionMatrix(data=factor(train_cellid_labels, levels=celltypes), reference = factor(colData(fca_train)$celltype, levels=celltypes))
confusion_test <- confusionMatrix(data=factor(test_cellid_labels, levels=celltypes), reference = factor(colData(fca_test)$celltype, levels=celltypes))

corrplot(confusion_train$table, is.corr=F, tl.pos='lt', tl.cex=0.75, method='color')
corrplot(confusion_test$table, is.corr=F, tl.pos='lt', tl.cex=0.75, method='color')

training_performance <- as_tibble(confusion_train$overall, rownames="criterion") %>%
  mutate(dataset="train")
testing_performance <- as_tibble(confusion_test$overall, rownames="criterion") %>%
  mutate(dataset="test")

bind_rows(training_performance, testing_performance) %>%
  mutate(dataset=factor(dataset, levels=c("train", "test")))
  ggplot(aes(x=criterion, y=value, fill=dataset)) +
  geom_bar(stat="identity", position=position_dodge())


save.image("/project2/gilad/jpopp/ebQTL/data/fca/train_test_split_analysis.Rdata")

