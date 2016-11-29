geno.score <- read.csv("geno_score.csv", row.names = 1)
pheno <- read.csv("pheno.csv",row.names = 1)
map <- read("map.csv", row.names = 1)
 エラー:  関数 "read" を見つけることができませんでした 
map <- read.csv("map.csv", row.names = 1)
head(pheno)
  HULLED.SEED.LENGTH Amylose.Content
1               8.06           16.00
2               7.57           23.51
3               7.71           23.00
4               8.24           24.00
5               9.71           20.00
6               7.12           24.00
head(map)
          chr     pos cum.pos
id1000223   1  420422  420422
id1000556   1  655693  655693
id1000673   1  740153  740153
id1000830   1  913806  913806
id1000955   1 1041748 1041748
id1001073   1 1172387 1172387
x <- as.matrix(geno.score) - 1
y <- pheno[,1]
selector <- !is.na(y)
y <- y[selector]
x <- x[selector,]
library(rrBLUP)
 警告メッセージ: 
 パッケージ ‘rrBLUP’ はバージョン 3.3.2 の R の下で造られました  
library(BLR)
 要求されたパッケージ SuppDists をロード中です 
Package 'BLR', 1.4 (2014-12-03). 
Type 'help(BLR)' for summary information
 警告メッセージ: 
1:  パッケージ ‘BLR’ はバージョン 3.3.2 の R の下で造られました  
2:  パッケージ ‘SuppDists’ はバージョン 3.3.2 の R の下で造られました  
library(glmnet)
 library(glmnet) でエラー: 
   ‘glmnet’ という名前のパッケージはありません 
p.naive <- rep(NA, ncol(x))
for(i in 1:ncol(x)){
+   data <- data.frame(y, x[,i])
+   model <- lm(y ~ ., data = data)
+   p.naive[i] <- summary(model)$coefficients[2,4]}
plot(map$cum.pos, -log(p.naive, base = 10), col = map$chr + 1, type ="h", main = "Naive model",
+ xlab = "Position (bp)", ylab = "-log10(p)")
