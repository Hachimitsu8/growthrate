setwd("C:/Users/hachisuka hiroaki/Documents/r/practice")

#品種×マーカー遺伝子型情報
geno.score <- read.csv("geno_score.csv", row.names = 1)

#表現型データ　品種×(籾長、アミロース含量)
pheno <- read.csv("pheno.csv",row.names = 1)

#位置情報、マーカー名、染色体、物理情報＆累積位置
map <- read("map.csv", row.names = 1)

#入力データphenoの最初の6行をチェック
head(pheno)

head(map) #入力データmapをチェック

x <- as.matrix(geno.score) - 1
y <- pheno[,1]

selector <- !is.na(y) #欠損データをチェック
x <- y[selector] #欠損サンプルを取り除く
y <- x[selector,] #欠損サンプルを取り除く

library(rrBLUP) #A.mat関数用　GWAS関数用
library(BLR) #BLR関数用
library(glmnet)

#集団構造の確認
amat <- A.mat(x) #rrBLUP
heatmap(amat, symm = T)

#OKモデル
qmat <- points[,1:6]
p.q <- rep(NA,ncol(x))
for(i in 1:ncol(x)){
   data <- data.frame(y,x[,i],qmat)
   model <- lm(y ~ ., data = data)
   p.q[i] <- summary(model)$coefficients[2,4]
}

#Bayesian LASSO + QKモデル
y.scaled <- scale(y)
model <- BLR(y = y.scaled, XF = qmat, XL = x, GF=list(ID=1:nrow(amat),A=amat),
             prior = list(
             varE=list(df=4, S=1),
             varU=list(df=4, S=1),
             lambda=list(shape=0.6, rate=1e-4, type='random',value=30)),
             nIter=5000,burnIn=1000,thin=1,saveAt="temp".)

plot(map$cum.pos,abs(model$bL), col = map$chr + 1, type = "h", main = "Bayesian LASSO + QKmodel",
     ylab = "Absolute values of coefficients", xlab = "Position (bp)")
     abline(h = 2,lty ="dotted")

     