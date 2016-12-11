setwd("C:/Users/hachisuka hiroaki/Documents/r/practice")

#�i��~�}�[�J�[��`�q�^���
geno.score <- read.csv("geno_score.csv", row.names = 1)

#�\���^�f�[�^�@�i��~(�����A�A�~���[�X�ܗ�)
pheno <- read.csv("pheno.csv",row.names = 1)

#�ʒu���A�}�[�J�[���A���F�́A������񁕗ݐψʒu
map <- read("map.csv", row.names = 1)

#���̓f�[�^pheno�̍ŏ���6�s���`�F�b�N
head(pheno)

head(map) #���̓f�[�^map���`�F�b�N

x <- as.matrix(geno.score) - 1
y <- pheno[,1]

selector <- !is.na(y) #�����f�[�^���`�F�b�N
x <- y[selector] #�����T���v������菜��
y <- x[selector,] #�����T���v������菜��

library(rrBLUP) #A.mat�֐��p�@GWAS�֐��p
library(BLR) #BLR�֐��p
library(glmnet)

#�W�c�\���̊m�F
amat <- A.mat(x) #rrBLUP
heatmap(amat, symm = T)

#OK���f��
qmat <- points[,1:6]
p.q <- rep(NA,ncol(x))
for(i in 1:ncol(x)){
   data <- data.frame(y,x[,i],qmat)
   model <- lm(y ~ ., data = data)
   p.q[i] <- summary(model)$coefficients[2,4]
}

#Bayesian LASSO + QK���f��
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

     