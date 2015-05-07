# Rでガウス分布の平均・精度パラメータを逐次ベイズ推定
# by data9824
# http://qiita.com/data9824/items/be042ecf5eafda735d8e
# PRML 2.3.6に記載の通り、ガウス分布に従う観測値xが与えられたとき、そのガウス分布の平均・精度パラメータの共役事前分布を用いて、xの観測毎に事後分布のパラメータを逐次求めます。
# 以下のベイズ推定をそれぞれ行います。
# 
#     xの精度が既知として、共役事前分布としてガウス分布p(u)=N(u|nu,nv)を用いて、xの平均パラメータの事後分布を求めます。
#     xの平均が既知として、共役事前分布としてガンマ分布p(λ)=Gam(λ|ga,gb)を用いて、xの精度パラメータの事後分布を求めます。
#     xの平均と精度が未知として、共役事前分布としてガウス-ガンマ分布（正規-ガンマ分布）p(u,λ)=N(u|ngu,(ngbeta*λ)^(-1))Gam(λ|nga,ngb)を用いて、xの平均・精度パラメータの同時事後分布を求めます。
# 
# xのプロットに重ねて、以下を示します。
# 
#     xの精度が既知として、xの平均パラメータの事後分布の期待値と標準偏差
#     xの平均が既知として、xの精度パラメータの事後分布の期待値と標準偏差を、xにおける標準偏差に変換して既知のxの平均パラメータに足したもの
#     xの平均と精度が未知として、xの平均・精度パラメータの事後分布の平均周辺分布の期待値と標準偏差
#     xの平均と精度が未知として、xの平均・精度パラメータの事後分布の精度周辺分布の期待値と標準偏差を、xにおける標準偏差に変換して、事後分布の平均周辺分布の期待値に足したもの
# 
# また、一定の観測毎に、その時点の平均パラメータ事前分布、精度パラメータ事前分布、平均・精度パラメータ同時事前分布を、その事前分布の期待値と標準偏差、およびxを生成した真の値とともに示します。
# 
# なお、事前分布には、無情報事前分布に近くなるよう、大きな分散を持つ分布を設定しています。

frame()
set.seed(0)
par(mfrow=c(4, 6))
layout(matrix(c(rep(1, 6), t(matrix(2:(3 * 6 + 1), 1))), 4, byrow=T), heights=c(3, 1, 1, 1))
par(mar=c(2.3, 2.5, 1, 0.1))
par(mgp=c(1.3, .5, 0))
urange <- c(148, 152)
lrange <- c(0, 0.5)
x1 <- seq(urange[1], urange[2], 0.1)
x2 <- seq(lrange[1], lrange[2], 0.01)

# サンプルの生成
N <- 100
U <- 150
V <- 4
x <- rnorm(N, U, sqrt(V))
xrange <- c(min(x), max(x))
# 平均パラメータの事前分布
# nu <- 0
# nv <- 1.0E+6  # 平均パラメータの事前分布の分散nvを大きく取る
# 精度パラメータの事前分布
ga <- 1
gb <- 1.0E-6  # 精度パラメータの事前分布の分散ga/gbを大きく取る
# 平均・精度パラメータの事前分布
ngu <- 0
ngbeta <- 1.0E-6  # 平均パラメータの事前分布の分散(ngbeta*λ)^-1を大きく取る
nga <- 1
ngb <- 1.0E-6  # 精度パラメータの事前分布の分散nga/ngbを大きく取る

d <- data.frame(nu=nu, nv=nv, ga=ga, gb=gb, ngu=ngu, ngbeta=ngbeta, nga=nga, ngb=ngb)
for (n in 1:N) {
    Nml <- 1
    uml <- x[n] / Nml
    # vml <- (x[n] - U) ^ 2 / Nml
    # 平均・精度パラメータ(u,λ)の事後分布N(u|ngu,(ngbeta*λ)^(-1))Gam(λ|nga,ngb)を更新
    nga <- nga + Nml / 2
    ngb <- ngb + (Nml * 0 + ngbeta * Nml * (uml - ngu) ^ 2 / (ngbeta + Nml)) / 2
    ngu <- (Nml * uml + ngbeta * ngu) / (Nml + ngbeta)
    ngbeta <- ngbeta + Nml
    d <- rbind(d, c(nu, nv, ga, gb, ngu, ngbeta, nga, ngb))
}

# 平均・精度パラメータ事前分布の平均周辺分布の分散
d$nguvar <- d$ngb / (d$ngbeta * (d$nga - 1))
# 平均・精度パラメータ事前分布の精度周辺分布の期待値
d$nglexp <- d$nga / d$ngb
# 平均・精度パラメータ事前分布の精度周辺分布の分散
d$nglvar <- d$nga / d$ngb ^ 2

# 観測値
plot(x, pch=19, cex=1, xlim=c(1, N), ylim=xrange, xlab="n")
abline(h=U)
abline(h=U + sqrt(V) * c(-1, 1), lty=2)
# xの平均・精度パラメータの事後分布の平均周辺分布の期待値と標準偏差
par(new=T)
plot(d$ngu[-1], type="o", col=4, xlim=c(1, N), ylim=xrange, xlab="", ylab="", axes=F)
par(new=T)
plot((d$ngu - sqrt(d$nguvar))[-1], type="l", lty=2, col=4, xlim=c(1, N), ylim=xrange, xlab="", ylab="", axes=F)
par(new=T)
plot((d$ngu + sqrt(d$nguvar))[-1], type="l", lty=2, col=4, xlim=c(1, N), ylim=xrange, xlab="", ylab="", axes=F)
# xの平均・精度パラメータの事後分布の精度周辺分布の期待値と標準偏差
par(new=T)
plot((d$ngu + sqrt(1 / d$nglexp))[-1], type="o", col=6, xlim=c(1, N), ylim=xrange, xlab="", ylab="", axes=F)
par(new=T)
plot((d$ngu + sqrt(1 / (d$nglexp + sqrt(d$nglvar))))[-1], type="l", lty=2, col=6, xlim=c(1, N), ylim=xrange, xlab="", ylab="", axes=F)
par(new=T)
plot((d$ngu + sqrt(1 / (d$nglexp - sqrt(d$nglvar))))[-1], type="l", lty=2, col=6, xlim=c(1, N), ylim=xrange, xlab="", ylab="", axes=F)
legend("bottomright", legend=c(
    "p(u,λ):E(u|X)±√(Var(u|X))", 
    "p(u,λ):E(u|X)+√(1/(E(λ|X)±√(Var(λ|X))))"), 
    col=c(4, 6), pch=1)

for (n in c(1, (1:10) * N / 10)) {
    # 平均・精度パラメータ事前分布とその各周辺分布の期待値と標準偏差
    p <- outer(x1, x2, Vectorize(function(x1, x2) 
        dnorm(x1, d$ngu[n], sqrt(1 / (d$ngbeta[n] * x2))) * dgamma(x2, d$nga[n], rate=d$ngb[n])
        ))
    image(x1, x2, p, xlab=expression(mu), ylab=expression(lambda), col=rainbow(450)[256:1])
    abline(v=d$ngu[n], col=4)
    abline(v=d$ngu[n] + sqrt(d$nguvar[n]) * c(-1, 1), col=4, lty=2)
    abline(h=d$nglexp[n], col=6)
    abline(h=d$nglexp[n] + sqrt(d$nglvar[n]) * c(-1, 1), col=6, lty=2)
    points(U, 1/V)
    title(paste0("p(u,λ)#", n), cex.main=1.0)
}

