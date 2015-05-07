# Rでガウス分布の平均・精度パラメータを逐次ベイズ推定
# by data9824
# http://qiita.com/data9824/items/5378950c8a1c6437909e

frame()
set.seed(0)

# 試行回数
ite <- 10
g_cols <- 4
g_rows <- (ite + 2) / g_cols
par(mfrow=c(g_rows, g_cols))
par(mar=c(2.3, 2.5, 1, 0.1))
par(mgp=c(1.3, .5, 0))
murange <- c(-1, 1)
taurange <- c(0, 2)
mu <- seq(murange[1], murange[2], 0.01)
tau <- seq(taurange[1], taurange[2], 0.01)

# 観測点を生成する
N <- 15
MU <- 0
TAU <- 1
x <- rnorm(N, MU, sqrt(TAU ^ -1))
plot(x, xlab="n", main="samples")

# 平均・精度パラメータの事前分布を設定する
ngmu0 <- 0
nglambda0 <- 1.0E-6  # 平均パラメータの事前分布の分散(nglambda*λ)^-1を大きく取る
nga0 <- 1
ngb0 <- 1.0E-6  # 精度パラメータの事前分布の分散nga/ngbを大きく取る

# 真の事後分布であるガウスガンマ分布のパラメータを求める
uml <- mean(x)
sigma2ml <- var(x) * (N - 1) / N
nga <- nga0 + N / 2
ngb <- ngb0 + (N * sigma2ml + nglambda0 * N * (uml - ngmu0) ^ 2 / (nglambda0 + N)) / 2
ngmu <- (N * uml + nglambda0 * ngmu0) / (N + nglambda0)
nglambda <- nglambda0 + N

# 適当な初期値に設定する
muN <- 0.5
lambdaN <- 20
Emu <- muN
Emu2 <- lambdaN ^ -1 + muN ^ 2
aN <- 24
bN <- 16
Etau <- aN / bN

for (iteration in 0:ite) {
    if (iteration >= 1) {
        if (iteration %% 2 == 1) {
            # 因子 qmu(mu)=N(mu|muN,lambdaN) を求める
            muN <- (nglambda0 * ngmu0 + sum(x)) / (nglambda0 + N)
            lambdaN <- (nglambda0 + N) * Etau
            Emu <- muN
            Emu2 <- lambdaN ^ -1 + muN ^ 2
        } else {
            # 因子 qtau(tau)=Gam(tau|aN,bN) を求める
            aN <- nga0 + (N + 1) / 2
            bN <- ngb0 + (N + nglambda0) / 2 * Emu2 - (sum(x) + nglambda0 * ngmu0) * Emu +
                (sum(x ^ 2) + nglambda0 * ngmu0 ^ 2) / 2
            Etau <- aN / bN
        }
    }

    # 近似した分布を描画する
    q <- outer(mu, tau, function(mu, tau) {
        dnorm(mu, muN, sqrt(lambdaN ^ -1)) * dgamma(tau, aN, rate=bN)
        });
    contour(mu, tau, q, xlim=murange, ylim=taurange, xlab=expression(mu), ylab=expression(tau), col=4)

    # 真の分布を描画する
    p <- outer(mu, tau, Vectorize(function(mu, tau) 
        dnorm(mu, ngmu, sqrt(1 / (nglambda * tau))) * dgamma(tau, nga, rate=ngb)
        ))
    contour(mu, tau, p, xlab="", ylab="", col=3, add=T)

    legend("bottomright", c(
        expression(p(mu,tau)), 
        expression(q(mu,tau))
        ), col=c(3, 4), lty=1, bg="gray")
    if (iteration >= 1) {
        if (iteration %% 2 == 1) {
            title(bquote(paste("updated ", q[mu](mu), " #", .(iteration))))
        } else {
            title(bquote(paste("updated ", q[tau](tau), " #", .(iteration))))
        }
    } else {
        title("initial")
    }
}

