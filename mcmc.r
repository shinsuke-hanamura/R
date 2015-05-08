### load function
# setwd("C:/var/work/R")
# source("mcmc.r", encoding="UTF-8")


# 正規分布をMCMCでサンプリングする
r_mean <- 10
r_lambda <- 10
r_sd <- 1 / r_lambda
N <- 1000
# 正規分布を作成
r_x <- rnorm(N, r_mean, r_sd)
# hist(r_x)

# x <- c(970:1030) * 0.01
# y <- dnorm(x, r_mean, r_sd)

# 対象とする確率分布
mean <- r_mean
lambda <- r_lambda

sam_func <- function(arg_x, mean, lambda) {
  return(dnorm(arg_x, mean, sqrt(1/lambda)))
}

# サンプリング数
s_size <- 10000
x <- rep(0, s_size)
p <- rep(0, s_size)

# 初期値
x[1] <- 0
p[1] <- sam_func(x[1], mean, lambda)


# 候補値との距離
det <- 0.01

for(i in 1:s_size){
  # 採用判定
  reg_pro <- -1
  
  while (runif(1,0,1) > reg_pro) { # 候補値が採用されるまで
  
    # 候補の値を算出する
    x[i + 1] <- x[i] + (-1 + 2 * floor(runif(1,0,2))) * det
    p[i + 1] <- sam_func(x[i + 1], mean, lambda)
    reg_pro <- (p[i + 1] / p[i])
  }
}

# サンプルからBurn-in（最初から500）を除く
#x <- tail(x, length(x) - 500)
#p <- tail(p, length(p) - 500)

# 推定結果を表示

print(sprintf("mean:%1.5f variance:%1.5f", mean(x), var(x)))

# グラフ描画
g_cols <- 1
g_rows <- 1
par(mfrow=c(g_rows, g_cols))

x_range <- c(min(x), max(x))
# plot(x)


hist(x, xlim=x_range, xlab="" , ylab = "freq", main = "density of estimated distribution", axes = F, ann = F)  # axes, annは軸に何も描かないようにFalse設定
axis(4)
par(new = T) # 重ね合わせ
plot(x, p, xlim=x_range, ylim=c(0, max(p)),
     xlab="x" , ylab = "freq" , main="" , col = 3 , lty = 1, type = "o")







# ベイズ統計入門 ＆ モンテカルロ法と逆問題 I. ベイズ推定
# http://www.ism.ac.jp/~iba/kougi_2006_ism/c20062.pdf

# http://d.hatena.ne.jp/uncorrelated/20120128/1327752869

