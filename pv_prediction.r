#

t <- c(10:30)
g <- c(1:140)*0.001
# 2015-04-30 00:00:00,-3.92166000,0.56664800,-5.14571000
a <- -3.92166000
b <- 0.56664800
c <- -5.14571000

p_tg <- function(a, b, c, t, g) {
  g_size <- length(g)
  t_size <- length(t)
  p <- matrix(0, nrow=g_size, ncol=t_size)  
  for (i in c(1 : g_size)) {
    for (j in c(1 : t_size)) {
      p[i,j] <- (a + b*t[j] + c*g[i]) * g[i]
    }
  }
  return(p)
}

p <- p_tg(a, b, c, t, g)

g_cols <- 1
g_rows <- 1
par(mfrow=c(g_rows, g_cols))

library(scatterplot3d)
library(rgl)

plot3d(g, t, p)

# scatterplot3d(g, t, p,
#               #xlim = c(min(g), max(g)), ylim = c(min(t), max(t)), zlim = c(min(p), max(p)),
#               color = 1,
#               col.axis="blue", col.grid="lightblue",
#               main = titles, type="h", pch=20) 

t_index <- 10
plot3d(g, t[t_index], p[,t_index])


# 太陽高度および太陽方位角の計算

d <- c(1 : 365)

# 日の係数
chi <- function (d) {
  return((2 * pi * (d - 1) ) / 365)
}

# 均時差
e_t <- function (chi) {
  return((0.0172 + 0.4281 * cos(chi) - 7.3515 * sin(chi) - 3.3495 * cos(2 * chi) - 9.3619 * sin(2 * chi) ) / 60)

}

# d <- c(1 : 365)

e_t_d <- e_t(chi(d))

plot(e_t_d, type="l")

# e_t_d_fft <- fft( e_t_d - mean( e_t_d ) )
# hz <- c(1:length(e_t_d_fft))
# plot( hz, abs(e_t_d_fft), type="l" )

# 経度
lambda <- 135

# 緯度
phi <- (pi / 180) * 35

# 時角
omega <- function(h, lambda, e_t) {
  return(15 * h + lambda - 135 + 15 * e_t)
}

h_noon <- 12
omega_d <- omega(h_noon, lambda, e_t_d)

# 太陽赤緯
delta <- function(chi) {
  return(0.006918 - 0.399912 * cos(chi) + 0.070257 * sin(chi) - 0.006758 * cos(2 * chi) + 0.000908 * sin(2 * chi))
}

delta_d <- delta(chi(d))


# 太陽高度
s_h <- function(phi, delta, omega) {
  return(asin(cos(phi) * cos(delta) * cos((pi / 180) * omega) + sin(delta) * sin(phi)))
}

s_h_d <- s_h(phi, delta_d, omega_d)


# 太陽方位角
psi <- function(delta, omega, s_h) {
  return(asin(cos(delta) * sin((pi / 180) * omega) / cos(s_h)))
}


psi_d <- psi(delta_d, omega_d, s_h_d)

plot(s_h_d, type="l" )
plot(psi_d, type="l" )
plot(s_h_d, psi_d)
plot(s_h_d, psi_d, type="l" )

plot((180 / pi) * s_h_d, type="l" )
plot((180 / pi) * psi_d, type="l" )


# simulation

# 経度
lambda <-135
# 緯度
phi <- (pi / 180) * 35

# gsub(" ", "", paste("000", d, seq=""))

fillzero <- function(zeros, str) {
  mode(str) <- "character"
  if (zeros < 0) return(str)
  zerostr <- rep(paste(rep("0", zeros), sep = "", collapse = ""), length(str))
  substr <- paste(zerostr, str, sep = "")
  return(gsub(" ", "", substr))
}

# fillzero(3,5)



# 時角を変化させる

d <- c(1 : 365)
dirPath <- "C:/var/work/data/solar"
dirPath2d <- "C:/var/work/data/solar/2d"
#pngFileName <- "fft_p_pv_d20.png"
#pngFileName <- paste("solar_", fillzero(3, d), ".png", sep = "")
#pngFilePath <- paste(dirPath, pngFileName, sep = "/")

#d <- 1

for (d in c(1 : 365)) {
  h <- c(0 : 23)
  e_t_d <- e_t(chi(d))
  omega_d <- omega(h, lambda, e_t_d)
  delta_d <- delta(chi(d))
  # 太陽高度
  s_h_d <- s_h(phi, delta_d, omega_d)
  # 太陽方位角
  psi_d <- psi(delta_d, omega_d, s_h_d)
  
  x <- sin(psi_d)
  y <- cos(psi_d)
  z <- sin(s_h_d)
  c <- (s_h_d > 0) + 1
  titles <- paste("Solar Altitude ", d, "th day", sep = "")
  
  
  # 描画デバイスを開く
  pngFileName <- paste("solar_", fillzero(3, d), ".png", sep = "")
  pngFilePath <- paste(dirPath, pngFileName, sep = "/")

  png(pngFilePath, width = 360, height = 360)  
  scatterplot3d(x, y, z,
                xlim = c(-2, 2), ylim = c(-2, 2), zlim = c(-2, 2),
                color = c,
                col.axis="blue", col.grid="lightblue",
                main = titles, type="h", pch=20) 
  # 描画デバイスを閉じる
  dev.off() 

  # 描画デバイスを開く
  pngFileName <- paste("solar_", fillzero(3, d), "_2d.png", sep = "")
  pngFilePath <- paste(dirPath2d, pngFileName, sep = "/")

  titles <- paste("2d Solar Altitude ", d, "th day", sep = "")
  png(pngFilePath, width = 360, height = 360)  
  plot(x, y,
                xlim = c(-2, 2), ylim = c(-2, 2), 
                col = c,
                main = titles, pch=20) 
  # 描画デバイスを閉じる
  dev.off() 

}


h <- c(0 : 23)
  d <- 1

x <- cos(psi_d)
y <- sin(psi_d)
z <- sin(s_h_d)

library(scatterplot3d)
scatterplot3d(x, y, z, highlight.3d=TRUE, col.axis="blue",
              col.grid="lightblue", main="Solar Altitude", pch=20) 
t1 <- 10
t2 <- 20

x_t <- x[t1:t2]
y_t <- y[t1:t2]
z_t <- z[t1:t2]

scatterplot3d(x_t, y_t, z_t, highlight.3d=TRUE, col.axis="blue",
              col.grid="lightblue", main="Solar Altitude", pch=20) 


# s_h_d_fft <- fft( s_h_d - mean( s_h_d ) )
# hz <- c(1:length(s_h_d_fft))
# plot( hz, abs(s_h_d_fft), type="l" )



# sin(pi / 6) = 0.5
# asin(0.5) = pi / 6

d <- c(1 : 365)

for (d in ()) {





x <- 60 * 24 * 365
f <- 1:x
hz <- 1:( x/2 )
t_all <- seq( (1 / x), 1.00,  1 / length( f ) )
t <- t_all[1:60 * 24 * 20]

# t <- seq( (1 / x), 1.00,  1 / length( f ) )



# Daily 360Hzの正弦波
data_waveA <- sin( 2 * pi * t * 360 )

# Anually 1Hzの正弦波
data_waveB <- sin( 2 * pi * t * 1 )

# 上記2つの合成波
data_waveC <- data_waveA + data_waveB

# 波を描きます。
par( "mfrow" = c( 3, 1 ) )
plot( t, data_waveA, type="l" )
plot( t, data_waveB, type="l" )
plot( t, data_waveC, type="l" )

# plot( t, data_waveA, xlim = c(0, 2/365), type = "l" )
# plot( t, data_waveC, xlim = c(0, 22/365), type = "l" )



# 高速フーリエ変換を行ないます。
data_freqA <- fft( data_waveA - mean( data_waveA ) )
data_freqB <- fft( data_waveB - mean( data_waveB ) )
data_freqC <- fft( data_waveC - mean( data_waveC ) )

# 変換結果を表示します。
# 合成波の周波数も表示されています。
par( "mfrow" = c( 3, 1 ) )
plot( hz, abs( data_freqA[ hz ] ), xlim = c(0, 400), type="l" )
plot( hz, abs( data_freqB[ hz ] ), xlim = c(0, 400), type="l" )
plot( hz, abs( data_freqC[ hz ] ), xlim = c(0, 400), type="l" )

# plot( hz, abs( data_freqC[ hz ] ), xlim = c(0, 400), type="l" )

# グラフをファイルに保存する
dirPath <- "C:/var/work/data/fft"
#pngFileName <- "fft_p_pv_d20.png"
pngFileName <- paste("fft_p_pv_d20_", c("A", "B", "C"), ".png", sep = "")
pngFileName <- paste("fft_p_pv_d365_", c("A", "B", "C"), ".png", sep = "")
pngFilePath <- paste(dirPath, pngFileName, sep = "/")

# 描画デバイスを開く
png(pngFilePath[1], width = 480, height = 270)  
plot( hz, abs( data_freqA[ hz ] ), xlim = c(0, 400), type="l" )
# 描画デバイスを閉じる
dev.off() 

# 描画デバイスを開く
png(pngFilePath[2], width = 480, height = 270)  
plot( hz, abs( data_freqB[ hz ] ), xlim = c(0, 400), type="l" )
# 描画デバイスを閉じる
dev.off() 

# 描画デバイスを開く
png(pngFilePath[3], width = 480, height = 270)  
plot( hz, abs( data_freqC[ hz ] ), xlim = c(0, 400), type="l" )
# 描画デバイスを閉じる
dev.off() 

