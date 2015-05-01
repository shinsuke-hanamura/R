# �K�E�X���z�ɑ΂���x�C�Y���_


#���肷��K�E�X���z
g_mean <- 0.8
g_var <- 0.1
g_sd <- sqrt(g_var)
# ���x ��:1/��^2
g_lambda <- 1/ g_var

#
x_range <- c(-100, 200) * 0.01
g_x <- c(-100:200) * 0.01
g_p <- dnorm(g_x, g_mean, g_sd)

par(new=F)
plot(g_x, g_p, col = "red", xlim=x_range, xlab="" , ylab = "", type="l", ann = F, axes = F)
axis(side=4)



## 0. ���O���z�̐ݒ�
g_sd_0 <- sqrt(0.1)
g_mean_0 <- 0

## �x�C�Y�X�V
ite <- 1000
e_mean <- rep(0, ite)
mean_ml <- rep(0, ite)

for (i in c(1:ite)) {

    # �K�E�X���z�ɏ]���f�[�^�𐶐�����
    n_size <- 100
    g_y <- rnorm(n_size, g_mean, g_sd)

    #par(new=T)
    #hist(g_y, xlab="" , ylab = "freq of data(hist)", main = "density of distribution(line)")

    # �������ꂽ�f�[�^����A���z���x�C�Y���肷��


    ## 1. ���U�����m�Ƃ��A���ς𐄒肷��
    mean_ml[i] <- mean(g_y)
    e_mean[i] <- ((g_sd^2 * g_mean_0) + (n_size * g_sd_0^2 * mean_ml[i])) / (n_size * g_sd_0^2 + g_sd^2)
    e_sd <- 1 / sqrt((1 / g_sd_0^2) + (n_size / g_sd^2))

    sprintf("e_mean[%d] %3.3f", i, e_mean)
    sprintf("e_sd %3.3f", e_sd)

    #par(new=T)
    g_x <- c(-100:200) * 0.01
    g_e_p <- dnorm(g_x, e_mean[i], g_sd)
    #plot(g_x, g_e_p, xlim=x_range, xlab="" , ylab = "", type="l", ann = F, axes = F)

    ## 2. ���O���z�̍X�V
    g_mean_0 <- e_mean[i]

}
# �Ŗސ���l �� �x�C�X����l�̔�r

comp <- abs(mean_ml - g_mean) - abs(e_mean - g_mean)
# �Ŗސ���l ��� �x�C�X����l �̕����^�̒l�ɋ߂��ꍇ�A�ԐF�\��
c_comp <- (comp > 0) + 1
plot(comp, col=c_comp)

# �Ŗސ���l ��� �x�C�X����l �̕����^�̒l�ɋ߂��ꍇ�ATRUE(comp > 0)
table(comp > 0)

plot(e_mean, mean_ml, col=c_comp)

# �O���t�`��
par(new=F)
plot(g_x, g_p, xlim=x_range, col = "green", xlab="" , ylab = "", type="l", ann = F, axes = F)
axis(side=4)
par(new=T)
hist(g_y, xlim=x_range, xlab="" , ylab = "freq of data(hist)", main = "density of distribution(line)")
par(new=T)
plot(g_x, xlim=x_range, g_e_p, xlab="" , ylab = "", type="l", ann = F, axes = F)











d <- c(1 : 365)

# ���̌W��
chi <- function (d) {
  return((2 * pi * (d - 1) ) / 365)
}

# �ώ���
e_t <- function (chi) {
  return((0.0172 + 0.4281 * cos(chi) - 7.3515 * sin(chi) - 3.3495 * cos(2 * chi) - 9.3619 * sin(2 * chi) ) / 60)

}

# d <- c(1 : 365)

e_t_d <- e_t(chi(d))

plot(e_t_d, type="l")

# e_t_d_fft <- fft( e_t_d - mean( e_t_d ) )
# hz <- c(1:length(e_t_d_fft))
# plot( hz, abs(e_t_d_fft), type="l" )

# �o�x
lambda <- 135

# �ܓx
phi <- (pi / 180) * 35

# ���p
omega <- function(h, lambda, e_t) {
  return(15 * h + lambda - 135 + 15 * e_t)
}

h_noon <- 12
omega_d <- omega(h_noon, lambda, e_t_d)

# ���z�Ԉ�
delta <- function(chi) {
  return(0.006918 - 0.399912 * cos(chi) + 0.070257 * sin(chi) - 0.006758 * cos(2 * chi) + 0.000908 * sin(2 * chi))
}

delta_d <- delta(chi(d))


# ���z���x
s_h <- function(phi, delta, omega) {
  return(asin(cos(phi) * cos(delta) * cos((pi / 180) * omega) + sin(delta) * sin(phi)))
}

s_h_d <- s_h(phi, delta_d, omega_d)


# ���z���ʊp
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

# �o�x
lambda <-135
# �ܓx
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



# ���p��ω�������

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
  # ���z���x
  s_h_d <- s_h(phi, delta_d, omega_d)
  # ���z���ʊp
  psi_d <- psi(delta_d, omega_d, s_h_d)
  
  x <- sin(psi_d)
  y <- cos(psi_d)
  z <- sin(s_h_d)
  c <- (s_h_d > 0) + 1
  titles <- paste("Solar Altitude ", d, "th day", sep = "")
  
  
  # �`��f�o�C�X���J��
  pngFileName <- paste("solar_", fillzero(3, d), ".png", sep = "")
  pngFilePath <- paste(dirPath, pngFileName, sep = "/")

  png(pngFilePath, width = 360, height = 360)  
  scatterplot3d(x, y, z,
                xlim = c(-2, 2), ylim = c(-2, 2), zlim = c(-2, 2),
                color = c,
                col.axis="blue", col.grid="lightblue",
                main = titles, type="h", pch=20) 
  # �`��f�o�C�X�����
  dev.off() 

  # �`��f�o�C�X���J��
  pngFileName <- paste("solar_", fillzero(3, d), "_2d.png", sep = "")
  pngFilePath <- paste(dirPath2d, pngFileName, sep = "/")

  titles <- paste("2d Solar Altitude ", d, "th day", sep = "")
  png(pngFilePath, width = 360, height = 360)  
  plot(x, y,
                xlim = c(-2, 2), ylim = c(-2, 2), 
                col = c,
                main = titles, pch=20) 
  # �`��f�o�C�X�����
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



# Daily 360Hz�̐����g
data_waveA <- sin( 2 * pi * t * 360 )

# Anually 1Hz�̐����g
data_waveB <- sin( 2 * pi * t * 1 )

# ��L2�̍����g
data_waveC <- data_waveA + data_waveB

# �g��`���܂��B
par( "mfrow" = c( 3, 1 ) )
plot( t, data_waveA, type="l" )
plot( t, data_waveB, type="l" )
plot( t, data_waveC, type="l" )

# plot( t, data_waveA, xlim = c(0, 2/365), type = "l" )
# plot( t, data_waveC, xlim = c(0, 22/365), type = "l" )



# �����t�[���G�ϊ����s�Ȃ��܂��B
data_freqA <- fft( data_waveA - mean( data_waveA ) )
data_freqB <- fft( data_waveB - mean( data_waveB ) )
data_freqC <- fft( data_waveC - mean( data_waveC ) )

# �ϊ����ʂ�\�����܂��B
# �����g�̎��g�����\������Ă��܂��B
par( "mfrow" = c( 3, 1 ) )
plot( hz, abs( data_freqA[ hz ] ), xlim = c(0, 400), type="l" )
plot( hz, abs( data_freqB[ hz ] ), xlim = c(0, 400), type="l" )
plot( hz, abs( data_freqC[ hz ] ), xlim = c(0, 400), type="l" )

# plot( hz, abs( data_freqC[ hz ] ), xlim = c(0, 400), type="l" )

# �O���t���t�@�C���ɕۑ�����
dirPath <- "C:/var/work/data/fft"
#pngFileName <- "fft_p_pv_d20.png"
pngFileName <- paste("fft_p_pv_d20_", c("A", "B", "C"), ".png", sep = "")
pngFileName <- paste("fft_p_pv_d365_", c("A", "B", "C"), ".png", sep = "")
pngFilePath <- paste(dirPath, pngFileName, sep = "/")

# �`��f�o�C�X���J��
png(pngFilePath[1], width = 480, height = 270)  
plot( hz, abs( data_freqA[ hz ] ), xlim = c(0, 400), type="l" )
# �`��f�o�C�X�����
dev.off() 

# �`��f�o�C�X���J��
png(pngFilePath[2], width = 480, height = 270)  
plot( hz, abs( data_freqB[ hz ] ), xlim = c(0, 400), type="l" )
# �`��f�o�C�X�����
dev.off() 

# �`��f�o�C�X���J��
png(pngFilePath[3], width = 480, height = 270)  
plot( hz, abs( data_freqC[ hz ] ), xlim = c(0, 400), type="l" )
# �`��f�o�C�X�����
dev.off() 

