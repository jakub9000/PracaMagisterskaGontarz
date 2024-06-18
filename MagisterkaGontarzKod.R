## DANE SYMULOWANE

# biblioteka igraph z wbudowanym stochastycznym modelem blokowym
# install.packages("igraph")
library("igraph")

# biblioteki do operacji na macierzach
# install.packages("expm")
# install.packages("matlib")
library("expm")
library("matlib")

# bibliotek do szybszych wyznaczen rozkladow osobliwych i spektralnych
# install.packages("RSpectra")
library("RSpectra")

# install.packages("gtools")
library("gtools")

# funkcja do algorytmu na laplasjanie skalowanym
laplasjan <- function(P, rho, n, K, d, ilosc_prob = 5){
  wyniki_pods_laplas <- rep(0,length(n))
  
  for (l in 1:ilosc_prob){
    for (k in 1:length(n)){
      
      graph <- sample_sbm(n[k],P,n[k]*rho, directed = F ,loops = F)
      
      bloki <- rep(1,n[k]*rho[1])
      
      for(i in 2:length(rho)){
        
        bloki <- c(bloki, rep(i,n[k]*rho[i]))
      }
      
      V(graph)$color <- bloki
      
      L <- diag(rep(1,n[k])) - laplacian_matrix(graph, normalized = T)
      
      X <- eigs(L, d)[[2]]
      
      wyniki <- kmeans(X,K)[[1]]
      
      minimum <- 1

      wyniki_org <- wyniki

      permutacje <- permutations(K,K)
      
      for (i in 1:nrow(permutacje)){
        wyniki <- wyniki_org
        for(j in 1:K){
          wyniki[wyniki_org==j] <- permutacje[i,j]
        }
        minimum <- min(minimum, 1 - sum(wyniki == V(graph)$color)/n[k])
      }
      
      wyniki_pods_laplas[k] <- wyniki_pods_laplas[k] + minimum
      
    }
  }
  
  # wyniki_pods_laplas_unnorm <- wyniki_pods_laplas_unnorm/ilosc_prob
  wyniki_pods_laplas <- wyniki_pods_laplas/ilosc_prob
  
  return(wyniki_pods_laplas)
}

# funkcja do algorytmu na macierzy sasiedztwa dla grafu nieskierowanego 
sasiedztwa_nieskierowany <- function(P, rho, n, K, d, ilosc_prob = 5){
  wyniki_pods_sasiedz <- rep(0,length(n))
  wyniki_pods_sasiedz_skalowane <- rep(0,length(n))
  
  for (l in 1:ilosc_prob){
    for (k in 1:length(n)){
      
      graph <- sample_sbm(n[k],P,n[k]*rho, directed = F ,loops = F)
      
      bloki <- rep(1,n[k]*rho[1])
      
      for(i in 2:length(rho)){
        
        bloki <- c(bloki, rep(i,n[k]*rho[i]))
      }
      
      V(graph)$color <- bloki
      
      A <- as.matrix(as_adjacency_matrix(graph))
      
      svd <- svds(A,d)
      
      sigma <- diag(svd[[1]])
      U <- svd[[2]]
      # V <- svd[[3]]
      
      U_skalowane <- U %*% sqrtm(sigma)
      # V_skalowane <- V %*% sqrtm(sigma)
      
      # Z <- cbind(U_skalowane, V_skalowane)
      
      wyniki <- kmeans(U,K)[[1]]
      wyniki_skalowane <- kmeans(U_skalowane,K)[[1]]
      
      minimum <- 1
      minimum_skalowane <- 1
      
      wyniki_org <- wyniki
      wyniki_skalowane_org <- wyniki_skalowane
      
      permutacje <- permutations(K,K)
      
      for (i in 1:nrow(permutacje)){
        wyniki <- wyniki_org
        wyniki_skalowane <- wyniki_skalowane_org
        for(j in 1:K){
          wyniki[wyniki_org==j] <- permutacje[i,j]
          wyniki_skalowane[wyniki_skalowane_org==j] <- permutacje[i,j]
        }
        minimum <- min(minimum, 1 - sum(wyniki == V(graph)$color)/n[k])
        minimum_skalowane <- min(minimum_skalowane, 1 - sum(wyniki_skalowane == V(graph)$color)/n[k])
      }
      
      wyniki_pods_sasiedz[k] <- wyniki_pods_sasiedz[k] + minimum
      wyniki_pods_sasiedz_skalowane[k] <- wyniki_pods_sasiedz_skalowane[k] + minimum_skalowane
      
    }
  }
  
  wyniki_pods_sasiedz <- wyniki_pods_sasiedz/ilosc_prob
  wyniki_pods_sasiedz_skalowane <- wyniki_pods_sasiedz_skalowane/ilosc_prob
  
  out <- matrix(c(wyniki_pods_sasiedz,wyniki_pods_sasiedz_skalowane), nrow = 2, ncol = length(n), byrow = T)
  return(out)
}

# funkcja do algorytmu na macierzy sasiedztwa dla grafu skierowanego 
sasiedztwa_skierowany <- function(P, rho, n, K, d, ilosc_prob = 5){
  wyniki_pods_sasiedz <- rep(0,length(n))
  wyniki_pods_sasiedz_skalowane <- rep(0,length(n))
  
  for (l in 1:ilosc_prob){
    for (k in 1:length(n)){
      
      graph <- sample_sbm(n[k],P,n[k]*rho, directed = T ,loops = F)
      
      bloki <- rep(1,n[k]*rho[1])
      
      for(i in 2:length(rho)){
        
        bloki <- c(bloki, rep(i,n[k]*rho[i]))
      }
      
      V(graph)$color <- bloki
      
      A <- as.matrix(as_adjacency_matrix(graph))
      
      svd <- svds(A,d)
      
      sigma <- diag(svd[[1]])
      U <- svd[[2]]
      V <- svd[[3]]
      
      U_skalowane <- U %*% sqrtm(sigma)
      V_skalowane <- V %*% sqrtm(sigma)
      
      Z <- cbind(U, V)
      Z_skalowane <- cbind(U_skalowane, V_skalowane)
      
      wyniki <- kmeans(Z,K)[[1]]
      wyniki_skalowane <- kmeans(Z_skalowane,K)[[1]]
      
      minimum <- 1
      minimum_skalowane <- 1
      
      wyniki_org <- wyniki
      wyniki_skalowane_org <- wyniki_skalowane
      
      permutacje <- permutations(K,K)
      
      for (i in 1:nrow(permutacje)){
        wyniki <- wyniki_org
        wyniki_skalowane <- wyniki_skalowane_org
        for(j in 1:K){
          wyniki[wyniki_org==j] <- permutacje[i,j]
          wyniki_skalowane[wyniki_skalowane_org==j] <- permutacje[i,j]
        }
        minimum <- min(minimum, 1 - sum(wyniki == V(graph)$color)/n[k])
        minimum_skalowane <- min(minimum_skalowane, 1 - sum(wyniki_skalowane == V(graph)$color)/n[k])
      }
      
      wyniki_pods_sasiedz[k] <- wyniki_pods_sasiedz[k] + minimum
      wyniki_pods_sasiedz_skalowane[k] <- wyniki_pods_sasiedz_skalowane[k] + minimum_skalowane
      
    }
  }
  
  wyniki_pods_sasiedz <- wyniki_pods_sasiedz/ilosc_prob
  wyniki_pods_sasiedz_skalowane <- wyniki_pods_sasiedz_skalowane/ilosc_prob
  
  out <- matrix(c(wyniki_pods_sasiedz,wyniki_pods_sasiedz_skalowane), nrow = 2, ncol = length(n), byrow = T)
  return(out)
}

# wykres z oryginalnymi wynikami
n <- seq(500,2000, by = 100)
K <- 2
d <- 2

P <- matrix(data = c(0.42,0.42,0.42,0.5),nrow = K, ncol=K, byrow= T)
rho <- matrix(data = c(0.6,0.4), nrow = K)

ilosc_prob <- 100

wyniki_pods_laplas_norm <- laplasjan(P,rho,n,K,d, ilosc_prob)

test <- sasiedztwa_nieskierowany(P, rho, n, K, d, ilosc_prob)
wyniki_pods_sasiedz <- test[1,]
wyniki_pods_sasiedz_skalowane <- test[2,]

plot(n,wyniki_pods_sasiedz, type = "o", col = "blue", pch = 15, xlab = "n - liczba wierzchołków", ylab = "Proporcja błędów", xaxt='n', yaxt='n',
     ylim = c(0,0.5))
axis(1, at = seq(400,2000, by= 200), labels = seq(400,2000, by= 200), las = 0)
axis(2, at = seq(0,0.5, by= 0.05), labels = seq(0,0.5, by= 0.05), las = 0)

lines(n,wyniki_pods_sasiedz_skalowane, type = "o", col = "red", pch = 16)

lines(n,wyniki_pods_laplas_norm, type = "o", col = "black", pch = 15)

legend("topright", legend = c("Nieskalowana sąsiedztwa", "Skalowana sąsiedztwa", "Skalowany laplasjan"), col = c("blue", "red","black"),
       lty = c(1,1,1), cex = 0.8)

# wykres dla wiekszych grafow
n <- seq(2000,5000, by = 500)
K <- 2
d <- 2

P <- matrix(data = c(0.42,0.42,0.42,0.5),nrow = K, ncol=K, byrow= T)
rho <- matrix(data = c(0.6,0.4), nrow = K)

ilosc_prob <- 100

wyniki_pods_laplas_norm <- laplasjan(P,rho,n,K,d, ilosc_prob)

test <- sasiedztwa_nieskierowany(P, rho, n, K, d, ilosc_prob)
wyniki_pods_sasiedz <- test[1,]
wyniki_pods_sasiedz_skalowane <- test[2,]

plot(n,wyniki_pods_sasiedz, type = "o", col = "blue", pch = 15, xlab = "n - liczba wierzchołków", ylab = "Proporcja błędów", xaxt='n', yaxt='n',
     ylim = c(0,0.12))
axis(1, at = seq(2000,5000, by = 500), labels = seq(2000,5000, by = 500), las = 0)
axis(2, at = seq(0,0.12, by= 0.03), labels = seq(0,0.12, by= 0.03), las = 0)

lines(n,wyniki_pods_sasiedz_skalowane, type = "o", col = "red", pch = 16)

lines(n,wyniki_pods_laplas_norm, type = "o", col = "black", pch = 15)

legend("topright", legend = c("Nieskalowana sąsiedztwa", "Skalowana sąsiedztwa", "Skalowany laplasjan"), col = c("blue", "red","black"),
       lty = c(1,1,1), cex = 0.8)

# wykres dla pozostalych macierzy 2x2
n <- seq(500,2000, by = 100)
K <- 2
d <- 2

P <- matrix(data = c(0.15,0.11,0.11,0.13),nrow = K, ncol=K, byrow= T)
P <- matrix(data = c(0.4,0.5,0.5,0.6),nrow = K, ncol=K, byrow= T)
P <- matrix(data = c(0.83,0.82,0.82,0.9),nrow = K, ncol=K, byrow= T)
rho <- matrix(data = c(0.6,0.4), nrow = K)

ilosc_prob <- 100

wyniki_pods_laplas_norm <- laplasjan(P,rho,n,K,d, ilosc_prob)

test <- sasiedztwa_nieskierowany(P, rho, n, K, d, ilosc_prob)
wyniki_pods_sasiedz <- test[1,]
wyniki_pods_sasiedz_skalowane <- test[2,]

plot(n,wyniki_pods_sasiedz_skalowane, type = "o", col = "red", pch = 15, xlab = "n - liczba wierzchołków", ylab = "Proporcja błędów", xaxt='n', yaxt='n',
     ylim = c(0,0.5))
axis(1, at = seq(400,2000, by= 200), labels = seq(400,2000, by= 200), las = 0)
axis(2, at = seq(0,0.5, by= 0.05), labels = seq(0,0.5, by= 0.05), las = 0)

lines(n,wyniki_pods_laplas_norm, type = "o", col = "black", pch = 15)

legend("topright", legend = c("Skalowana sąsiedztwa", "Skalowany laplasjan"), col = c("red","black"),
       lty = c(1,1), cex = 0.8)

#wykres dla macierzy 3x3
n <- seq(500,2000, by = 100)
K <- 3
d <- 2

P <- matrix(data = c(0.4, 0.3 , 0.35,
                     0.3, 0.35, 0.27,
                     0.35, 0.27, 0.5),nrow = K, ncol=K, byrow= T)

P <- matrix(data = c(0.4, 0.35 , 0.35,
                     0.35, 0.39, 0.35,
                     0.35, 0.35, 0.42),nrow = K, ncol=K, byrow= T)

rho <- matrix(data = c(0.3,0.4,0.3), nrow = K)

ilosc_prob <- 100

wyniki_pods_laplas_norm <- laplasjan(P,rho,n,K,d, ilosc_prob)

test <- sasiedztwa_nieskierowany(P, rho, n, K, d, ilosc_prob)
wyniki_pods_sasiedz <- test[1,]
wyniki_pods_sasiedz_skalowane <- test[2,]

plot(n,wyniki_pods_sasiedz_skalowane, type = "o", col = "red", pch = 15, xlab = "n - liczba wierzchołków", ylab = "Proporcja błędów", xaxt='n', yaxt='n',
     ylim = c(0,0.5))
axis(1, at = seq(400,2000, by= 200), labels = seq(400,2000, by= 200), las = 0)
axis(2, at = seq(0,0.5, by= 0.05), labels = seq(0,0.5, by= 0.05), las = 0)

lines(n,wyniki_pods_laplas_norm, type = "o", col = "black", pch = 15)

legend("topright", legend = c("Nieskalowana sąsiedztwa", "Skalowana sąsiedztwa", "Skalowany laplasjan"), col = c("red","black"),
       lty = c(1,1), cex = 0.8)

#wykres dla macierzy 4x4
n <- seq(500,3000, by = 100)
K <- 4
d <- 2

P <- matrix(data = c(0.4, 0.35 , 0.35, 0.35,
                     0.35, 0.39, 0.35, 0.35,
                     0.35, 0.35, 0.42, 0.35,
                     0.35, 0.35, 0.35, 0.41),nrow = K, ncol=K, byrow= T)

rho <- matrix(data = c(0.3,0.2,0.3,0.2), nrow = K)


ilosc_prob <- 100

wyniki_pods_laplas_norm <- laplasjan(P,rho,n,K,d, ilosc_prob)

test <- sasiedztwa_nieskierowany(P, rho, n, K, d, ilosc_prob)
wyniki_pods_sasiedz <- test[1,]
wyniki_pods_sasiedz_skalowane <- test[2,]

plot(n,wyniki_pods_sasiedz_skalowane, type = "o", col = "red", pch = 15, xlab = "n - liczba wierzchołków", ylab = "Proporcja błędów", xaxt='n', yaxt='n',
     ylim = c(0.4,0.8))
axis(1, at = seq(400,3000, by= 200), labels = seq(400,3000, by= 200), las = 0)
axis(2, at = seq(0.4,0.8, by= 0.04), labels = seq(0.4,0.8, by= 0.04), las = 0)

lines(n,wyniki_pods_laplas_norm, type = "o", col = "black", pch = 15)

legend("topright", legend = c("Skalowana sąsiedztwa", "Skalowany laplasjan"), col = c("red","black"),
       lty = c(1,1), cex = 0.8)

#wykres dla macierzy 5x5
n <- seq(500,3000, by = 100)
K <- 5
d <- 2

P <- matrix(data = c(0.4, 0.35 , 0.35, 0.35, 0.35, 
                     0.35, 0.39, 0.35, 0.35, 0.35,
                     0.35, 0.35, 0.42, 0.35, 0.35,
                     0.35, 0.35, 0.35, 0.41, 0.35,
                     0.35, 0.35, 0.35, 0.35, 0.38),nrow = K, ncol=K, byrow= T)

rho <- matrix(data = c(0.3,0.2,0.2,0.2,0.1), nrow = K)

ilosc_prob <- 100

wyniki_pods_laplas_norm <- laplasjan(P,rho,n,K,d, ilosc_prob)

test <- sasiedztwa_nieskierowany(P, rho, n, K, d, ilosc_prob)
wyniki_pods_sasiedz <- test[1,]
wyniki_pods_sasiedz_skalowane <- test[2,]

plot(n,wyniki_pods_sasiedz_skalowane, type = "o", col = "red", pch = 15, xlab = "n - liczba wierzchołków", ylab = "Proporcja błędów", xaxt='n', yaxt='n',
     ylim = c(0.4,0.8))
axis(1, at = seq(400,3000, by= 200), labels = seq(400,3000, by= 200), las = 0)
axis(2, at = seq(0.4,0.8, by= 0.04), labels = seq(0.4,0.8, by= 0.04), las = 0)

lines(n,wyniki_pods_laplas_norm, type = "o", col = "black", pch = 15)

legend("topright", legend = c("Skalowana sąsiedztwa", "Skalowany laplasjan"), col = c("red","black"),
       lty = c(1,1), cex = 0.8)

# rozne d, K = 4
n <- seq(500,3000, by = 100)

K <- 4
P <- matrix(data = c(0.4, 0.35 , 0.35, 0.35,
                     0.35, 0.39, 0.35, 0.35,
                     0.35, 0.35, 0.42, 0.35,
                     0.35, 0.35, 0.35, 0.41),nrow = K, ncol=K, byrow= T)

rho <- matrix(data = c(0.3,0.2,0.3,0.2), nrow = K)

ilosc_prob <- 100

wyniki_sasiedz_d2 <- sasiedztwa_nieskierowany(P,rho,n,K,2,ilosc_prob)
wyniki_sasiedz_d3 <- sasiedztwa_nieskierowany(P,rho,n,K,3,ilosc_prob)
wyniki_sasiedz_d4 <- sasiedztwa_nieskierowany(P,rho,n,K,4,ilosc_prob)

wyniki_pods_sasiedz_skalowane_d2 <- wyniki_sasiedz_d2[2,]
wyniki_pods_sasiedz_skalowane_d3 <- wyniki_sasiedz_d3[2,]
wyniki_pods_sasiedz_skalowane_d4 <- wyniki_sasiedz_d4[2,]

plot(n,wyniki_pods_sasiedz_skalowane_d2, type = "o", col = "red", pch = 15, xlab = "n - liczba wierzchołków", ylab = "Proporcja błędów", xaxt='n', yaxt='n',
     ylim = c(0.2,0.8))
axis(1, at = seq(400,3000, by= 200), labels = seq(400,3000, by= 200), las = 0)
axis(2, at = seq(0.2,0.8, by= 0.06), labels = seq(0.2,0.8, by= 0.06), las = 0)

lines(n,wyniki_pods_sasiedz_skalowane_d3, type = "o", col = "black", pch = 15)
lines(n,wyniki_pods_sasiedz_skalowane_d4, type = "o", col = "blue", pch = 15)

legend("topright", legend = c("d = 2", "d = 3", "d = 4"), col = c("red","black","blue"),
       lty = c(1,1,1), cex = 0.8)

# rozne d, K = 5
n <- seq(500,3000, by = 100)

K <- 5
P <- matrix(data = c(0.4, 0.35 , 0.35, 0.35, 0.35, 
                     0.35, 0.39, 0.35, 0.35, 0.35,
                     0.35, 0.35, 0.42, 0.35, 0.35,
                     0.35, 0.35, 0.35, 0.41, 0.35,
                     0.35, 0.35, 0.35, 0.35, 0.38),nrow = K, ncol=K, byrow= T)

rho <- matrix(data = c(0.3,0.2,0.2,0.2,0.1), nrow = K)

ilosc_prob <- 100

wyniki_sasiedz_d2 <- sasiedztwa_nieskierowany(P,rho,n,K,2,ilosc_prob)
wyniki_sasiedz_d3 <- sasiedztwa_nieskierowany(P,rho,n,K,3,ilosc_prob)
wyniki_sasiedz_d4 <- sasiedztwa_nieskierowany(P,rho,n,K,4,ilosc_prob)

wyniki_pods_sasiedz_skalowane_d2 <- wyniki_sasiedz_d2[2,]
wyniki_pods_sasiedz_skalowane_d3 <- wyniki_sasiedz_d3[2,]
wyniki_pods_sasiedz_skalowane_d4 <- wyniki_sasiedz_d4[2,]

plot(n,wyniki_pods_sasiedz_skalowane_d2, type = "o", col = "red", pch = 15, xlab = "n - liczba wierzchołków", ylab = "Proporcja błędów", xaxt='n', yaxt='n',
     ylim = c(0.2,0.8))
axis(1, at = seq(400,3000, by= 200), labels = seq(400,3000, by= 200), las = 0)
axis(2, at = seq(0.2,0.8, by= 0.06), labels = seq(0.2,0.8, by= 0.06), las = 0)

lines(n,wyniki_pods_sasiedz_skalowane_d3, type = "o", col = "black", pch = 15)
lines(n,wyniki_pods_sasiedz_skalowane_d4, type = "o", col = "blue", pch = 15)

legend("topright", legend = c("d = 2", "d = 3", "d = 4"), col = c("red","black","blue"),
       lty = c(1,1,1), cex = 0.8)


## grafy skierowany a nieskierowany
n <- seq(500,2000, by = 100)
K <- 2
d <- 2

P <- matrix(data = c(0.42,0.42,0.42,0.5),nrow = K, ncol=K, byrow= T)
P <- matrix(data = c(0.82,0.82,0.82,0.9),nrow = K, ncol=K, byrow= T)
P <- matrix(data = c(0.15,0.11,0.11,0.13),nrow = K, ncol=K, byrow= T)
rho <- matrix(data = c(0.6,0.4), nrow = K)


P <- matrix(data = c(0.4, 0.3 , 0.35,
                     0.3, 0.35, 0.27,
                     0.35, 0.27, 0.5),nrow = K, ncol=K, byrow= T)

P <- matrix(data = c(0.4, 0.35 , 0.35,
                     0.35, 0.39, 0.35,
                     0.35, 0.35, 0.42),nrow = K, ncol=K, byrow= T)
rho <- matrix(data = c(0.3,0.4,0.3), nrow = K)

K <- 4
P <- matrix(data = c(0.4, 0.35 , 0.35, 0.35,
                     0.35, 0.39, 0.35, 0.35,
                     0.35, 0.35, 0.42, 0.35,
                     0.35, 0.35, 0.35, 0.41),nrow = K, ncol=K, byrow= T)

rho <- matrix(data = c(0.3,0.2,0.3,0.2), nrow = K)

ilosc_prob <- 100

oba_sasiedz_nieskierowany <- sasiedztwa_nieskierowany(P,rho,n,K,d,ilosc_prob)
wyniki_sasiedz_nieskierowany <- oba_sasiedz_nieskierowany[1,]
wyniki_sasiedz_nieskierowany_skalowane <- oba_sasiedz_nieskierowany[2,]

oba_sasiedz_skierowany <- sasiedztwa_skierowany(P,rho,n,K,d,ilosc_prob)
wyniki_sasiedz_skierowany <- oba_sasiedz_skierowany[1,]
wyniki_sasiedz_skierowany_skalowane <- oba_sasiedz_skierowany[2,]

dev.new(width=960, height=480, unit="px")
plot(n,wyniki_sasiedz_nieskierowany, type = "o", col = "red", pch = 16, xlab = "n - liczba wierzchołków", ylab = "Proporcja błędów", xaxt='n', yaxt='n',
     ylim = c(0,0.6))
axis(1, at = seq(400,2000, by= 200), labels = seq(400,2000, by= 200), las = 0)
axis(2, at = seq(0,0.6, by= 0.06), labels = seq(0,0.6, by= 0.06), las = 0)

lines(n,wyniki_sasiedz_nieskierowany_skalowane, type = "o", col = "blue", pch = 16)

lines(n,wyniki_sasiedz_skierowany, type = "o", col = "black", pch = 15)
lines(n,wyniki_sasiedz_skierowany_skalowane, type = "o", col = "purple", pch = 15)

legend("topright", legend = c("Nieskierowany nieskalowany", 
                              "Nieskierowany skalowany", 
                              "Skierowny nieskalowany",
                              "Skierowany skalowany"), col = c("red","blue","black","purple"),lty = c(1,1,1,1), cex = 0.7)


## blogi 
## graf skierowany z linkami miedzy blogami ktore odegraly wazna role w wybrorach w stanach w 2004 roku
## wartosci wierzcholkow sa opisane jako przynaleznosc do pratii politycznej jednej z dwoch
blogi <- read.graph("C:\\Users\\user\\Desktop\\studia\\MGR4\\polblogs.gml", format = c("gml"))
vcount(blogi)

V(blogi)$value <- V(blogi)$value + 1 

A <- as.matrix(as_adjacency_matrix(blogi))
isSymmetric.matrix(A)

wyniki_pods_sasiedz <- rep(0,length(n))

#algorytm na macierzy sasiedztwa
d <- 2 #osadzenie w dwoch wymiarach chcemy zrobic
K <- 2 #mamy 2 partie polityczne

svd <- svds(A,d)

sigma <- diag(svd[[1]])
U <- svd[[2]]
V <- svd[[3]]

U_skalowane <- U %*% sqrtm(sigma)
V_skalowane <- V %*% sqrtm(sigma)

Z <- cbind(U_skalowane, V_skalowane)

wyniki <- kmeans(Z,2)[[1]]

min(1 - sum(wyniki == V(blogi)$value)/vcount(blogi),sum(wyniki == V(blogi)$value)/vcount(blogi))
