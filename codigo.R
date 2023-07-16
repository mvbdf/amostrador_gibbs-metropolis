set.seed(123)

# Razão do núcleo da condicional completa de k nos pontos do observação
# proposta e do último valor aceito
rk = function(y, k, lambda1, lambda2, X) {
  lambda1^(sum(X[1:y]) - sum(X[1:k])) *
    lambda2^(sum(X[(y + 1):112]) - sum(X[(k + 1):112])) *
    exp((k - y) * (lambda1 - lambda2))
}

# Amostrador
# shape1: parâmetro de forma das distribuições a priori de lambda1 e lambda2
# shape2: parâmetro de forma da distribuição a priori de alpha
amostrador = function(shape1, shape2, X, R) {
  lambda1 = lambda2 = alpha = k = numeric(R)
  
  # Valores iniciais
  lambda1[1] = 1
  lambda2[1] = 5
  alpha[1]   = 5
  k[1]       = 30
  
  nac = 0 # número de sucessos na geração de k
  
  for (j in 2:R) {
    lambda1[j] = rgamma(1, shape1 + sum(X[1:k[j-1]]), k[j - 1] + alpha[j - 1])
    lambda2[j] = rgamma(1, shape1 + sum(X[(k[j-1]+1):112]), 112 - k[j - 1] + alpha[j - 1])
    alpha[j]   = rgamma(1, shape2 + 6, lambda1[j] + lambda2[j] + 10)
    
    # Passo de Metropolis-Hastings para k
    y = sample(111, 1) # gera o candidato
    
    if (runif(1) > 1 -  min(1, rk(y, k[j-1], lambda1[j], lambda2[j], X))) {
      k[j] = y
      nac = nac + 1
    } else
      k[j] = k[j - 1]
  }
  
  # Taxa de aceitação de k
  print(paste("Aceitação de k: ", nac / R))
  
  data.frame(lambda1 = lambda1, lambda2 = lambda2, alpha = alpha, k = k)
}

X = c(4, 5, 4, 1, 0, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6, 3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5,
      3, 4, 2, 5, 2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 1, 1, 3, 0, 0, 1, 0, 1, 1, 0, 0, 3, 1,
      0, 3, 2, 2, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2, 3, 3, 1, 1,
      2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 0, 1, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1)
R = 10000

cadeias = amostrador(3, 10, X, R)

# Gráficos das cadeias geradas
par(mfrow = c(2,2))
plot(cadeias$lambda1, type = 'l', xlab = "Índice", ylab = expression(lambda[1]))
plot(cadeias$lambda2, type = 'l', xlab = "Índice", ylab = expression(lambda[2]))
plot(cadeias$alpha, type = 'l', xlab = "Índice", ylab = expression(alpha))
plot(cadeias$k, type = 'l', xlab = "Índice", ylab = expression(k))

# Estimativas
apply(cadeias, 2, mean)

# Distribuições a posteriori
par(mfrow = c(1, 3))
hist(cadeias$lambda1, ylab = "Frequência", xlab = expression(lambda[1]), main = ""); box();
hist(cadeias$lambda2, ylab = "Frequência", xlab = expression(lambda[2]), main = ""); box();
hist(cadeias$k, ylab = "Frequência", main = ""); box()

# Médias ergódicas
par(mfrow = c(2,2))
plot(cumsum(cadeias$lambda1) / (1:R), type = "l", xlab = "Índice",
     ylab = expression(paste("Médias ergódicas de ", lambda[1])))
plot(cumsum(cadeias$lambda2) / (1:R), type = "l", xlab = "Índice",
     ylab = expression(paste("Médias ergódicas de ", lambda[2])))
plot(cumsum(cadeias$alpha) / (1:R), type = "l", xlab = "Índice",
     ylab = expression(paste("Médias ergódicas de ", alpha)))
plot(cumsum(cadeias$k) / (1:R), type = "l", xlab = "Índice", ylab = "Médias ergódicas de k")

# Estudo da sensibilidade em relação à escolha dos parâmetros
cadeias = amostrador(30, 100, X, R)

# Gráficos das cadeias geradas
par(mfrow = c(2,2))
plot(cadeias$lambda1, type = 'l', xlab = "Índice", ylab = expression(lambda[1]))
plot(cadeias$lambda2, type = 'l', xlab = "Índice", ylab = expression(lambda[2]))
plot(cadeias$alpha, type = 'l', xlab = "Índice", ylab = expression(alpha))
plot(cadeias$k, type = 'l', xlab = "Índice", ylab = expression(k))

# Estimativas
apply(cadeias, 2, mean)