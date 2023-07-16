import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(123)

# Razão do núcleo da condicional completa de k nos pontos do observação
# proposta e do último valor aceito
def rk(y, k, lambda1, lambda2, X):
    return lambda1**(sum(X[0:y]) - sum(X[0:k])) *\
                lambda2**(sum(X[y:112]) - sum(X[k:112])) *\
                np.exp((k-y)*(lambda1-lambda2))

# Amostrador
# shape1: parâmetro de forma das distribuições a priori de lambda1 e lambda2
# shape2: parâmetro de forma da distribuição a priori de alpha
def amostrador(shape1, shape2, X, R):
    lambda1 = np.zeros(R)
    lambda2 = np.zeros(R)
    alpha   = np.zeros(R)
    k       = np.zeros(R, dtype = int)

    # Valores iniciais
    lambda1[0] = 1
    lambda2[0] = 5
    alpha[0]   = 5
    k[0]       = 30

    nac = 0 # número de sucessos na geração de k

    for i in range(1, R):
        lambda1[i] = np.random.gamma(shape1 + sum(X[0:k[i-1]]), 1/(k[i - 1] + alpha[i - 1]))
        lambda2[i] = np.random.gamma(shape1 + sum(X[k[i-1]:112]), 1/(112 - k[i - 1] + alpha[i - 1]))
        alpha[i]   = np.random.gamma(shape2 + 6, 1/(lambda1[i] + lambda2[i] + 10))
    
        # Passo de Metropolis-Hastings para k
        y = np.random.randint(1, 111)  # gera o candidato

        if np.random.uniform() > 1 - min([1, rk(y, k[i-1], lambda1[i], lambda2[i], X)]):
            k[i] = y
            nac += 1
        else:
            k[i] = k[i-1]
    
    print("Taxa de aceitação de k: ", nac/R)
    return(pd.DataFrame({'lambda1':lambda1, 'lambda2':lambda2, 'k':k, 'alpha':alpha}))


X = [4, 5, 4, 1, 0, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6, 3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5,
        3, 4, 2, 5, 2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 1, 1, 3, 0, 0, 1, 0, 1, 1, 0, 0, 3, 1,
        0, 3, 2, 2, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2, 3, 3, 1, 1,
        2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 0, 1, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1]
R = 10000

cadeias = amostrador(3, 10, X, R)
print(cadeias.apply(np.mean))

names = [r'$\lambda_1$', r'$\lambda_2$', r'$k$', r'$\alpha$']

# Gráficos das cadeias geradas
for i in range(len(names)):
    plt.subplot(2,2, i+1)
    sns.lineplot(cadeias.iloc[:,i])
    plt.ylabel(names[i])
plt.show()

# Distribuições a posteriori
for i in range(len(names)-1):
    plt.subplot(1,3, i+1)
    sns.histplot(cadeias.iloc[:,i], bins = 'sturges')
    plt.xlabel(names[i])
    plt.ylabel("Frequência")
plt.show()

# Médias ergódicas
for i in range(len(names)):
    plt.subplot(2,2, i+1)
    sns.lineplot(np.cumsum(cadeias.iloc[:,i]) / np.arange(1, R+1))
    plt.ylabel(names[i])
plt.show()

# Estudo da sensibilidade em relação à escolha dos parâmetros
cadeias = amostrador(30, 100, X, R)

for i in range(len(names)):
    plt.subplot(2,2, i+1)
    sns.lineplot(cadeias.iloc[:,i])
    plt.ylabel(names[i])
plt.show()

# Gráficos das cadeias geradas
print(cadeias.apply(np.mean))