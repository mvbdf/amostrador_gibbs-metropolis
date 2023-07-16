using Random
using Distributions
using Plots

Random.seed!(123)

# Razão do núcleo da condicional completa de k nos pontos do observação
# proposta e do último valor aceito
function rk(y, k, lambda1, lambda2, X)
    lambda1^(sum(X[1:y]) - sum(X[1:k])) *
      lambda2^(sum(X[(y+1):112]) - sum(X[(k+1):112])) *
      exp((k-y)*(lambda1-lambda2))
end

# Amostrador
# shape1: parâmetro de forma das distribuições a priori de lambda1 e lambda2
# shape2: parâmetro de forma da distribuição a priori de alpha
function amostrador(shape1, shape2, X, R)
    lambda1 = zeros(Float64, R)
    lambda2 = zeros(Float64, R)
    alpha   = zeros(Float64, R)
    k       = zeros(Int, R)

    # Valores iniciais
    lambda1[1] = 1
    lambda2[1] = 5
    alpha[1]   = 5
    k[1]       = 30

    nac = 0 # número de sucessos na geração de k

    for j in 2:R
        lambda1[j] = rand(Gamma(shape1 + sum(X[1:k[j-1]]), 1/(k[j - 1] + alpha[j - 1])), 1)[1]
        lambda2[j] = rand(Gamma(shape1 + sum(X[(k[j-1]+1):112]), 1/(112 - k[j - 1] + alpha[j - 1])), 1)[1]
        alpha[j]   = rand(Gamma(shape2 + 6, 1/(lambda1[j] + lambda2[j] + 10)), 1)[1]
    
        # Passo de Metropolis-Hastings para k
        y = rand(DiscreteUniform(1, 111), 1)[1] # gera o candidato

        if rand(Uniform(0,1), 1)[1] > 1 - minimum([1, rk(y, k[j-1], lambda1[j], lambda2[j], X)])
            k[j] = y
            nac += 1
        else
            k[j] = k[j - 1]
        end
    end
    
    println("Taxa de aceitação de k: ", nac/R)

    return [lambda1, lambda2, alpha, k]
end

function main()
    X = [4, 5, 4, 1, 0, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6, 3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5,
         3, 4, 2, 5, 2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 1, 1, 3, 0, 0, 1, 0, 1, 1, 0, 0, 3, 1,
         0, 3, 2, 2, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2, 3, 3, 1, 1,
         2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 0, 1, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1]
    R = 10000

    (lambda1, lambda2, alpha, k) = amostrador(3, 10, X, R)
    println("Estimativas: ", mean.([lambda1, lambda2, alpha, k]))

    # Gráficos das cadeias geradas
    plot([lambda1 lambda2 alpha k], xlabel = "Índice", ylabel = ["λ₁" "λ₂" "α" "k"], labels = "", layout = 4)
    savefig("plot1.pdf")

    # Distribuições a posteriori
    histogram([lambda1 lambda2 k], layout = (1,3), xlabel = ["λ₁" "λ₂" "k"], bins = :sturges, labels = "")
    savefig("plot2.pdf")

    # Médias ergódicas
    plot([(cumsum(lambda1) ./ (1:R)) (cumsum(lambda2) ./ (1:R)) (cumsum(alpha) ./ (1:R)) (cumsum(k) ./ (1:R))],
        layout = 4, labels = "", ylabel = ["λ₁" "λ₂" "α" "k"])
    savefig("plot3.pdf")

    # Estudo da sensibilidade em relação à escolha dos parâmetros
    (lambda1, lambda2, alpha, k) = amostrador(30, 100, X, R)
    println("Estimativas: ", mean.([lambda1, lambda2, alpha, k]))

    # Gráficos das cadeias geradas
    plot([lambda1 lambda2 alpha k], xlabel = "Índice", ylabel = ["λ₁" "λ₂" "α" "k"], labels = "", layout = 4)
    savefig("plot4.pdf")
end

main()