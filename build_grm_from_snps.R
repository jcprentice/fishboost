# https://zjuwhw.github.io/2021/08/20/GRM.html

# Simulate the genotype matrix

n <- 5 # individuals
m <- 3 # no. of SNPs

set.seed(10)

# Allele frequency are drawn from a uniform distribution
p <- runif(m, min = 0.2, max = 0.5)

# For the plink ped file
x_A1 <- matrix(rbinom(n * m, 1, p), n, m, byrow = TRUE)
x_A2 <- matrix(rbinom(n * m, 1, p), n, m, byrow = TRUE)

x <- x_A1 + x_A2
rownames(x) <- str_c("indi", 1:n)
colnames(x) <- str_c("SNP", 1:m)
x

# Now we save the matrix x, convert to plink bfile, and then calculate GRM using GCTA
x_A1 <- ifelse(x_A1 == 0, "a", "A")
x_A2 <- ifelse(x_A2 == 0, "a", "A")

if (FALSE) {
    message(" - mkdir toy")
    dir.create("toy")
    write.table(cbind(data.frame(FID = rownames(x), IID = rownames(x), 0, 0, 0, 0),
                      cbind(x_A1, x_A2)[, order(c(1:m, 1:m + 0.5))]),
                "toy/toy.ped",
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(data.frame(CHR = 1, SNP = colnames(x), cM = 0, BP = 10000 + 1:ncol(x)),
                "toy/toy.map",
                sep = "\t", quote = F, row.names = F, col.names = F)
}

# After standardization
p_hat <- apply(x, 2, sum) / (2 * n)
w <- apply(rbind(x, p_hat), 2, \(x) {
    xl <- x[[length(x)]]
    (x - 2 * xl) / sqrt(2 * xl * (1 - xl))
})[1:n, ]
w

# Calculate the GRM
A <- w %*% t(w) / m
A
diag(A)
A[row(A) < col(A)]
A[row(A) > col(A)]
