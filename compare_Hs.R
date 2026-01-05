library(data.table)
library(ggplot2)
library(Matrix)
library(cowplot)
library(viridisLite)

# Manually step through import_Hinv.R and import_Hinv2.R to get these matrices
Ricardo   <- fread("Hinv_matrices/ricardo.tsv")
Christina <- fread("Hinv_matrices/christina.tsv")
Pedigree  <- fread("Hinv_matrices/pedigree.tsv")

Ricardo <- rbind(Ricardo, Ricardo[i != j, .(i = j, j = i, x)])
setorder(Ricardo, i, j)

{
    plt1 <- ggplot(Ricardo) +
        geom_tile(aes(x = i, y = j, fill = x)) +
        labs(title = "Ricardo") +
        theme_classic()
    plt2 <- ggplot(Christina) +
        geom_tile(aes(x = i, y = j, fill = x)) +
        labs(title = "Christina") +
        theme_classic()
    plt3 <- ggplot(Pedigree) +
        geom_tile(aes(x = i, y = j, fill = x)) +
        labs(title = "Pedigree") +
        theme_classic()
    heatmap <- plot_grid(plt1, plt2, plt3, nrow = 1)
    ggsave("heatmaps.png", heatmap, width = 12, height = 3)
}

get_inverse <- function(x_dt) {
    x_sm <- with(x_dt, sparseMatrix(i, j, x = x))
    ix_sm <- solve(x_sm)
    ix_dt <- as.data.table(summary(ix_sm))
    #setnames(ix_dt, c("i", "j", "x"))
    ix_dt
}

iR <- get_inverse(Ricardo)
iC <- get_inverse(Christina)
iP <- get_inverse(Pedigree)

{
    plt1 <- ggplot(iR) +
        geom_tile(aes(x = i, y = j, fill = x)) +
        labs(title = "Ricardo")
    plt2 <- ggplot(iC) +
        geom_tile(aes(x = i, y = j, fill = x)) +
        labs(title = "Christina")
    plt3 <- ggplot(iP) +
        geom_tile(aes(x = i, y = j, fill = x)) +
        labs(title = "Christina")
    inverses <- plot_grid(plt1, plt2, plt3, nrow = 1)
    ggsave("heatmaps.png", inverses, width = 12, height = 6)
}


# Ricardo <- fread("fb_data/Hinv.txt", skip = 1)
# setnames(Ricardo, c("i", "j", "x"))
# 
# Christina <- fread("HinvOrig.txt")
# setnames(Christina, c("i", "j", "x"))

dt <- rbind(
    data.table(name = "Ricardo",   y = sort(Ricardo$x)),
    data.table(name = "Christina", y = sort(Christina$x)),
    data.table(name = "Pedigree",  y = sort(Pedigree$x))
)

{
    dt[, x := seq(0, 1, length = .N), name]
    dt[, section := factor(fcase(x < 0.06, "low",
                                 0.984 < x, "high",
                                 default = "mid"),
                           levels = c("low", "mid", "high"))]
    
    ht <- ggplot(dt[section != "mid"]) +
        geom_line(aes(x = x, y = y, colour = name),
                  linewidth = 1.2) +
        facet_wrap(. ~ section, scales = "free")
    ggsave("head_tail.png", ht, width = 12, height = 6)
}
