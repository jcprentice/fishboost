{
    source("libraries.R")
    source("source_files.R")
    source("crop.R")
    source("set_cors.R")
    library(purrr)
    library(ggplot2)
}


## Set parameters ----
{
    params <- make_parameters(
        model_type = "SEIDR", # "SIR", "SEIR", "SIDR", or "SEIDR"
        name = "expt1",
        setup = "single", # chris, small, fishboost, fb1, fb2, single
        use_traits = "sit", # "all", "none", "sit", "si" etc.
        vars = 0.5, # c(0.5, 1.2, 0.5)
        cors = 0.2, # 0.2, "positive", "negative", "mixed", "sit_only" etc.
        group_layout = "family", # "random", "family", "striped", "fishboost"
        trial_fe = "ildt",
        donor_fe = "ildt",
        txd_fe = "ildt",
        weight_fe = "sildt",
        weight_is_nested = TRUE,
        sim_new_data = "r"
    )
    
    params$traits_source <- "pedigree"
    params$patch_dataset <- "fb-final"
    params$patch_name <- "scen-1-1"
    params$patch_state <- TRUE
    params$msgs <- FALSE
}


params_base <- copy(params)

# Make sure we keep the right number of initial infectives, in case some are
# removed during the cropping stage
reset_donors <- function(popn, params) {
    I0 <- params$I0
    popn2 <- copy(popn)
    popn2[sdp == "progeny", donor := sample(c(rep(1, I0), rep(0, .N - I0)))]
    popn2
}


N <- 10L

R0s    <- numeric(N)
R0s_s  <- numeric(N)
R0s_r0 <- numeric(N)

prop <- 0.8


{
    tic()
    for (i in seq_len(N)) {
        message("\n", str_glue("Pop {i} / {N}..."), "\n")
        
        params <- params_base |>
            patch_params() |>
            fix_weight_fes() |>
            apply_links()
        
        popn <- make_pedigree(params) |>
            set_traits(params) |>
            set_groups(params) |>
            apply_fixed_effects(params)
        
        # Add in R0
        popn[, r0 := sus_g + inf_g + tol_g]
        
        # crop, set groups, and apply FEs
        popn_s  <- crop(popn, trait = "sus_g", prop = prop) |>
            reset_donors(params)
        popn_r0 <- crop(popn, trait = "r0", prop = prop) |>
            reset_donors(params)
        popn <- popn |>
            reset_donors(params)
        
        R0s[[i]]    <- simulate_epidemic(popn,    params) |> get_R0()
        R0s_s[[i]]  <- simulate_epidemic(popn_s,  params) |> get_R0()
        R0s_r0[[i]] <- simulate_epidemic(popn_r0, params) |> get_R0()
    }
    time_taken <- toc()
}


summary(R0s)
summary(R0s_s)
summary(R0s_r0)

R0dt = data.table(R0s, R0s_s, R0s_r0) |>
    melt(measure.vars = 1:3)

ylabs <- c("None", "Susceptibility", "R0")
setattr(R0dt$variable, "levels", ylabs)


plt <- ggplot(R0dt) +
    # geom_density(aes(x = value, fill = variable),
    #              alpha = 0.5) +
    # xlim(0, 10) +
    geom_boxplot(aes(x = value, y = variable, fill = variable),
                 staplewidth = 0.5,
                 outliers = FALSE) +
    geom_vline(xintercept = 1,
               linetype = "dashed") +
    scale_x_continuous(bquote(R[0]),
                       breaks = 0:20) +
    scale_y_discrete(labels = ylabs) +
    guides(fill = guide_legend(position = "inside")) +
    labs(x = bquote(R[0]),
         y = "Selection on",
         title = bquote("Distribution of" ~ R[0] ~ "values")) +
    theme(legend.position = "none")

print(plt)


ggsave(str_glue("{params$patch_dataset}/gfx/{ds}-R0-sim.pdf"),
       plot = plt, width = 6, height = 4)


saveRDS(mget(c("R0s", "R0s_s", "R0s_r0", "time_taken", "plt")),
        file = "r0selection.rds")
# readRDS("r0selection.rds")
# plt <- plot_model(pops[[1]], params)
# print(plt)

