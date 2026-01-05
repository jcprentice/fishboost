{
    library(stringr)
    source("flatten_bici_states.R")
}

run_bici_sim <- function(dataset = "fb-test",
                         name = "scen-1-1",
                         bici_cmd = "post-sim",
                         nreps = 0) {
    
    # dataset <- "sim-base-inf"; name <- "scen-5-1"; bici_cmd <- "post-sim"; nreps <- 0
    
    bici_file <- str_glue("datasets/{dataset}/data/{name}.bici")
    
    if (!file.exists(bici_file)) {
        message(str_glue("Missing BICI file for '{name}'"))
        return(NULL)
    }
    
    # Open config file and get no. of reps (cf_reps)
    lines <- readLines(bici_file)
    line_n <- str_which(lines, switch(bici_cmd, "sim" = "^simulation.*number", "^post-sim.*number"))
    parts <- str_split_1(lines[[line_n]], " ")
    part_n <- str_which(parts, "number")
    cf_reps <- str_split_i(parts[[part_n]], "=", 2) |> as.numeric()
    
    if (nreps %notin% c(NA, 0, cf_reps)) {
        message(str_glue("modifying '{name}.bici' to change reps"))
        parts[[part_n]] <- str_glue("number={nreps}")
        line <- str_flatten(parts, " ")
        lines[[line_n]] <- line
        writeLines(lines, bici_file)
    }
    
    cmd <- str_glue(
        "mpirun -n {cores} --output :raw --oversubscribe",
        "../BICI/bici-{platform}", bici_file, bici_cmd,
        .sep = " ",
        platform = Sys.info()[["sysname"]],
        cores = 2 #primes::gcd(parallel::detectCores(), cf_reps)
    )
    message(cmd)
    
    out <- system(cmd)
    if (out != 0) return(NULL)
    
    flatten_bici_states(dataset, name, bici_cmd)
}

# etc <- run_bici_sim(dataset = "fb-test", "scen-1-1", "post-sim", 50)
# etc <- run_bici_sim(dataset = "sim-test2", "scen-1-1", "post-sim", 50)
