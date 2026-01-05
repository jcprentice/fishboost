# Base set of libraries needed to run pipeline without plotting

{
    library(data.table) # Used for all data table manipulations
    library(purrr)      # Functional programming is great
    library(stringr)    # Better string manipulation
    library(MASS)       # Needed for mvnorm
    library(tictoc)     # Check performance
    library(lubridate)  # Tools for dates
}

# Useful but not loaded by default

# library(here)       # Useful according to Hadley Wickham
# library(xml2)       # Used to generate XML file
# library(tidyverse)  # Modern way of working with tables
# library(magrittr)   # More pipes for FP
# library(ggplot2)    # Plotting
# library(ggthemes)   # Plotting
# library(gganimate)  # Plotting
# library(cowplot)    # Grid plotting
# library(pROC)       # Working with ROC curves.
# library(codetools)  # Check for global variables

# Other libraries other people have found useful:

# corpcor           # efficient estimation of Cov and (partial) Corr
# MCMCglmm

# In case any are missing
install_base <- function() {
    install.packages(c("data.table", "tidyverse", "Matrix", "MASS", "coda",
                       "xml2", "tictoc", "gtools", "codetools", "janitor",
                       "fitdistrplus", "lubridate", "ggplot2", "ggthemes",
                       "gganimate", "ggcorrplot", "cowplot", "pROC",
                       "microbenchmark"))
}
