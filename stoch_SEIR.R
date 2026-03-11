library(tidyverse)

# Set up ----

## Population ----

N0 <- 1e4L # Initial population size
E0 <- 5L   # Initial number of exposed
X0 <- list(
    S = N0 - E0, # Susceptible
    E = E0,      # Exposed
    I = 0L,      # Infectious
    R = 0L       # Exposed
)

# Population is sampled at these times
times <- seq(0, 100, by = 0.5)

## Model parameters ----
pars <- list(
    beta = 2.0,     # infection coefficient
    eta = 1/14,     # incubation period (~ 14 days)
    gamma = 1/7,    # recovery rate (~ 7 days)
    N0 = N0,        # population equilibrium size
    mu = 0.1 / 365, # growth and mortality rate (~ 10 years)
    nu = 0.1 / 14,  # disease induced mortality (~ 10% case fatality)
    c = 2,          # culling rate of infectious (2 / day)
    # times of culling (in pairs)
    culls = c(20, 22,
              30, 32,
              40, 42,
              50, 52)
)

# Main function ----

stoch_seir <- function(X0, times, pars) {
    # Stochastic SEIR model with demography and time dependent parameter
    # using Gillespie SSA algorithm
    
    compartments <- names(X0)
    
    # Storing all the population data in a matrix is very efficient, can convert
    # to data.frame or tibble later
    X <- matrix(NA_integer_,
                nrow = length(times),
                ncol = length(compartments),
                dimnames = list(NULL, compartments))
    
    ## Matrix of events ----
    event_mat <- matrix(
        #  S   E   I   R
        c(+1, +0, +0, +0, # birth
          -1, +0, +0, +0, # death_S
          +0, -1, +0, +0, # death_E
          +0, +0, -1, +0, # death_I
          +0, +0, +0, -1, # death_R
          -1, +1, +0, +0, # infection S -> E
          +0, -1, +1, +0, # incubation E -> I
          +0, +0, -1, +1  # recovery I -> R
        ) |> as.integer(),
        byrow = TRUE,
        ncol = length(compartments),
        dimnames = list(c("birth", "death_S", "death_E", "death_I", "death_R",
                          "infection", "incubation", "recovery"),
                        compartments))
    
    # Meta variables
    x <- unlist(X0)
    t <- first(times)
    t_end <- last(times)
    
    X[1L, compartments] <- x
    rec <- 2L
    
    # Main loop ----
    
    while (TRUE) {
        ## Calculate event rates ----
        event_rates <- with(c(as.list(x), pars), {
            N <- sum(x)
            
            # Culling of infectives at c/day when t in a culling interval
            c1 <- t >= culls[c(TRUE, FALSE)]
            c2 <- t <  culls[c(FALSE, TRUE)]
            ct <- if (any(c1 & c2)) c else 0
            
            c(birth   = mu * N * (1 - N / N0),
              death_S = mu * S,
              death_E = mu * E,
              death_I = (mu + nu + ct) * I,
              death_R = mu * R,
              # Mind the divide by 0 case when N == 0
              infection  = if (N == 0) 0 else beta * S * I / N,
              incubation = eta * E,
              recovery   = gamma * I)
        })
        
        ## Check for negative event rates ----
        if (any(event_rates < 0)) {
            print(event_rates)
            stop("Negative event detected!")
        }
        
        ## Calculate dx/dt ----
        # In case we implement the tau leaping algorithm
        # dX <- with(as.list(event_rates),
        #            c(dS = + birth - death_S - infection,
        #              dE = - death_E + infection - incubation,
        #              dI = - death_I + incubation - recovery,
        #              dR = - death_R + recovery))
        
        ## Calculate dt and select random event ----
        total_event_rate <- sum(event_rates)
        dt <- if (total_event_rate > 0) {
            rexp(1, rate = total_event_rate)
        } else {
            t_end
        }
        event <- sample(names(event_rates), 1L,
                        prob = event_rates)
        
        ## Update model ----
        x <- x + event_mat[event, ]
        t <- t + dt
        
        ## Record data ----
        while (rec <= length(times) && t > times[rec]) {
            X[rec, ] <- x
            rec <- rec + 1L
        }
        
        # Test for end of simulation
        if (t > t_end) break
    }
    
    # Combine times and compartments into data.frame (remove # for tibble)
    data.frame(time = times, X) # |> tibble()
}

# Run simulation ----
X <- stoch_seir(X0, times, pars)

# In case we want to add N
# X <- bind_cols(X, X |> rowwise() |> summarise(N = sum(c_across(is.integer))))


# Plot ----

# Need long format for ggplot
X1 <- pivot_longer(X, cols = -time, names_to = "compartment")

plt <- ggplot(X1,
              aes(time, value, colour = compartment)) +
    geom_line() +
    geom_vline(xintercept = pars$culls,
               linetype = "dashed",
               linewidth = 0.2) +
    scale_colour_manual("Compartment",
                        breaks = c("S", "E", "I", "R"),
                        values = c("blue", "red4", "red", "green4")) +
    labs(x = "Time (days)",
         y = "Popn",
         colour = "Compartment",
         title = "Stochastic SEIR model") +
    theme_classic()

print(plt)
