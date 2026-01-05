#include <cmath.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// Find the no. of individuals capable of infecting susceptibles at some point
// (so includes those currently exposed and not yet infectious).
int get_seidr_infectives(DataFrame pop) {
    NumericVector status = pop["status"];
    
    int n = 0;
    for (auto s : status)
        if (s == 'E' || s == 'I' || s == 'D')
            ++n;
    return n;
}


// A Susceptible individual's future disease trajectory is fixed at the point of
// exposure.
List generate_seidr_path(double epi_time, DataFrame pop, int id, List params) {
    double LP_shape = params["LP_shape"];
    double LP_scale = params["LP_scale"];
    double DP_shape = params["DP_shape"];
    double DP_scale = params["DP_scale"];
    double RP_shape = params["RP_shape"];
    double RP_scale = params["RP_scale"];
    
    NumericVector lat = pop["lat"];
    NumericVector det = pop["det"];
    NumericVector recoverability = pop["recoverability"];
    
    double Tinf = epi_time;
    double Tinc = Tinf + R::rgamma(LP_shape, 1.0 / LP_scale) * lat[id];
    double Tsym = Tinc + R::rgamma(DP_shape, 1.0 / DP_scale) * det[id];
    double Trec = Tsym + R::rgamma(RP_shape, 1.0 / RP_scale) * recoverability[id];
    
    return List::create("E", Tinf, Tinc, Tsym, Trec);
}


// Find the time of the next non-infection event, and the id of the individual
List next_seidr_ni_event(DataFrame pop, double epi_time) {
    double t = inf;
    int id = -1;
    
    NumericVector pop_Tinc = pop["Tinc"];
    NumericVector pop_Tsym = pop["Tsym"];
    NumericVector pop_Trec = pop["Trec"];
    NumericVector pop_id = pop["id"];
    
    for (int i = 0; i < pop.size(); ++i) {
        double tinc = pop_Tinc[i];
        double tsym = pop_Tsym[i];
        double trec = pop_Trec[i];
        int pid = pop_id[i];
        
        if (tinc > epi_time && tinc < t) {
            t = tinc; id = pid;
        } else if (tsym > epi_time && tsym < t) {
            t = tinc; id = pid;
        } else if (trec > epi_time && trec < t) {
            t = tinc; id = pid;
        }
    }
    
    /* as.list(pop[, .(.I,
                    Tinc2 = fifelse(Tinc > epi_time, Tinc, inf),
                    Tsym2 = fifelse(Tsym > epi_time, Tsym, inf),
                    Trec2 = fifelse(Trec > epi_time, Trec, inf))
    ][, .(I, Tmin = pmin(Tinc2, Tsym2, Trec2))
    ][, list(t_next_event = min(Tmin, na.rm = TRUE), id_next_event = which.min(Tmin))])
     */
    return List(t, id);
        
}


double signif(double x, int n) {
    return round(x * pow(10.0, n)) * pow(10.0, -n);
}


// Main model ----
// [[Rcpp::export]]
DataFrame model_SEIDR(DataFrame traits, List params) {
    Rcout << "Simulating an SEIDR epidemic ..." << endl;

    // copy necessary parameters
    double r_beta = params["r_beta"];
    bool DEBUG  = params["DEBUG"];

    // initialise populations ----
    DataFrame pop = init_pop(traits, params);
    NumericVector group_inf 

    // Start epidemic simulation loop ----
    double epi_time = 0.0;
    while (get_seidr_infectives(pop) > 0) {
        if (DEBUG) Rcout << "time = " << signif(epi_time, 5) << endl;

        // Calculate infection rates in each group ----
        // this is the sum of the log infectivitities
        // pop[, group_inf := mean(inf * r_beta * (status == "I")), by = group]
        pop[, group_inf := r_beta * GE * mean(fifelse(status %in% c("I", "D"), inf, 0.0)), by = group]

        // if S, infection at rate beta SI
        pop[, event_rate := fifelse(status == "S", sus * group_inf, 0.0)]

        // id and time of next non-infection event
        int ni_event        = next_seidr_ni_event(pop, epi_time);
        double t_next_event = ni_event["t_next_event"];
        int id_next_event   = ni_event["id_next_event"];

        if (is.na(t_next_event)) t_next_event = inf;

        if (DEBUG) Rcout << "next NI event id = " << id_next_event << " at t = " << signif(t_next_event, 5)) endl;

        // generate random timestep ----
        double total_event_rate = sum(pop["event_rate"]);

        // calculate dt if infections event rate > 0
        double dt = (total_event_rate > 0.0) ? R::rexp(rate = total_event_rate) : inf;

        if (DEBUG) Rcout << "Total infections event rate = " << signif(total_event_rate, 5)) << endl;

        // check if next event is infection or non-infection ----
        if (epi_time + dt < t_next_event) {
            if (DEBUG) Rcout << "next event is infection at t = " << signif(epi_time, 5)) << endl;

            epi_time = epi_time + dt;

            // randomly select individual
            id_next_event = sample(nrow(pop),
                                    size = 1L,
                                    prob = pop$event_rate)

            group_id = pop$group[id_next_event]
            infectives = pop[, .(.I, group, status, inf)][group == group_id & status %in% c("I", "D")]
            infd_by = safe_sample(x = infectives$I,
                                   size = 1L,
                                   prob = infectives$inf)
            next_gen = pop$generation[infd_by] + 1L

            set(pop, id_next_event, c("status", "Tinf", "Tinc", "Tsym", "Trec"),
                generate_seidr_path(epi_time, pop, id_next_event, params))
            set(pop, id_next_event, c("generation", "infected_by"), list(next_gen, infd_by))

            if (DEBUG) Rcout << "ID " << id_next_event << ": S -> E, infected by ID " << infd_by << endl;
        } else {
            if (DEBUG) Rcout << "next event is non-infection at t = " << signif(t_next_event, 5)) << endl;

            epi_time = t_next_event;

            status = pop$status[id_next_event];

            if (status == "E") {
                if (DEBUG) Rcout << "ID " << id_next_event << ": E -> I" << endl;
                set(pop, id_next_event, "status", "I");
            } else if (status == "I") {
                if (DEBUG) Rcout << "ID " << id_next_event << ": I -> D" << endl;
                set(pop, id_next_event, "status", "D");
            } else if (status == "D") {
                if (DEBUG) Rcout << "ID " << id_next_event << ": D -> R" << endl;
                set(pop, id_next_event, "status", "R");
            } else {
                Rcout << "status = " << status << endl;
                print(pop[, .(group, donor, status, Tinf, Tinc, Tsym, Trec, group_inf, event_rate)]);
                print(pop[id_next_event, .(group, donor, status, Tinf, Tinc, Tsym, Trec, group_inf, event_rate)]);
                stop("selected ID ", id_next_event, "... unexpected event!");
                break;
            }
        }

        if (DEBUG) {
            print(pop[, .(group, donor, status, Tinf, Tinc, Tsym, Trec, group_inf, event_rate)]);
        }
    }
    Rcout << " - Final t = " << signif(epi_time, 5) << ", values are:" << endl;
    print(table(pop$status));

    // tidy up pop
    pop[, c("group_inf", "event_rate") := NULL];

    // fix generation
    pop[status == "R" & donor == 0L, generation := 2L];

    traits2 = rbind(traits[sdp != "progeny"], pop, fill = TRUE);

    return traits2;
}
