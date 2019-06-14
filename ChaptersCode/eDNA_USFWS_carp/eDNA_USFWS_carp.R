## @knitr stanSettings
library(rstan)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## @knitr occ_model

comb = function(n, r){ factorial(n)/(factorial(r) * factorial(n -r ))}

sample_detect_one <-
    Vectorize(function(J = 50,
                       K = 8,
                       theta = 0.06,
                       p_detection  = 0.3
                       ){
        j_index = J:0
	prob = 	sum(
            comb( J, j_index) *
            (1 - theta) ^ j_index *
            (theta * (1 - p_detection)^K )  ^ rev(j_index))
	return(1 - prob)
    })

## @knitr occ_model_run
results <-
    expand.grid(J = 1:120,
                theta = c(0.05, 0.1, 0.2, 0.4, 0.8, 1.0), 
                p_detection = c(0.05, 0.1, 0.2, 0.4, 0.8, 1.0),
                K = c(2, 4, 8, 16)) %>%
    as_tibble()

results <- 
    results %>%
    mutate(prob_detect = sample_detect_one(J = J, K = K,
                                           theta = theta,
                                           p_detection = p_detection),
           theta_plot = factor(paste0("theta = ", theta)),
           p_detection_plot = factor(paste0("p = ", p_detection)),
           K_plot = factor(paste0("K = ", K))) %>%
    mutate(K_plot = fct_reorder(K_plot, K))


## @knitr occ_code
detect_one <-
    ggplot(data = results, aes(x = J, y = prob_detect, color = K_plot)) +
    geom_line() +
    facet_grid( p_detection_plot ~ theta,
               labeller = label_bquote(cols = theta == .(theta))) +
    theme_minimal() +
    ylab("Probabiltiy of detecting species at site") +
    xlab("Number of samples per site (J)") +
    scale_color_manual("Molecular\nreplicates",
                       values = c("red", "blue", "black", "seagreen",
                                  "orange", "grey50")) +
    scale_x_continuous(breaks = seq(0,125, by = 30)) +
    scale_y_continuous(breaks = seq(0,1, by = .5)) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0),
          panel.spacing = unit(1, "lines"))
print(detect_one)

## @knitr occ_prob

samples_detect <-
    function(
             n_sims = 2,
             J_in = c(10, 100),
             theta_in = c(0.1, 0.9),
             K_in = 8,
             p_detection_in = c(0.25, 0.75)
             ){
        
        results <-
            expand.grid(simulation = 1:n_sims,
                        J = J_in, theta = theta_in,
                        K = K_in, p_detection = p_detection_in) %>%
            as_tibble()

        results <-
            results %>%
            mutate(samples_DNA = rbinom(n = nrow(results),
                                        size = J,
                                        prob = theta)) %>%
            mutate(samples_detect = rbinom(n = nrow(results),
                                         size = samples_DNA,
                                         prob = 1 - (1 - p_detection)^K) ) %>%
            mutate(samples_p_positive = samples_detect /  J,
                   Samples_plot = factor(paste0("J = ", J)),
                   theta_plot = factor(paste0("theta = ", theta)),
                   p_detection_plot = factor(paste0("p = ", p_detection)),
                   K_plot = paste0("K = ", K)) %>%
            mutate(Samples_plot = fct_reorder(Samples_plot, J),
                   theta_plot = fct_reorder(theta_plot, theta),
                   p_detection_plot = fct_reorder(p_detection_plot, p_detection),
                   K_plot = fct_reorder(K_plot, K))

	return(results)
}

## @knitr occ_prob_sim

sample_results <- samples_detect(n_sims = 500,
                                 theta = c(0.05, 0.1, 0.2, 0.4, 0.8, 1.0),
                                 p_detection = c(0.05, 0.1, 0.2, 0.4, 0.8, 1.0),
                                 J = c(5, 10, 20, 40, 80, 120),
                                 K = c(2, 4, 8, 16))

compare_sites <-
    ggplot(sample_results,
           aes(x = K_plot, y = samples_p_positive, color = factor(theta))) +
	geom_boxplot(outlier.size = 0.5) +
	facet_grid(  Samples_plot ~  p_detection_plot ) +
    theme_minimal() +
    ylab(expression("Simulated "*theta~bgroup("(",
                                              over("Number of simulated positive samples",
                                                   "Total number of simulated samples"), ")"))) +
    xlab("Number of molecular replicates (K)") +
    scale_color_manual(expression("Generating "* theta),
                       values = c("red", "blue", "black", "seagreen", "orange", "grey50"))
print(compare_sites)
