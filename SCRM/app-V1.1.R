#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(parallel)
library(lhs)
library(mc2d)
library(bslib)
library(bsicons)
library(shiny)
library(shinyalert)
library(shinyjs)
library(shinyWidgets)
library(DT)
library(beepr)

step_base <- function(model_step, microbial_process, distribution="Normal (log10)", basis, freq, param_1=0, param_2=0, param_3=0, carried_contamination=0, lhs_row, field_mass, unit = NA) {
    # Did the event occur?
    if (is.na(freq)) {
        occurance <- 1
    } else {
        # occurance <- mc2d::rbern(1, freq)
        occurance <- mc2d::qbern(lhs_row[1], freq)
        if (length(lhs_row)!=1) {
            lhs_row <- lhs_row[-1] 
        }
    }
    # Modeled Value
    if (microbial_process == "Risk Output Test" | microbial_process == "Product Test" | microbial_process == "Select a microbial process" | is.na(microbial_process)) {
        modeled_value <- NA
    } else {
        # "Contamination/Removal" | "Increase/Reduction"
        if (grepl("Normal", distribution)) {
            if (occurance == 1) {
                modeled_value <- qnorm(lhs_row[1], param_1, param_2)
            } else {
                modeled_value <- NA
            }
            if (length(lhs_row)!=1) {
                lhs_row <- lhs_row[-1]
            }
        } else if (grepl("Uniform", distribution)) {
            if (occurance == 1) {
                modeled_value <- qunif(lhs_row[1], param_1, param_2)
            } else {
                modeled_value <- NA
            }
            if (length(lhs_row)!=1) {
                lhs_row <- lhs_row[-1]
            }
        } else if (grepl("Pert", distribution)) {
            if (occurance == 1) {
                modeled_value <- mc2d::qpert(lhs_row[1], param_1, param_2, param_3)
            } else {
                modeled_value <- NA
            }
            if (length(lhs_row)!=1) {
                lhs_row <- lhs_row[-1]
            }
        } else if (grepl("Triangular", distribution)) {
            if (occurance == 1) {
                modeled_value <- mc2d::qtriang(lhs_row[1], param_1, param_2, param_3)
            } else {
                modeled_value <- NA
            }
            if (length(lhs_row)!=1) {
                lhs_row <- lhs_row[-1]
            }
        } else {
            modeled_value <- NA
        }
    }
    # Probability of a Positive Grab
    if (microbial_process == "Risk Output Test" | microbial_process == "Product Test") {
        p_positive_grab <- 1-exp(-carried_contamination*(param_1/param_2))
    } else {
        p_positive_grab <- NA
    }
    # Probability of Positive Test
    if (microbial_process == "Risk Output Test" | microbial_process == "Product Test") {
        p_positive_test <- 1-stats::pbinom(0, param_2, p_positive_grab, lower.tail = TRUE)
    } else {
        p_positive_test <- NA
    }
    # FOR OUTPUT Probability of at least one of the test being positive Regulatory
    if (microbial_process == "Risk Output Test") {
        p_at_least_one_test_being_positive_regulatory <- 1-stats::pbinom(0, param_3, p_positive_test, lower.tail = TRUE)
    } else {
        p_at_least_one_test_being_positive_regulatory <- 0
    }
    # Probability of at least one of the test being positive Product Testing
    if (microbial_process == "Product Test") {
        p_at_least_one_test_being_positive_product <- 1-stats::pbinom(0, param_3, p_positive_test, lower.tail = TRUE)
    } else {
        p_at_least_one_test_being_positive_product <- 0
    }
    # Testing Outcome Regulatory (1 = Reject, 0 = Accept)
    if (microbial_process == "Risk Output Test") {
        # testing_outcome_regulatory <- mc2d::rbern(1, p_at_least_one_test_being_positive_regulatory)
        testing_outcome_regulatory <- mc2d::qbern(lhs_row[1], p_at_least_one_test_being_positive_regulatory)
        if (length(lhs_row)!=1) {
            lhs_row <- lhs_row[-1] 
        }
    } else {
        testing_outcome_regulatory <- NA
    }
    # Testing Outcome Testing (1 = Reject, 0 = Accept)
    if (microbial_process == "Product Test") {
        # testing_outcome_product <- mc2d::rbern(1, p_at_least_one_test_being_positive_product)
        testing_outcome_product <- mc2d::qbern(lhs_row[1], p_at_least_one_test_being_positive_product)
        if (length(lhs_row)!=1) {
            lhs_row <- lhs_row[-1] 
        }
    } else {
        testing_outcome_product <- NA
    }
    # Carried Cont
    if (!is.na(testing_outcome_product) & testing_outcome_product == 1) {
        carried_contamination <- 0
    } else {
        if (is.na(modeled_value)) {
            carried_contamination <- carried_contamination
        } else {
            if (microbial_process == "Contamination/Removal") {
                if (occurance == 0) {
                    carried_contamination <- carried_contamination
                } else {
                    if (grepl("log10", distribution)) {
                        modeled_value <- 10^modeled_value
                    } else {
                        modeled_value <- modeled_value
                    }
                    if (basis == "absolute") {
                        if (unit == "g") {
                            carried_contamination <- ((modeled_value)/field_mass)+carried_contamination
                        } else if (unit == "lb") {
                            carried_contamination <- ((modeled_value)/(field_mass*454))+carried_contamination
                        }
                    } else if (basis == "per g") {
                        carried_contamination <- (modeled_value)+carried_contamination
                    } else if (basis == "per lb") {
                        carried_contamination <- ((modeled_value)/454)+carried_contamination
                    }
                    if (carried_contamination < 0) {
                        carried_contamination <- 0
                    }
                }
            } else {
                # If Increase/Reduction
                if (carried_contamination == 0) {
                    carried_contamination <- 0
                } else {
                    carried_contamination <- 10^(log10(carried_contamination)+modeled_value)
                }
            }
        }
    }
    # Total tests product testing for cost analysis
    if (model_step == "Product Test") {
        total_tests <- param_3
    } else {
        total_tests <- 0
    }
    output <- list(
        model_step = model_step,
        microbial_process = microbial_process,
        distribution = distribution,
        basis = basis,
        field_mass = field_mass,
        unit = unit,
        freq = freq,
        occurance = occurance,
        param_1 = param_1,
        param_2 = param_2,
        param_3 = param_3,
        modeled_value = modeled_value,
        p_positive_grab = p_positive_test,
        p_at_least_one_test_being_positive_regulatory = p_at_least_one_test_being_positive_regulatory,
        p_at_least_one_test_being_positive_product = p_at_least_one_test_being_positive_product,
        testing_outcome_regulatory = testing_outcome_regulatory,
        testing_outcome_product = testing_outcome_product,
        carried_contamination = carried_contamination,
        total_tests = total_tests
    )
    return(output)
}

scrm_core <- function(i, total_iterations, steps, lhs_matrix, field_mass, serving_mass_g, alpha, beta, unit=NA) {
    lhs_row <- lhs_matrix[i,]
    
    current_step <- 1
    while (current_step <= dim(steps)[1]) {
        x <- 1
        
        if (!is.na(steps$freq[current_step])) {
            x <- x + 1
        }
        if (steps$microbial_process[current_step] == "Contamination/Removal" | steps$microbial_process[current_step] == "Increase/Reduction") {
            x <- x + 1
        }  else if (steps$microbial_process[current_step] == "Risk Output Test") {
            x <- x + 1
        }  else if (steps$microbial_process[current_step] == "Product Test") {
            x <- x + 1
        }
        
        if (current_step==1) {
            results <- data.frame(step_base(steps[current_step,1], steps[current_step,2], steps[current_step,3], steps[current_step,4], steps[current_step,5], steps[current_step,6], steps[current_step,7], steps[current_step,8], steps[current_step,9], lhs_row, field_mass, unit))
        } else {
            results <- rbind(results, data.frame(step_base(steps[current_step,1], steps[current_step,2], steps[current_step,3], steps[current_step,4], steps[current_step,5], steps[current_step,6], steps[current_step,7], steps[current_step,8], as.numeric(tail(as.data.frame(results)$carried_contamination, n = 1)), lhs_row = lhs_row, field_mass, unit=unit)))
        }
        if (current_step <= dim(steps)[1]) {
            lhs_row <- lhs_row[x:length(lhs_row)] 
        }
        current_step <- current_step + 1
    }
    
    results <- as.data.frame(results)
    rownames(results) <- seq(1, dim(results)[1])
    
    # field_mass <- 160000
    # field_mass_g <- field_mass*454
    # serving_mass_g <- 85
    # total_servings <- field_mass_g/serving_mass_g
    # alpha <- 0.267
    # beta <- 229.2928
    dose_per_serving <- as.numeric(tail(results$carried_contamination, 1))*serving_mass_g
    p_illness_per_serving <- (1-(1+dose_per_serving/beta)^-alpha)*1
    
    # Regulatory Testing
    p_positive_regulatory_test <- 1-(prod(1-as.numeric(results$p_at_least_one_test_being_positive_regulatory)))
    
    # Regulatory Testing on Contaminated Runs only
    if (results$occurance[which(results$microbial_process == "Risk Output Test")] == 1) {
        p_positive_regulatory_test_given_contamination <- p_positive_regulatory_test
    } else {
        p_positive_regulatory_test_given_contamination <- NA
    }
    
    # Product Testing Outcomes
    if (sum(as.numeric(results$testing_outcome_product[which(!is.na(results$testing_outcome_product))]))>0) {
        overall_outcome <- 1
    } else {
        overall_outcome <- 0
    }
    
    if (all(is.na(results$testing_outcome_regulatory))) {
        testing_outcome_regulatory_sum <- NA
    } else {
        testing_outcome_regulatory_sum <- sum(as.numeric(results$testing_outcome_regulatory[which(!is.na(results$testing_outcome_regulatory))]))
    }
    
    if (all(is.na(results$testing_outcome_product))) {
        testing_outcome_product_sum <- NA
    } else {
        testing_outcome_product_sum <- sum(as.numeric(results$testing_outcome_product[which(!is.na(results$testing_outcome_product))]))
    }
    
    single_run <- list(
        iteration = i,
        # steps = steps,
        # results = results,
        dose_per_serving = dose_per_serving,
        p_illness_per_serving = p_illness_per_serving,
        p_positive_regulatory_test = p_positive_regulatory_test,
        p_positive_regulatory_test_given_contamination = p_positive_regulatory_test_given_contamination,
        overall_outcome = overall_outcome,
        testing_outcome_regulatory_sum = sum(as.numeric(results$testing_outcome_regulatory[which(!is.na(results$testing_outcome_regulatory))])),
        testing_outcome_product_sum = sum(as.numeric(results$testing_outcome_product[which(!is.na(results$testing_outcome_product))]))
    )
    return(single_run)
}

scrm_clustered <- function(total_iterations, total_simulations, cores=detectCores(), steps, lhs_matrix, field_mass, serving_mass_g, alpha, beta, cl, parallel = TRUE, unit) {
    if (parallel == TRUE) {
        clusterEvalQ(cl, {
            library(mc2d)
        })
        model_runs <- parLapply(cl, 1:total_iterations, scrm_core, total_iterations, steps, lhs_matrix, field_mass, serving_mass_g, alpha, beta, unit)
        stopCluster(cl)
        # write.csv(model_runs, "model_runs.csv")
        model_runs_tidy <- data.table::rbindlist(model_runs)
        rm(model_runs)
        gc()
    } else {
        j <- 1
        while (j <= total_iterations) {
            if (j == 1) {
                model_runs_tidy <- as.data.frame(scrm_core(j, total_iterations, steps, lhs_matrix, field_mass, serving_mass_g, alpha, beta, unit))
            } else {
                model_runs_tidy <- rbind(model_runs_tidy, as.data.frame(scrm_core(j, total_iterations, steps, lhs_matrix, field_mass, serving_mass_g, alpha, beta, unit)))
            }
            j <- j+1
        }
    }
    
    # model_runs_tidy <- model_runs_tidy %>%
    #   distinct()
    
    model_runs_tidy <- as.data.frame(model_runs_tidy)
    rownames(model_runs_tidy) <- seq(1, dim(model_runs_tidy)[1])
    # write.csv(model_runs_tidy, "model_runs.csv", row.names = FALSE)
    return(model_runs_tidy)
}

risk_binner <- function(model_runs, iterations) {
    bins <- as.data.frame(as.numeric(model_runs$p_positive_regulatory_test))
    colnames(bins) <- c("p_positive")
    bins$fraction_p <- 1/bins$p_positive
    bins$fraction_p[which(is.infinite(bins$fraction_p))] <- 0
    
    positive_risk_bins <- data.frame(
        counts = c(
            baseline = sum(bins$fraction_p>1),
            greater_than_1_in_10 = sum(bins$fraction_p>10),
            greater_than_1_in_100 = sum(bins$fraction_p>100),
            greater_than_1_in_1000 = sum(bins$fraction_p>1000),
            greater_than_1_in_10000 = sum(bins$fraction_p>10000),
            greater_than_1_in_100000 = sum(bins$fraction_p>100000),
            greater_than_1_in_1000000 = sum(bins$fraction_p>1000000),
            smaller_than_1_in_10000 = NA
        ),
        difference = c(
            rep(NA, 8)
        )
    )
    
    positive_risk_bins$difference[1] <- NA
    positive_risk_bins$difference[2] <- positive_risk_bins$counts[1]-positive_risk_bins$counts[2]
    positive_risk_bins$difference[3] <- positive_risk_bins$counts[2]-positive_risk_bins$counts[3]
    positive_risk_bins$difference[4] <- positive_risk_bins$counts[3]-positive_risk_bins$counts[4]
    positive_risk_bins$difference[5] <- positive_risk_bins$counts[4]-positive_risk_bins$counts[5]
    positive_risk_bins$difference[6] <- positive_risk_bins$counts[5]-positive_risk_bins$counts[6]
    positive_risk_bins$difference[7] <- positive_risk_bins$counts[6]-positive_risk_bins$counts[7]
    positive_risk_bins$difference[8] <- iterations-sum(positive_risk_bins$difference[2:6])
    
    positive_risk_bins$difference_pct <- positive_risk_bins$difference/iterations*100
    
    # detected by retail?|means retail
    positive_risk_bins$detected_by_retail_counts <- NA
    positive_risk_bins$detected_by_retail_counts[1] <- sum(as.logical(bins$fraction_p>1)*as.logical(model_runs$testing_outcome_regulatory_sum))
    positive_risk_bins$detected_by_retail_counts[2] <- sum(as.logical(bins$fraction_p>10)*as.logical(model_runs$testing_outcome_regulatory_sum))
    positive_risk_bins$detected_by_retail_counts[3] <- sum(as.logical(bins$fraction_p>100)*as.logical(model_runs$testing_outcome_regulatory_sum))
    positive_risk_bins$detected_by_retail_counts[4] <- sum(as.logical(bins$fraction_p>1000)*as.logical(model_runs$testing_outcome_regulatory_sum))
    positive_risk_bins$detected_by_retail_counts[5] <- sum(as.logical(bins$fraction_p>10000)*as.logical(model_runs$testing_outcome_regulatory_sum))
    positive_risk_bins$detected_by_retail_counts[6] <- sum(as.logical(bins$fraction_p>100000)*as.logical(model_runs$testing_outcome_regulatory_sum))
    
    # detected by fp testing?|means product testing
    positive_risk_bins$detected_by_fp_testing_counts <- NA
    positive_risk_bins$detected_by_fp_testing_counts[1] <- sum(as.logical(bins$fraction_p>1)*as.logical(model_runs$testing_outcome_product_sum))
    positive_risk_bins$detected_by_fp_testing_counts[2] <- sum(as.logical(bins$fraction_p>10)*as.logical(model_runs$testing_outcome_product_sum))
    positive_risk_bins$detected_by_fp_testing_counts[3] <- sum(as.logical(bins$fraction_p>100)*as.logical(model_runs$testing_outcome_product_sum))
    positive_risk_bins$detected_by_fp_testing_counts[4] <- sum(as.logical(bins$fraction_p>1000)*as.logical(model_runs$testing_outcome_product_sum))
    positive_risk_bins$detected_by_fp_testing_counts[5] <- sum(as.logical(bins$fraction_p>10000)*as.logical(model_runs$testing_outcome_product_sum))
    positive_risk_bins$detected_by_fp_testing_counts[6] <- sum(as.logical(bins$fraction_p>100000)*as.logical(model_runs$testing_outcome_product_sum))
    
    positive_risk_bins$difference_retail <- 0
    positive_risk_bins$difference_retail[2] <- positive_risk_bins$detected_by_retail_counts[1]-positive_risk_bins$detected_by_retail_counts[2] # Greater than 1 in 10 lots
    positive_risk_bins$difference_retail[3] <- positive_risk_bins$detected_by_retail_counts[2]-positive_risk_bins$detected_by_retail_counts[3] # Greater than 1 in 100 lots
    positive_risk_bins$difference_retail[4] <- positive_risk_bins$detected_by_retail_counts[3]-positive_risk_bins$detected_by_retail_counts[4] # Greater than 1 in 1,000 lots
    positive_risk_bins$difference_retail[5] <- positive_risk_bins$detected_by_retail_counts[4]-positive_risk_bins$detected_by_retail_counts[5] # Greater than 1 in 10,000 lots
    positive_risk_bins$difference_retail[6] <- positive_risk_bins$detected_by_retail_counts[5]-positive_risk_bins$detected_by_retail_counts[6] # Greater than 1 in 100,000 lots
    positive_risk_bins$difference_retail[8] <- NA
    
    positive_risk_bins$difference_product_testing <- 0
    positive_risk_bins$difference_product_testing[2] <- positive_risk_bins$detected_by_fp_testing_counts[1]-positive_risk_bins$detected_by_fp_testing_counts[2] # Greater than 1 in 10 lots
    positive_risk_bins$difference_product_testing[3] <- positive_risk_bins$detected_by_fp_testing_counts[2]-positive_risk_bins$detected_by_fp_testing_counts[3] # Greater than 1 in 100 lots
    positive_risk_bins$difference_product_testing[4] <- positive_risk_bins$detected_by_fp_testing_counts[3]-positive_risk_bins$detected_by_fp_testing_counts[4] # Greater than 1 in 1,000 lots
    positive_risk_bins$difference_product_testing[5] <- positive_risk_bins$detected_by_fp_testing_counts[4]-positive_risk_bins$detected_by_fp_testing_counts[5] # Greater than 1 in 10,000 lots
    positive_risk_bins$difference_product_testing[6] <- positive_risk_bins$detected_by_fp_testing_counts[5]-positive_risk_bins$detected_by_fp_testing_counts[6] # Greater than 1 in 100,000 lots
    positive_risk_bins$difference_product_testing[8] <- NA
    
    return(positive_risk_bins)
}

risk_binner_02 <- function(model_runs) {
    bins <- as.data.frame(as.numeric(model_runs$p_positive_regulatory_test))
    colnames(bins) <- c("p_positive")
    bins$fraction_p <- 1/bins$p_positive
    bins$fraction_p[which(is.infinite(bins$fraction_p))] <- 0
    
    positive_risk_bins <- data.frame(
        counts = c(
            x_gt_1 = length(which(bins$fraction_p > 1)),
            x_bt_1_10 = length(which(10 > bins$fraction_p & bins$fraction_p >= 1)),
            x_bt_10_100 = length(which(100 > bins$fraction_p & bins$fraction_p >= 10)),
            x_bt_100_1000 = length(which(1000 > bins$fraction_p & bins$fraction_p >= 100)),
            x_bt_1000_10000 = length(which(10000 > bins$fraction_p & bins$fraction_p >= 1000)),
            x_gt_10000 = length(which(bins$fraction_p >= 10000)),
            x_et_0 = length(which(bins$fraction_p == 0))
        ),
        regulatory_tests = c(
            x_gt_1 = length(which(bins$fraction_p > 1 & model_runs$testing_outcome_regulatory_sum > 0)),
            x_bt_1_10 = length(which(10 > bins$fraction_p & bins$fraction_p >= 1 & model_runs$testing_outcome_regulatory_sum > 0)),
            x_bt_10_100 = length(which(100 > bins$fraction_p & bins$fraction_p >= 10 & model_runs$testing_outcome_regulatory_sum > 0)),
            x_bt_100_1000 = length(which(1000 > bins$fraction_p & bins$fraction_p >= 100 & model_runs$testing_outcome_regulatory_sum > 0)),
            x_bt_1000_10000 = length(which(10000 > bins$fraction_p & bins$fraction_p >= 1000 & model_runs$testing_outcome_regulatory_sum > 0)),
            x_gt_10000 = length(which(bins$fraction_p >= 10000 & model_runs$testing_outcome_regulatory_sum > 0)),
            x_et_0 = length(which(bins$fraction_p == 0 & model_runs$testing_outcome_regulatory_sum == 0))
        )
    )
}

lhs_matrix_width <- function(steps) {
    i <- 1
    lhs_width <- 0
    while (i <= dim(steps)[1]) {
        if (!is.na(steps$freq[i]))  {
            lhs_width <- lhs_width + 1
        }
        if (steps$microbial_process[i] == "Contamination/Removal" | steps$microbial_process[i] == "Increase/Reduction") {
            lhs_width <- lhs_width + 1
        }  else if (steps$microbial_process[i] == "Risk Output Test") {
            lhs_width <- lhs_width + 1
        }  else if (steps$microbial_process[i] == "Product Test") {
            lhs_width <- lhs_width + 1
        }
        i <- i+1
    }
    return(lhs_width)
}

step_cleaner <- function(steps) {
    steps <- data.frame(
        "model_step" = steps[,1], 
        "microbial_process" = steps[,2],
        "distribution" = steps[,3],
        "freq" = as.numeric(steps[,4]),
        "param_1" = as.numeric(steps[,5]),
        "param_2" = as.numeric(steps[,6]),
        "param_3" = as.numeric(steps[,7]),
        "carried_contamination" = as.numeric(steps[,8]),
        stringsAsFactors = FALSE
    )
    rownames(steps) <- seq(1, dim(steps)[1])
    return(steps)
}

p_positive_regulatory_summarizer <- function(model_runs) {
    summary <- data.frame(
        mean = mean(as.numeric(model_runs$p_positive_regulatory_test), na.rm = TRUE),
        mean_geo = exp(mean(log(model_runs$p_positive_regulatory_test[model_runs$p_positive_regulatory_test>0]))),
        mean_arithmetic_of_nonzeroes = mean(model_runs$p_positive_regulatory_test[model_runs$p_positive_regulatory_test>0]),
        sd = sd(as.numeric(model_runs$p_positive_regulatory_test), na.rm = TRUE),
        percentile_2.5th = as.numeric(quantile(as.numeric(model_runs$p_positive_regulatory_test), 0.025)),
        percentile_97.5th = as.numeric(quantile(as.numeric(model_runs$p_positive_regulatory_test), 0.975)),
        mean_given_contamination = mean(as.numeric(model_runs$p_positive_regulatory_test_given_contamination), na.rm = TRUE),
        sd_given_contamination = sd(as.numeric(model_runs$p_positive_regulatory_test_given_contamination), na.rm = TRUE),
        percentile_2.5th_given_contamination = as.numeric(quantile(as.numeric(model_runs$p_positive_regulatory_test_given_contamination), 0.025)),
        percentile_97.5th_given_contamination = as.numeric(quantile(as.numeric(model_runs$p_positive_regulatory_test_given_contamination), 0.975))
    )
    return(summary)
}

# Lists ----
list_contamination_distributions <- c("Normal (log10)", "Uniform (log10)", "Pert (log10)", "Triangular (log10)", "Normal (linear)", "Uniform (linear)", "Pert (linear)", "Triangular (linear)")
list_distributions <- c("Normal (log10)", "Uniform (log10)", "Pert (log10)", "Triangular (log10)")
list_processes <- c("Contamination/Removal", "Increase/Reduction", "Risk Output Test", "Product Test")
list_contamination_types <- c("absolute", "per g", "per lb")

list_scenarios <- list.files(paste0(getwd(), "/scenarios"), full.names = TRUE)
list_scenarios_data <- lapply(list_scenarios, read.csv)
list_scenario_names <- gsub(".csv", "", list.files(paste0(getwd(), "/scenarios"), full.names = FALSE))
if (length(list_scenario_names) == 0) {
    list_scenario_names <- "NA"
    list_scenario_position <- NULL
} else {
    list_scenario_position <- "Baseline"
}

# Glossary ----
glossary <- list(
    accordion_panel("Total Mass", "The size of the lot (user-defined) you would like the iteration to represent."),
    accordion_panel("Mass Unit", "lb or g. This is the unit tied to the field mass."),
    accordion_panel("Iterations", "This is the number of times the tool will run the scenario you have defined."),
    accordion_panel("Step", "A modular component of the supply chain model. Microbial contamination is sequentially transferred between steps by converting it to the number of cells per gram (CFU/g)."),
    accordion_panel("Step Name", "A free-text description of the step."),
    accordion_panel("Step Type", HTML("Selection of the mechanic for the step at its respective point in the supply chain. Options include:<ul><li><b>Contamination/Removal</b></li><li><b>Increase/Reduction</b></li><li><b>Risk Output Test</b></li><li><b>Product Test</b></li></ul>")),
    accordion_panel("Contamination/Removal ", "A step type where microbial contamination is added to or removed from the system. This can be entered as an absolute number of cells (CFUs) or as cells per unit mass (CFU/g; CFU/lb) in linear or log scale. Cells are applied uniformly across the defined field mass (i.e., \"lot\")."),
    accordion_panel("Increase/Reduction", "A step type where any microbial contamination existing in the system at that point can be increased (representing growth) or reduced (representing die-off or reductions from processing steps). Increases and reductions must be entered per gram on a log scale."),
    accordion_panel("Product Test", "A step type where a mass in grams, number of grabs, and number of tests are defined. Each test will yield a probability of being positive based on the contamination level in the lot at that point. The result of the test is also determined based on this probability. If the test is positive, the lot is \"rejected\" and now has 0 probability of testing positive in subsequent steps. This can occur anywhere but must precede the risk output test."),
    accordion_panel("Risk Output Test", "This is mathematically identical to a product test, but also the defined model output. The probability of the risk output test being positive defines the risk category the lot is placed in."),
    accordion_panel("Probability of Occurrence (%)", "This is input as a percentage (0% to 100%) and represents the likelihood of a step occurring during an iteration. If the step does not occur, it is skipped, and the model will move to the next step."),
    accordion_panel("Distribution", HTML("The user selects the distribution of their input values, options include:<ul><li><b>Normal</b></li><li><b>Uniform</b></li><li><b>Pert</b></li><li><b>Triangular</b></ul>")),
    accordion_panel("Absolute or per mass", "The user defines if the input values of a distribution are per unit mass (e.g., CFU/g, CFU/lb) or if the input value is an absolute number (CFU). If absolute is selected, then the value drawn from that distribution when the model is run will be divided by the Field Mass provided by the user."),
    accordion_panel("Mean [CFU/g or log (CFU/g)] - Normal Distribution only", HTML("For <b>Contamination/Removal</b> step type, this is the mean microbial concentration that is added to or removed from the system.<br>For <b>Increase/Reduction</b> step type, this is the change in the existing microbial concentration that results from applying this step type.")),
    accordion_panel("Standard Deviation (SD) [CFU/g or log (CFU/g)] - Normal Distribution only", HTML("For <b>Contamination/Removal</b> step type, this is the variation about the microbial concentration which is added or removed from the system.<br>For <b>Increase/Reduction</b> step type, this is the variation about the change to the microbial concentration.")),
    accordion_panel("Minimum [CFU/g or log (CFU/g)] - Pert or Triangular Distributions only", HTML("For <b>Contamination/Removal</b> step type, this is the minimum microbial concentration which is added or removed from the system.<br>For <b>Increase/Reduction</b> step type, this is the minimum change in the existing microbial concentration that results from applying this step type.")),
    accordion_panel("Mode [CFU/g or log (CFU/g)] - Pert or Triangular Distributions only", HTML("For <b>Contamination/Removal</b> step type, this is the most likely microbial concentration which is added or removed from the system.<br>For <b>Increase/Reduction</b> step type, this is the most likely change in the existing microbial concentration that results from applying this step type.")),
    accordion_panel("Maximum [CFU/g or log (CFU/g)] - Pert or Triangular Distributions only", HTML("For <b>Contamination/Removal</b> step type, this is the maximum microbial concentration which is added or removed from the system.<br>For <b>Increase/Reduction</b> step type, this is the maximum change in the existing microbial concentration that results from applying this step type.")),
    accordion_panel("Lower Bound [CFU/g or log (CFU/g)] - Uniform Distribution only", HTML("For <b>Contamination/Removal</b> step type, this is the lower bound microbial concentration which is added or removed from the system.<br>For <b>Increase/Reduction</b> step type, this is the lower bound change in the existing microbial concentration that results from applying this step type.")),
    accordion_panel("Upper Bound [CFU/g or log (CFU/g)] - Uniform Distribution only", HTML("For <b>Contamination/Removal</b> step type, this is the upper bound microbial concentration which is added or removed from the system.<br>For <b>Increase/Reduction</b> step type, this is the upper bound change in the existing microbial concentration that results from applying this step type.")),
    accordion_panel("Composite Mass (g)", HTML("This is the mass of the sample taken for a <b>Product Test</b> or <b>Risk Output Test</b>. The sample is assumed to be composited.")),
    accordion_panel("Tests (N)", HTML("This is the number of tests that occur for a <b>Product Test</b> or <b>Risk Output Test</b>.")),
    accordion_panel("Grabs (N)", HTML("This is the number of grabs that occur for a <b>Product Test</b> or <b>Risk Output Test</b>. Currently, because contamination is homogenously applied across the field mass, number of grabs will not change the result of the test.")),
    accordion_panel("Risk Category", HTML("This is the category an iterated lot is placed in based on its probability of the risk output test being positive.<ul><li>Highest risk is defined as > 1 in 10 risk of a positive test</li><li>High risk is defined as between 1 in 10 and 1 in 100 risk of a positive test</li><li>Medium-High risk is defined as between 1 in 100 and 1 in 1,000 risk of a positive test</li><li>Medium-Low risk is defined as between 1 in 1,000 and 1 in 10,000 risk of a positive test</li><li>Low risk is defined as < 1 in 10,000 risk of a positive test.</li><li>Lowest risk is defined as any lots which were not contaminated or rejected by a product test before the <b>Risk Output Test</b></li></ul>")),
    accordion_panel("Overall Risk of a Positive Test ", HTML("This is the average probability of the <b>Risk Output Test</b> being positive across all iterated lots. Under the current preset baseline scenario, this is our proxy for recall risk as the <b>Risk Output Test</b> occurs after the <b>Presentation to Consumer (Retail)</b> step.")),
    accordion_panel("Count of Lots with Highest Risk", HTML("This is the number of lots in the highest risk bin, which have a greater than a 1-in-10 risk of a positive <b>Risk Output Test</b>. These lots are most likely to cause a public health event and are our proxy for public health risk."))
)

baseline_scenario_terminology <- list (
    accordion_panel("Initial Contamination", HTML("A <b>Contamination/Removal</b> step that represents the introduction of preharvest microbial contamination. (It is assumed that no additional <b>Contamination/Removal</b> or <b>Increase/Reduction</b> steps occur before the product reaches the <b>Processing</b> step)")),
    accordion_panel("Primary Raw Material Production", "A step that does not currently have an input on the model but is present as a stand-in for future modifications to baseline."),
    accordion_panel("Harvesting", "A step that does not currently have an input on the model but is present as a stand-in for future modifications to baseline."),
    accordion_panel("Processing", HTML("An <b>Increase/Reduction</b> step that represents the aggregate of multiple smaller reductions from a prewash and wash step. ")),
    accordion_panel("Presentation to Consumer (Retail)", HTML("An <b>Increase/Reduction</b> step that represents the aggregate of <i>E. coli</i> change from transportation to retail and at retail display (we assume no growth occurs due to proper temperature management).")),
    accordion_panel("Risk Output Testing", HTML("A <b>Risk Output Test</b> step that can represent a sample taken at a retail store by a regulatory agency, or other body. This captures the contamination state of the lot at this point in the supply chain and uses this to calculate the probability of the Risk Output Test being positive.")),
    accordion_panel("Consumer Handling", HTML("An <b>Increase/Reduction</b> step that represents the aggregate of <i>E. coli</i> change from transportation to consumer and at consumer home. This step currently does not affect the modeled output as it occurs after the <b>Risk Output Test</b>."))
)

full_accordion <- list (
    accordion_panel(HTML("<span style=\"font-size:21pt\">Glossary</span>"), accordion(!!!glossary, open = FALSE)),
    accordion_panel(HTML("<span style=\"font-size:21pt\">Baseline Scenario Terminology</span>"), accordion(!!!baseline_scenario_terminology, open = FALSE))
)

# Page 01 ----
page_01 <- fluidPage(
    shinyjs::useShinyjs(),
    sidebarLayout(
        sidebarPanel(width = 2,
            # style = "position:fixed;",
            fluidRow(
                column(
                    selectInput("slctin_scenario", "Preset Scenarios", choices = list_scenario_names, selected = list_scenario_position),
                    actionBttn("btn_load_scenario", "Load Scenario", size = "xs", color = "success"),
                    hr(),
                    textInput("txtin_scenario_name", "Scenario Name", "", width = "100%"),
                    numericInput("nmcin_field_mass", "Field Mass", 160000, width = "100%"),
                    selectInput("slctin_field_mass_unit", "Mass Unit", choices = c("lb", "g"), selected = "lb", width = "100%"),
                    # (numericInput("nmcin_serving_mass", "Serving Mass", 85, width = "50%")),
                    # fluidRow(numericInput("nmcin_alpha", "Alpha", 0.267)),
                    # fluidRow(numericInput("nmcin_beta", "Beta", 229.2928)),
                    numericInput("nmcin_n_iterations", "Iterations", 1E5, width = "100%"),
                    actionBttn("btn_reset", "Reset Inputs", size = "xs", color = "warning"),
                    actionBttn("btn_run", "Run", size = "xs", color = "success"),
                    hr(),
                    selectInput("slctin_save_slot", "Slot", choices = seq(1,25), selected = 1, width = "100%"),
                    actionBttn("btn_save_to_slot", HTML("Save to Slot"), size = "xs", color = "success", width = "40%"),
                    hr(),
                    checkboxInput("ckbxin_beep", "Beep?", value = FALSE, width = NULL),
                    width = 12, align="center"
                ),
            ),
        ),
        mainPanel(width = 10,
            # style = "position: absolute; width:90%;",
            tags$head(tags$style(HTML(".bslib-card, .tab-content, .tab-pane, .card-body {overflow: visible !important;}"))),
            card(
                class = "foo",
                card_header("Results"),
                card_body(class = "foo",
                          fluidRow(
                              h5("Summary Results for Main Risk Outputs"),
                              tableOutput("summary_table_counts"),
                              htmlOutput("summary_text")
                          )
                ),
                hr(),
                card_body(class = "foo",
                          fluidRow(
                              h5("Number of Lots with a Given Risk of a Positive Test Result for the Risk Output Test"),
                              tableOutput("lot_table_counts")
                          ),
                          column(
                              downloadBttn("dwnldbtwn_results_counts", label = "Download Results", color = "success", size = "xs"),
                              width = 4
                          ),
                ),
                hr(),
                card_body(class = "foo",
                          fluidRow(
                              h5("Number of Positive Tests for the Risk Output Test by Risk Category"),
                              tableOutput("lot_table_regulatory_tests")
                          ),
                          column(
                              downloadBttn("dwnldbtwn_results_regulatory_tests", label = "Download Results", color = "success", size = "xs"),
                              width = 4
                          ),
                ),
            ),
            ## Card 01 ----
            card(
                class = "foo",
                card_header(column(fluidRow(checkboxInput("ckbxin_crd01", "Step 01", FALSE), HTML("Step Name"), textInput("txtin_crd01_model_step", NULL, ""), style='padding:0px;'), width = 6)),
                card_body(
                    class = "foo",
                    layout_column_wrap(
                        selectInput(
                            "slctin_crd01_microbial_process",
                            label = tooltip(trigger = list("Step Type", bs_icon("info-circle")), "", id = "tltp_crd01_microbial_process"),
                            choices = list_processes,
                            selected = list_processes[1]
                        ),
                        numericInput(
                            "nmcin_crd01_p_occurrence",
                            tags$span(style="font-size:95%;", "Probability of Occurrence (%)"),
                            100,
                            0,
                            100
                        ),
                        selectInput(
                            "slctin_crd01_distribution",
                            "Distribution",
                            choices = list_distributions,
                            selected = list_distributions[1]
                            ),
                        selectInput(
                            "slctin_crd01_contamination_basis",
                            "Absolute or per mass",
                            choices = list_contamination_types,
                            selected = list_contamination_types[2]
                        )
                    ),
                ),
                card_body(
                    fluidRow(
                        numericInput("nmcin_crd01_parameter_1", "Param 1", 0, width = "25%"),
                        numericInput("nmcin_crd01_parameter_2", "Param 2", 0, width = "25%"),
                        numericInput("nmcin_crd01_parameter_3", "Param 3", 0, width = "25%")
                    )
                )
            ),
            ## Card 02 ----
            card(
                class = "foo",
                card_header(column(fluidRow(checkboxInput("ckbxin_crd02", "Step 02", FALSE), HTML("Step Name"), textInput("txtin_crd02_model_step", NULL, ""), style='padding:0px;'), width = 6)),
                card_body(
                    class = "foo",
                    layout_column_wrap(
                        selectInput(
                            "slctin_crd02_microbial_process",
                            label = tooltip(trigger = list("Step Type", bs_icon("info-circle")), "", id = "tltp_crd02_microbial_process"),
                            choices = list_processes,
                            selected = list_processes[2]
                        ),
                        numericInput(
                            "nmcin_crd02_p_occurrence",
                            tags$span(style="font-size:95%;", "Probability of Occurrence (%)"),
                            100,
                            0,
                            100
                        ),
                        selectInput(
                            "slctin_crd02_distribution",
                            "Distribution",
                            choices = list_distributions,
                            selected = list_distributions[1]
                        ),
                        selectInput(
                            "slctin_crd02_contamination_basis",
                            "Absolute or per mass",
                            choices = list_contamination_types,
                            selected = list_contamination_types[2]
                        )
                    ),
                ),
                card_body(
                    fluidRow(
                        numericInput("nmcin_crd02_parameter_1", "Param 1", 0, width = "25%"),
                        numericInput("nmcin_crd02_parameter_2", "Param 2", 0, width = "25%"),
                        numericInput("nmcin_crd02_parameter_3", "Param 3", 0, width = "25%")
                    )
                )
            ),
            ## Card 03 ----
            card(
                class = "foo",
                card_header(column(fluidRow(checkboxInput("ckbxin_crd03", "Step 03", FALSE), HTML("Step Name"), textInput("txtin_crd03_model_step", NULL, ""), style='padding:0px;'), width = 6)),
                card_body(
                    class = "foo",
                    layout_column_wrap(
                        selectInput(
                            "slctin_crd03_microbial_process",
                            label = tooltip(trigger = list("Step Type", bs_icon("info-circle")), "", id = "tltp_crd03_microbial_process"),
                            choices = list_processes,
                            selected = list_processes[2]
                        ),
                        numericInput(
                            "nmcin_crd03_p_occurrence",
                            tags$span(style="font-size:95%;", "Probability of Occurrence (%)"),
                            100,
                            0,
                            100
                        ),
                        selectInput(
                            "slctin_crd03_distribution",
                            "Distribution",
                            choices = list_distributions,
                            selected = list_distributions[1]
                        ),
                        selectInput(
                            "slctin_crd03_contamination_basis",
                            "Absolute or per mass",
                            choices = list_contamination_types,
                            selected = list_contamination_types[2]
                        )
                    ),
                ),
                card_body(
                    fluidRow(
                        numericInput("nmcin_crd03_parameter_1", "Param 1", 0, width = "25%"),
                        numericInput("nmcin_crd03_parameter_2", "Param 2", 0, width = "25%"),
                        numericInput("nmcin_crd03_parameter_3", "Param 3", 0, width = "25%")
                    )
                )
            ),
            ## Card 04 ----
            card(
                class = "foo",
                card_header(column(fluidRow(checkboxInput("ckbxin_crd04", "Step 04", FALSE), HTML("Step Name"), textInput("txtin_crd04_model_step", NULL, ""), style='padding:0px;'), width = 6)),
                card_body(
                    class = "foo",
                    layout_column_wrap(
                        selectInput(
                            "slctin_crd04_microbial_process",
                            label = tooltip(trigger = list("Step Type", bs_icon("info-circle")), "", id = "tltp_crd04_microbial_process"),
                            choices = list_processes,
                            selected = list_processes[2]
                        ),
                        numericInput(
                            "nmcin_crd04_p_occurrence",
                            tags$span(style="font-size:95%;", "Probability of Occurrence (%)"),
                            100,
                            0,
                            100
                        ),
                        selectInput(
                            "slctin_crd04_distribution",
                            "Distribution",
                            choices = list_distributions,
                            selected = list_distributions[1]
                        ),
                        selectInput(
                            "slctin_crd04_contamination_basis",
                            "Absolute or per mass",
                            choices = list_contamination_types,
                            selected = list_contamination_types[2]
                        )
                    ),
                ),
                card_body(
                    fluidRow(
                        numericInput("nmcin_crd04_parameter_1", "Param 1", 0, width = "25%"),
                        numericInput("nmcin_crd04_parameter_2", "Param 2", 0, width = "25%"),
                        numericInput("nmcin_crd04_parameter_3", "Param 3", 0, width = "25%")
                    )
                )
            ),
            ## Card 05 ----
            card(
                class = "foo",
                card_header(column(fluidRow(checkboxInput("ckbxin_crd05", "Step 05", FALSE), HTML("Step Name"), textInput("txtin_crd05_model_step", NULL, ""), style='padding:0px;'), width = 6)),
                card_body(
                    class = "foo",
                    layout_column_wrap(
                        selectInput(
                            "slctin_crd05_microbial_process",
                            label = tooltip(trigger = list("Step Type", bs_icon("info-circle")), "", id = "tltp_crd05_microbial_process"),
                            choices = list_processes,
                            selected = list_processes[2]
                        ),
                        numericInput(
                            "nmcin_crd05_p_occurrence",
                            tags$span(style="font-size:95%;", "Probability of Occurrence (%)"),
                            100,
                            0,
                            100
                        ),
                        selectInput(
                            "slctin_crd05_distribution",
                            "Distribution",
                            choices = list_distributions,
                            selected = list_distributions[1]
                        ),
                        selectInput(
                            "slctin_crd05_contamination_basis",
                            "Absolute or per mass",
                            choices = list_contamination_types,
                            selected = list_contamination_types[2]
                        )
                    ),
                ),
                card_body(
                    fluidRow(
                        numericInput("nmcin_crd05_parameter_1", "Param 1", 0, width = "25%"),
                        numericInput("nmcin_crd05_parameter_2", "Param 2", 0, width = "25%"),
                        numericInput("nmcin_crd05_parameter_3", "Param 3", 0, width = "25%")
                    )
                )
            ),
            ## Card 06 ----
            card(
                class = "foo",
                card_header(column(fluidRow(checkboxInput("ckbxin_crd06", "Step 06", FALSE), HTML("Step Name"), textInput("txtin_crd06_model_step", NULL, ""), style='padding:0px;'), width = 6)),
                card_body(
                    class = "foo",
                    layout_column_wrap(
                        selectInput(
                            "slctin_crd06_microbial_process",
                            label = tooltip(trigger = list("Step Type", bs_icon("info-circle")), "", id = "tltp_crd06_microbial_process"),
                            choices = list_processes,
                            selected = list_processes[2]
                        ),
                        numericInput(
                            "nmcin_crd06_p_occurrence",
                            tags$span(style="font-size:95%;", "Probability of Occurrence (%)"),
                            100,
                            0,
                            100
                        ),
                        selectInput(
                            "slctin_crd06_distribution",
                            "Distribution",
                            choices = list_distributions,
                            selected = list_distributions[1]
                        ),
                        selectInput(
                            "slctin_crd06_contamination_basis",
                            "Absolute or per mass",
                            choices = list_contamination_types,
                            selected = list_contamination_types[2]
                        )
                    ),
                ),
                card_body(
                    fluidRow(
                        numericInput("nmcin_crd06_parameter_1", "Param 1", 0, width = "25%"),
                        numericInput("nmcin_crd06_parameter_2", "Param 2", 0, width = "25%"),
                        numericInput("nmcin_crd06_parameter_3", "Param 3", 0, width = "25%")
                    )
                )
            ),
            ## Card 07 ----
            card(
                class = "foo",
                card_header(column(fluidRow(checkboxInput("ckbxin_crd07", "Step 07", FALSE), HTML("Step Name"), textInput("txtin_crd07_model_step", NULL, ""), style='padding:0px;'), width = 6)),
                card_body(
                    class = "foo",
                    layout_column_wrap(
                        selectInput(
                            "slctin_crd07_microbial_process",
                            label = tooltip(trigger = list("Step Type", bs_icon("info-circle")), "", id = "tltp_crd07_microbial_process"),
                            choices = list_processes,
                            selected = list_processes[2]
                        ),
                        numericInput(
                            "nmcin_crd07_p_occurrence",
                            tags$span(style="font-size:95%;", "Probability of Occurrence (%)"),
                            100,
                            0,
                            100
                        ),
                        selectInput(
                            "slctin_crd07_distribution",
                            "Distribution",
                            choices = list_distributions,
                            selected = list_distributions[1]
                        ),
                        selectInput(
                            "slctin_crd07_contamination_basis",
                            "Absolute or per mass",
                            choices = list_contamination_types,
                            selected = list_contamination_types[2]
                        )
                    ),
                ),
                card_body(
                    fluidRow(
                        numericInput("nmcin_crd07_parameter_1", "Param 1", 0, width = "25%"),
                        numericInput("nmcin_crd07_parameter_2", "Param 2", 0, width = "25%"),
                        numericInput("nmcin_crd07_parameter_3", "Param 3", 0, width = "25%")
                    )
                )
            ),
            ## Card 08 ----
            card(
                class = "foo",
                card_header(column(fluidRow(checkboxInput("ckbxin_crd08", "Step 08", FALSE), HTML("Step Name"), textInput("txtin_crd08_model_step", NULL, ""), style='padding:0px;'), width = 6)),
                card_body(
                    class = "foo",
                    layout_column_wrap(
                        selectInput(
                            "slctin_crd08_microbial_process",
                            label = tooltip(trigger = list("Step Type", bs_icon("info-circle")), "", id = "tltp_crd08_microbial_process"),
                            choices = list_processes,
                            selected = list_processes[2]
                        ),
                        numericInput(
                            "nmcin_crd08_p_occurrence",
                            tags$span(style="font-size:95%;", "Probability of Occurrence (%)"),
                            100,
                            0,
                            100
                        ),
                        selectInput(
                            "slctin_crd08_distribution",
                            "Distribution",
                            choices = list_distributions,
                            selected = list_distributions[1]
                        ),
                        selectInput(
                            "slctin_crd08_contamination_basis",
                            "Absolute or per mass",
                            choices = list_contamination_types,
                            selected = list_contamination_types[2]
                        )
                    ),
                ),
                card_body(
                    fluidRow(
                        numericInput("nmcin_crd08_parameter_1", "Param 1", 0, width = "25%"),
                        numericInput("nmcin_crd08_parameter_2", "Param 2", 0, width = "25%"),
                        numericInput("nmcin_crd08_parameter_3", "Param 3", 0, width = "25%")
                    )
                )
            ),
            ## Card 09 ----
            card(
                class = "foo",
                card_header(column(fluidRow(checkboxInput("ckbxin_crd09", "Step 09", FALSE), HTML("Step Name"), textInput("txtin_crd09_model_step", NULL, ""), style='padding:0px;'), width = 6)),
                card_body(
                    class = "foo",
                    layout_column_wrap(
                        selectInput(
                            "slctin_crd09_microbial_process",
                            label = tooltip(trigger = list("Step Type", bs_icon("info-circle")), "", id = "tltp_crd09_microbial_process"),
                            choices = list_processes,
                            selected = list_processes[2]
                        ),
                        numericInput(
                            "nmcin_crd09_p_occurrence",
                            tags$span(style="font-size:95%;", "Probability of Occurrence (%)"),
                            100,
                            0,
                            100
                        ),
                        selectInput(
                            "slctin_crd09_distribution",
                            "Distribution",
                            choices = list_distributions,
                            selected = list_distributions[1]
                        ),
                        selectInput(
                            "slctin_crd09_contamination_basis",
                            "Absolute or per mass",
                            choices = list_contamination_types,
                            selected = list_contamination_types[2]
                        )
                    ),
                ),
                card_body(
                    fluidRow(
                        numericInput("nmcin_crd09_parameter_1", "Param 1", 0, width = "25%"),
                        numericInput("nmcin_crd09_parameter_2", "Param 2", 0, width = "25%"),
                        numericInput("nmcin_crd09_parameter_3", "Param 3", 0, width = "25%")
                    )
                )
            ),
        )
    )
)

# Page 02 ----
page_02 <- fluidPage(
    h5("Lot Counts By Risk Category"),
    DT::dataTableOutput("results_lot_table_counts"),
    downloadBttn("dwnldbtwn_results_comparision_counts", label = "Download Table", color = "success", size = "xs"),
    h5("Regulatory Testing Counts"),
    DT::dataTableOutput("results_lot_table_regulatory_tests"),
    downloadBttn("dwnldbtwn_results_comparision_regulatory_tests", label = "Download Table", color = "success", size = "xs")
)

# Page 03 ----
page_03 <- fluidPage(
    useShinyjs(),
    h3("Overview"),
    h5(HTML("This document provides access to a risk management tool developed by The University of Illinois Urbana-Champaign and Cornell University.")),
    h5(HTML("The tool is available upon request (<a href='mailto:mstasie@illinois.edu'>mstasie@illinois.edu</a>, <a href='mailto:cecilwb2@illinois.edu'>cecilwb2@illinois.edu</a>, <a href='mailto:gnpinto2@illinois.edu'>gnpinto2@illinois.edu</a>)")),
    h5(HTML("The tool developed is a flexible supply chain risk model (SCRM) for fresh produce that can model 3 key mechanics, contamination, increase/reduction, and testing, within five main process stages representing farm to consumer supply chain.")),
    h5(HTML("The Supply Chain Risk Model (SCRM) consists of 5 generic produce supply chain process stages: (i) <i>Primary Raw Material Production</i>, (ii) <i>Harvesting</i>, (iii) <i>Processing</i>, (iv) <i>Presentation to the Consumer (Retail)</i>, and (v) <i>Consumer Handling</i>. In any given process stage, any of three mechanics built into this tool can be entered. (i) <i>Contamination/Removal</i>, which is a step where microbial contamination is added to or removed from the system. This can be entered as an absolute number of cells that is applied uniformly across the lot size or as cells per unit mass. Contamination or removal can be entered on log or linear scale. (ii) <i>Increase/Reduction</i> step, which is a step where any microbial contamination in the system at that point can be increased (representing growth) or reduced (representing die-off or reductions from processing steps). Increases and reductions must be entered per unit mass on a log scale. (iii) <i>Testing</i> step, which in this model, always includes the <i>Risk Output Test</i> but may also include one or more <i>Product Tests</i> upstream of the <i>Risk Output Test</i>. Testing steps must specify the number of tests, the mass of the tests, and the number of grabs. Each test will yield a probability of being positive based on the contamination level in the lot at that point. The result of the test is also determined based on this probability.")),
    h5(HTML("We present two model outputs, (i) the <i>Overall Risk</i> of a positive <i>Risk Output Test</i>, which averages the probabilities of positive <i>Risk Output Tests</i> across all iterated lots and is by default set to occur at the <i>Presentation to Consumer (Retail)</i> stage. In this context, we use this metric as a proxy for recall risk. (i) The number of lots categorized as \"highest\" risk, meaning lots with a greater than 1 in 10 chance of the <i>Risk Output Test</i> being positive. In this context, we use this metric as a proxy for public health risk, as lots with this probability of testing positive would be likely to cause a public health event.")),
    h5(HTML("The flexibility of this model allows the user to choose options such as introducing contamination at various stages, reducing contamination, and increasing contamination steps through the system. At any stage, users can define changes in the concentration of the product as it moves through the process. The changes at each stage may be parametrized by multiple distributions such as normal (used in the current baseline), uniform, pert, and triangular. The values for these distributions should be chosen based on peer-reviewed literature review findings or expert opinion.")),
    h5(HTML("Users of this tool should be familiar with generic produce supply chain processes and have a general knowledge of statistics.")),
    br(),
    h3("User Guide Resources"),
    h5("Demo videos coming soon"),
    br(),
    accordion(!!!full_accordion, open = FALSE, multiple = TRUE)
    
)

# Page 04 ----
page_04 <- fluidPage(
    useShinyjs(),
    tags$style(".uiucimg {
                            margin-left:65px;
                            margin-right:0px;
                            margin-top:65px;
                          }"),
    tags$style(".cornellimg {
                            margin-left:130px;
                            margin-right:65px;
                            margin-top:65px;
                          }"),
    titlePanel("Acknowledgements"),
    h3(""),
    h5("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."),
    h5(HTML("")),
    h3("Contacts"),
    h5("Please feel free to contact us if you have any questions."),
    h5(HTML("<b>Matthew Stasiewicz</b>: Principal Investigator, Associate Professor of Applied Food Safety, PhD | Email: "), a("mstasie@illinois.edu",
                                                                                                                               target = "_blank",
                                                                                                                               href = "mailto:mstasie@illinois.edu"
    )),
    h5(HTML("<b>Martin Wiedmann</b>: Co-Principal Investigator, Gellert Family Professor in Food Safety, PhD | Email: "), a("martin.wiedmann@cornell.edu",
                                                                                                                            target = "_blank",
                                                                                                                            href = "mailto:martin.wiedmann@cornell.edu"
    )),
    h5(HTML("<b>Cecil Barnett-Neefs</b>: App Creator, Model Author | Email: "), a("cecilwb2@illinois.edu",
                                                                                  target = "_blank",
                                                                                  href = "mailto:cecilwb2@illinois.edu"
    )),
    h5(HTML("<b>Gabriella Pinto</b>: Model Author | Email: "), a("gnpinto2@illinois.edu",
                                                                 target = "_blank",
                                                                 href = "mailto:gnpinto2@illinois.edu"
    )),
    h5(HTML("<b>Erin Kealey</b>: Scenario Developer | Email: "), a("ekealey@illinois.edu",
                                                             target = "_blank",
                                                             href = "mailto:ekealey@illinois.edu"
    )),
    h5(HTML("<b>Jorge Quintanilla Portillo</b>: Scenario Developer | Email: "), a("jfq@illinois.edu",
                                                                            target = "_blank",
                                                                            href = "mailto:jfq@illinois.edu"
    )),
    h5(HTML("<b>YeonJin Jung</b>: Scenario Developer | Email: "), a("yj354@cornell.edu",
                                                              target = "_blank",
                                                              href = "mailto:yj354@cornell.edu"
    )),
    h3("Acknowledgements"),
    fluidRow(
        # column(1, div(class="uiucimg", (imageOutput("img_uiuc")))),
        # column(1, div(class="cornellimg", (imageOutput("img_cornell"))))
        column(2, div(class="uiucimg", img(src = "University-Wordmark-Full-Color-RGB-TM.png", width = "250px", align = "left"))),
        column(2, div(class="cornellimg", img(src = "bold_cornell_logo_pms187_red.png", width = "250px", align = "left")))
    )
)

# UI ----
ui <- navbarPage(
    theme = bs_theme(bootswatch = "flatly", version = 5),
    title = "Supply Chain Risk Model (in R) - 1.1",
    tabPanel("Model", page_01),
    tabPanel("Results Comparisions", page_02),
    tabPanel("User Guide", page_03),
    tabPanel("Acknowledgements", page_04)
)

# Server ----
server <- function(input, output, session) {
    rvs <- reactiveValues()
    rvs$crd01 <- TRUE
    rvs$crd02 <- TRUE
    rvs$crd03 <- TRUE
    rvs$crd04 <- TRUE
    rvs$crd05 <- TRUE
    rvs$crd06 <- TRUE
    rvs$crd07 <- TRUE
    rvs$crd08 <- TRUE
    rvs$crd09 <- TRUE
    rvs$crd01_microbial_process <- list_contamination_distributions[1]
    rvs$crd02_microbial_process <- list_distributions[1]
    rvs$crd03_microbial_process <- list_distributions[1]
    rvs$crd04_microbial_process <- list_distributions[1]
    rvs$crd05_microbial_process <- list_distributions[1]
    rvs$crd06_microbial_process <- list_distributions[1]
    rvs$crd07_microbial_process <- list_distributions[1]
    rvs$crd08_microbial_process <- list_distributions[1]
    rvs$crd09_microbial_process <- list_distributions[1]
    
    rvs$beep <- TRUE
    
    rvs$alpha <- 0.267
    rvs$beta <- 229.2928
    # rvs$field_mass <- 160000
    rvs$serving_mass <- 85
    
    summary_table_counts <- data.frame(
        scenario = c(
            NA
        ),
        overall = c(
            NA
        ),
        x_bt_1_10 = c(
            NA
        )
    )
    colnames(summary_table_counts) <- c("Scenario Name", "Overall Risk of Positive Test (1 in ...)", "Count of Lots with Risk of >1 in 10 Positives")
    output$summary_table_counts <- renderTable(summary_table_counts, rownames = FALSE)
    
    output$summary_text <- NULL
    
    lot_table_counts <- data.frame(
        overall = c(
            NA
        ),
        x_bt_1_10 = c(
            NA
        ),
        x_bt_10_100 = c(
            NA
        ),
        x_bt_100_1000 = c(
            NA
        ),
        x_bt_1000_10000 = c(
            NA
        ),
        x_gt_10000 = c(
            NA
        ),
        x_et_0 = c(
            NA
        ),
        total = c(
            NA
        )
    )
    colnames(lot_table_counts) <- c("Overall Risk of Positive Test (1 in ...)", ">1 in 10", ">1 in 100", ">1 in 1,000", ">1 in 10,000", "< 1 in 10,000", "Non-contaminated/Rejected", "Total")
    output$lot_table_counts <- renderTable(lot_table_counts, rownames = FALSE)
    rvs$lot_table_counts <- lot_table_counts
    
    lot_table_regulatory_tests <- data.frame(
        x_bt_1_10 = c(
            NA
        ),
        x_bt_10_100 = c(
            NA
        ),
        x_bt_100_1000 = c(
            NA
        ),
        x_bt_1000_10000 = c(
            NA
        ),
        x_gt_10000 = c(
            NA
        ),
        x_et_0 = c(
            NA
        ),
        total = c(
            NA
        )
    )
    colnames(lot_table_regulatory_tests) <- c(">1 in 10", ">1 in 100", ">1 in 1,000", ">1 in 10,000", "< 1 in 10,000", "Non-contaminated/Rejected", "Total Positive")
    output$lot_table_regulatory_tests <- renderTable(lot_table_regulatory_tests, rownames = FALSE)
    rvs$lot_table_regulatory_tests <- lot_table_regulatory_tests
    
    blank_lot_table_counts <- data.frame(
        scenario = c(
            c(rep(NA, 25))  
        ),
        overall = c(
            c(rep(NA, 25))  
        ),
        x_bt_1_10 = c(
            c(rep(NA, 25))  
        ),
        x_bt_10_100 = c(
            c(rep(NA, 25))  
        ),
        x_bt_100_1000 = c(
            c(rep(NA, 25))  
        ),
        x_bt_1000_10000 = c(
            c(rep(NA, 25))  
        ),
        x_gt_10000 = c(
            c(rep(NA, 25))  
        ),
        x_et_0 = c(
            c(rep(NA, 25))  
        ),
        total = c(
            c(rep(NA, 25))  
        )
    )
    colnames(blank_lot_table_counts) <- c("Scenario", "Overall Risk of Positive Test (1 in ...)", ">1 in 10", ">1 in 100", ">1 in 1,000", ">1 in 10,000", "< 1 in 10,000", "Non-contaminated/Rejected", "Total")
    output$results_lot_table_counts <- DT::renderDataTable(blank_lot_table_counts, rownames = FALSE)
    rvs$results_lot_table_counts <- blank_lot_table_counts
    
    blank_lot_table_regulatory_tests <- data.frame(
        scenario = c(
            c(rep(NA, 25))  
        ),
        x_bt_1_10 = c(
            c(rep(NA, 25))  
        ),
        x_bt_10_100 = c(
            c(rep(NA, 25))  
        ),
        x_bt_100_1000 = c(
            c(rep(NA, 25))  
        ),
        x_bt_1000_10000 = c(
            c(rep(NA, 25))  
        ),
        x_gt_10000 = c(
            c(rep(NA, 25))  
        ),
        x_et_0 = c(
            c(rep(NA, 25))  
        ),
        total = c(
            c(rep(NA, 25))  
        )
    )
    colnames(blank_lot_table_regulatory_tests) <- c("Scenario", ">1 in 10", ">1 in 100", ">1 in 1,000", ">1 in 10,000", "< 1 in 10,000", "Non-contaminated/Rejected", "Total")
    output$results_lot_table_regulatory_tests <- DT::renderDataTable(blank_lot_table_regulatory_tests, rownames = FALSE)
    rvs$results_lot_table_regulatory_tests <- blank_lot_table_regulatory_tests
    
    rvs$alert <- FALSE
    
    output$dwnldbtwn_results_counts <- downloadHandler(
        filename = function() {
            paste0(rvs$timestamp, " ", as.character(rvs$scenario_name), "_lot_table_counts.csv")
        },
        content = function(con) {
            write.csv(rvs$lot_table_counts, con, row.names = FALSE)
        }
    )
    
    output$dwnldbtwn_results_regulatory_tests <- downloadHandler(
        filename = function() {
            paste0(rvs$timestamp, " ", as.character(rvs$scenario_name), "_lot_table_regulatory_tests.csv")
        },
        content = function(con) {
            write.csv(rvs$lot_table_regulatory_tests, con, row.names = FALSE)
        }
    )
    
    output$dwnldbtwn_results_comparision_counts <- downloadHandler(
        filename = function() {
            paste0(as.character(paste0(format(Sys.Date()),"-", format(Sys.time(), "%H-%M-%S"))), " comparision_table_counts.csv")
        },
        content = function(con) {
            write.csv(na.omit(rvs$results_lot_table_counts), con, row.names = FALSE)
        }
    )
    
    output$dwnldbtwn_results_comparision_regulatory_tests <- downloadHandler(
        filename = function() {
            paste0(as.character(paste0(format(Sys.Date()),"-", format(Sys.time(), "%H-%M-%S"))), " comparision_table_regulatory_tests.csv")
        },
        content = function(con) {
            write.csv(na.omit(rvs$results_lot_table_regulatory_tests), con, row.names = FALSE)
        }
    )
    
    toListenSelectinScenario <- reactive({
        list(
            input$slctin_scenario
        )
    })
    
    toListenTextinScenarioName <- reactive({
        list(
            input$txtin_scenario_name
        )
    })
    
    toListenNumericinFieldMassData <- reactive({
        list(
            input$nmcin_field_mass,
            input$slctin_field_mass_unit
        )
    })
    
    toListenNumericinFieldMassUnit <- reactive({
        list(
            input$slctin_field_mass_unit
        )
    })
    
    toListenNumericinServingMass <- reactive({
        list(
            input$nmcin_serving_mass
        )
    })
    
    toListenNumericinAlpha <- reactive({
        list(
            input$nmcin_alpha
        )
    })
    
    toListenNumericinBeta <- reactive({
        list(
            input$nmcin_beta
        )
    })
    
    toListenNumericinIterations <- reactive({
        list(
            input$nmcin_n_iterations
        )
    })
    
    toListenSelectinSaveSlot <- reactive({
        list(
            input$slctin_save_slot
        )
    })
    
    toListenCheckboxBeep <- reactive({
        list(
            input$ckbxin_beep
        )
    })
    
    toListenCheckboxCard01 <- reactive({
        list(
            input$ckbxin_crd01
        )
    })
    toListenTextinCard01ModelStep <- reactive({
        list(
            input$txtin_crd01_model_step
        )
    })
    toListenSelectinCard01MicrobialProcess <- reactive({
        list(
            input$slctin_crd01_microbial_process
        )
    })
    toListenSelectinCard01Distribution <- reactive({
        list(
            input$slctin_crd01_distribution
        )
    })
    toListenSelectinCard01ProcessAndDist <- reactive({
        list(
            input$slctin_crd01_microbial_process,
            input$slctin_crd01_distribution
        )
    })
    toListenNumberinCard01POccurrence <- reactive({
        list(
            input$nmcin_crd01_p_occurrence
        )
    })
    toListenSelectinCard01ContaminationBasis <- reactive({
        list(
            input$slctin_crd01_contamination_basis
        )
    })
    toListenNumberinCard01Parameter1 <- reactive({
        list(
            input$nmcin_crd01_parameter_1
        )
    })
    toListenNumberinCard01Parameter2 <- reactive({
        list(
            input$nmcin_crd01_parameter_2
        )
    })
    toListenNumberinCard01Parameter3 <- reactive({
        list(
            input$nmcin_crd01_parameter_3
        )
    })
    
    toListenCheckboxCard02 <- reactive({
        list(
            input$ckbxin_crd02
        )
    })
    toListenTextinCard02ModelStep <- reactive({
        list(
            input$txtin_crd02_model_step
        )
    })
    toListenSelectinCard02MicrobialProcess <- reactive({
        list(
            input$slctin_crd02_microbial_process
        )
    })
    toListenSelectinCard02Distribution <- reactive({
        list(
            input$slctin_crd02_distribution
        )
    })
    toListenSelectinCard02ProcessAndDist <- reactive({
        list(
            input$slctin_crd02_microbial_process,
            input$slctin_crd02_distribution
        )
    })
    toListenNumberinCard02POccurrence <- reactive({
        list(
            input$nmcin_crd02_p_occurrence
        )
    })
    toListenSelectinCard02ContaminationBasis <- reactive({
        list(
            input$slctin_crd02_contamination_basis
        )
    })
    toListenNumberinCard02Parameter1 <- reactive({
        list(
            input$nmcin_crd02_parameter_1
        )
    })
    toListenNumberinCard02Parameter2 <- reactive({
        list(
            input$nmcin_crd02_parameter_2
        )
    })
    toListenNumberinCard02Parameter3 <- reactive({
        list(
            input$nmcin_crd02_parameter_3
        )
    })
    
    toListenCheckboxCard03 <- reactive({
        list(
            input$ckbxin_crd03
        )
    })
    toListenTextinCard03ModelStep <- reactive({
        list(
            input$txtin_crd03_model_step
        )
    })
    toListenSelectinCard03MicrobialProcess <- reactive({
        list(
            input$slctin_crd03_microbial_process
        )
    })
    toListenSelectinCard03Distribution <- reactive({
        list(
            input$slctin_crd03_distribution
        )
    })
    toListenSelectinCard03ProcessAndDist <- reactive({
        list(
            input$slctin_crd03_microbial_process,
            input$slctin_crd03_distribution
        )
    })
    toListenNumberinCard03POccurrence <- reactive({
        list(
            input$nmcin_crd03_p_occurrence
        )
    })
    toListenSelectinCard03ContaminationBasis <- reactive({
        list(
            input$slctin_crd03_contamination_basis
        )
    })
    toListenNumberinCard03Parameter1 <- reactive({
        list(
            input$nmcin_crd03_parameter_1
        )
    })
    toListenNumberinCard03Parameter2 <- reactive({
        list(
            input$nmcin_crd03_parameter_2
        )
    })
    toListenNumberinCard03Parameter3 <- reactive({
        list(
            input$nmcin_crd03_parameter_3
        )
    })
    
    toListenCheckboxCard04 <- reactive({
        list(
            input$ckbxin_crd04
        )
    })
    toListenTextinCard04ModelStep <- reactive({
        list(
            input$txtin_crd04_model_step
        )
    })
    toListenSelectinCard04MicrobialProcess <- reactive({
        list(
            input$slctin_crd04_microbial_process
        )
    })
    toListenSelectinCard04Distribution <- reactive({
        list(
            input$slctin_crd04_distribution
        )
    })
    toListenSelectinCard04ProcessAndDist <- reactive({
        list(
            input$slctin_crd04_microbial_process,
            input$slctin_crd04_distribution
        )
    })
    toListenNumberinCard04POccurrence <- reactive({
        list(
            input$nmcin_crd04_p_occurrence
        )
    })
    toListenSelectinCard04ContaminationBasis <- reactive({
        list(
            input$slctin_crd04_contamination_basis
        )
    })
    toListenNumberinCard04Parameter1 <- reactive({
        list(
            input$nmcin_crd04_parameter_1
        )
    })
    toListenNumberinCard04Parameter2 <- reactive({
        list(
            input$nmcin_crd04_parameter_2
        )
    })
    toListenNumberinCard04Parameter3 <- reactive({
        list(
            input$nmcin_crd04_parameter_3
        )
    })
    
    toListenCheckboxCard05 <- reactive({
        list(
            input$ckbxin_crd05
        )
    })
    toListenTextinCard05ModelStep <- reactive({
        list(
            input$txtin_crd05_model_step
        )
    })
    toListenSelectinCard05MicrobialProcess <- reactive({
        list(
            input$slctin_crd05_microbial_process
        )
    })
    toListenSelectinCard05Distribution <- reactive({
        list(
            input$slctin_crd05_distribution
        )
    })
    toListenSelectinCard05ProcessAndDist <- reactive({
        list(
            input$slctin_crd05_microbial_process,
            input$slctin_crd05_distribution
        )
    })
    toListenNumberinCard05POccurrence <- reactive({
        list(
            input$nmcin_crd05_p_occurrence
        )
    })
    toListenSelectinCard05ContaminationBasis <- reactive({
        list(
            input$slctin_crd05_contamination_basis
        )
    })
    toListenNumberinCard05Parameter1 <- reactive({
        list(
            input$nmcin_crd05_parameter_1
        )
    })
    toListenNumberinCard05Parameter2 <- reactive({
        list(
            input$nmcin_crd05_parameter_2
        )
    })
    toListenNumberinCard05Parameter3 <- reactive({
        list(
            input$nmcin_crd05_parameter_3
        )
    })
    
    toListenCheckboxCard06 <- reactive({
        list(
            input$ckbxin_crd06
        )
    })
    toListenTextinCard06ModelStep <- reactive({
        list(
            input$txtin_crd06_model_step
        )
    })
    toListenSelectinCard06MicrobialProcess <- reactive({
        list(
            input$slctin_crd06_microbial_process
        )
    })
    toListenSelectinCard06Distribution <- reactive({
        list(
            input$slctin_crd06_distribution
        )
    })
    toListenSelectinCard06ProcessAndDist <- reactive({
        list(
            input$slctin_crd06_microbial_process,
            input$slctin_crd06_distribution
        )
    })
    toListenNumberinCard06POccurrence <- reactive({
        list(
            input$nmcin_crd06_p_occurrence
        )
    })
    toListenSelectinCard06ContaminationBasis <- reactive({
        list(
            input$slctin_crd06_contamination_basis
        )
    })
    toListenNumberinCard06Parameter1 <- reactive({
        list(
            input$nmcin_crd06_parameter_1
        )
    })
    toListenNumberinCard06Parameter2 <- reactive({
        list(
            input$nmcin_crd06_parameter_2
        )
    })
    toListenNumberinCard06Parameter3 <- reactive({
        list(
            input$nmcin_crd06_parameter_3
        )
    })
    
    toListenCheckboxCard07 <- reactive({
        list(
            input$ckbxin_crd07
        )
    })
    toListenTextinCard07ModelStep <- reactive({
        list(
            input$txtin_crd07_model_step
        )
    })
    toListenSelectinCard07MicrobialProcess <- reactive({
        list(
            input$slctin_crd07_microbial_process
        )
    })
    toListenSelectinCard07Distribution <- reactive({
        list(
            input$slctin_crd07_distribution
        )
    })
    toListenSelectinCard07ProcessAndDist <- reactive({
        list(
            input$slctin_crd07_microbial_process,
            input$slctin_crd07_distribution
        )
    })
    toListenNumberinCard07POccurrence <- reactive({
        list(
            input$nmcin_crd07_p_occurrence
        )
    })
    toListenSelectinCard07ContaminationBasis <- reactive({
        list(
            input$slctin_crd07_contamination_basis
        )
    })
    toListenNumberinCard07Parameter1 <- reactive({
        list(
            input$nmcin_crd07_parameter_1
        )
    })
    toListenNumberinCard07Parameter2 <- reactive({
        list(
            input$nmcin_crd07_parameter_2
        )
    })
    toListenNumberinCard07Parameter3 <- reactive({
        list(
            input$nmcin_crd07_parameter_3
        )
    })
    
    toListenCheckboxCard08 <- reactive({
        list(
            input$ckbxin_crd08
        )
    })
    toListenTextinCard08ModelStep <- reactive({
        list(
            input$txtin_crd08_model_step
        )
    })
    toListenSelectinCard08MicrobialProcess <- reactive({
        list(
            input$slctin_crd08_microbial_process
        )
    })
    toListenSelectinCard08Distribution <- reactive({
        list(
            input$slctin_crd08_distribution
        )
    })
    toListenSelectinCard08ProcessAndDist <- reactive({
        list(
            input$slctin_crd08_microbial_process,
            input$slctin_crd08_distribution
        )
    })
    toListenNumberinCard08POccurrence <- reactive({
        list(
            input$nmcin_crd08_p_occurrence
        )
    })
    toListenSelectinCard08ContaminationBasis <- reactive({
        list(
            input$slctin_crd08_contamination_basis
        )
    })
    toListenNumberinCard08Parameter1 <- reactive({
        list(
            input$nmcin_crd08_parameter_1
        )
    })
    toListenNumberinCard08Parameter2 <- reactive({
        list(
            input$nmcin_crd08_parameter_2
        )
    })
    toListenNumberinCard08Parameter3 <- reactive({
        list(
            input$nmcin_crd08_parameter_3
        )
    })
    
    toListenCheckboxCard09 <- reactive({
        list(
            input$ckbxin_crd09
        )
    })
    toListenTextinCard09ModelStep <- reactive({
        list(
            input$txtin_crd09_model_step
        )
    })
    toListenSelectinCard09MicrobialProcess <- reactive({
        list(
            input$slctin_crd09_microbial_process
        )
    })
    toListenSelectinCard09Distribution <- reactive({
        list(
            input$slctin_crd09_distribution
        )
    })
    toListenSelectinCard09ProcessAndDist <- reactive({
        list(
            input$slctin_crd09_microbial_process,
            input$slctin_crd09_distribution
        )
    })
    toListenNumberinCard09POccurrence <- reactive({
        list(
            input$nmcin_crd09_p_occurrence
        )
    })
    toListenSelectinCard09ContaminationBasis <- reactive({
        list(
            input$slctin_crd09_contamination_basis
        )
    })
    toListenNumberinCard09Parameter1 <- reactive({
        list(
            input$nmcin_crd09_parameter_1
        )
    })
    toListenNumberinCard09Parameter2 <- reactive({
        list(
            input$nmcin_crd09_parameter_2
        )
    })
    toListenNumberinCard09Parameter3 <- reactive({
        list(
            input$nmcin_crd09_parameter_3
        )
    })
    
    output$img_uiuc <- renderImage(
        {
            list(
                src = normalizePath(file.path('./www','University-Wordmark-Full-Color-RGB-TM.png')),
                contentType = 'image/png',
                alt = "Block I with University of Illinois Urbana-Champaign wordmark to the right",
                width = "250"
            )
        },
        deleteFile = FALSE
    )
    
    output$img_cornell <- renderImage(
        {
            list(
                src = normalizePath(file.path('./www','bold_cornell_logo_pms187_red.png')),
                contentType = 'image/png',
                alt = "Cornell University Logo",
                width = "250"
            )
        },
        deleteFile = FALSE
    )
    
    reset_steps <- function() {
        updateTextInput(session, "txtin_scenario_name", value = "")
        # Card 01
        updateCheckboxInput(session, "ckbxin_crd01", value = FALSE)
        updateTextInput(session, "txtin_crd01_model_step", value = "")
        updateSelectInput(session, "slctin_crd01_microbial_process", selected = list_processes[1])
        updateSelectInput(session, "slctin_crd01_distribution", selected = list_distributions[1])
        updateSelectInput(session, "slctin_crd01_contamination_basis", selected = list_contamination_types[2])
        updateNumericInput(session, "nmcin_crd01_p_occurrence", value = 100)
        updateNumericInput(session, "nmcin_crd01_parameter_1", value = 0)
        updateNumericInput(session, "nmcin_crd01_parameter_2", value = 0)
        updateNumericInput(session, "nmcin_crd01_parameter_3", value = 0)
        # Card 02
        updateCheckboxInput(session, "ckbxin_crd02", value = FALSE)
        updateTextInput(session, "txtin_crd02_model_step", value = "")
        updateSelectInput(session, "slctin_crd02_microbial_process", selected = list_processes[2])
        updateSelectInput(session, "slctin_crd02_distribution", selected = list_distributions[1])
        updateSelectInput(session, "slctin_crd02_contamination_basis", selected = list_contamination_types[2])
        updateNumericInput(session, "nmcin_crd02_p_occurrence",value = 1)
        updateNumericInput(session, "nmcin_crd02_parameter_1", value = 0)
        updateNumericInput(session, "nmcin_crd02_parameter_2", value = 0)
        updateNumericInput(session, "nmcin_crd02_parameter_3", value = 0)
        # Card 03
        updateCheckboxInput(session, "ckbxin_crd03", value = FALSE)
        updateTextInput(session, "txtin_crd03_model_step", value = "")
        updateSelectInput(session, "slctin_crd03_microbial_process", selected = list_processes[2])
        updateSelectInput(session, "slctin_crd03_distribution", selected = list_distributions[1])
        updateSelectInput(session, "slctin_crd03_contamination_basis", selected = list_contamination_types[2])
        updateNumericInput(session, "nmcin_crd03_p_occurrence", value = 100)
        updateNumericInput(session, "nmcin_crd03_parameter_1", value = 0)
        updateNumericInput(session, "nmcin_crd03_parameter_2", value = 0)
        updateNumericInput(session, "nmcin_crd03_parameter_3", value = 0)
        # Card 04
        updateCheckboxInput(session, "ckbxin_crd04", value = FALSE)
        updateTextInput(session, "txtin_crd04_model_step", value = "")
        updateSelectInput(session, "slctin_crd04_microbial_process", selected = list_processes[2])
        updateSelectInput(session, "slctin_crd04_distribution", selected = list_distributions[1])
        updateSelectInput(session, "slctin_crd04_contamination_basis", selected = list_contamination_types[2])
        updateNumericInput(session, "nmcin_crd04_p_occurrence", value = 100)
        updateNumericInput(session, "nmcin_crd04_parameter_1", value = 0)
        updateNumericInput(session, "nmcin_crd04_parameter_2", value = 0)
        updateNumericInput(session, "nmcin_crd04_parameter_3", value = 0)
        # Card 05
        updateCheckboxInput(session, "ckbxin_crd05", value = FALSE)
        updateTextInput(session, "txtin_crd05_model_step", value = "")
        updateSelectInput(session, "slctin_crd05_microbial_process", selected = list_processes[2])
        updateSelectInput(session, "slctin_crd05_distribution", selected = list_distributions[1])
        updateSelectInput(session, "slctin_crd05_contamination_basis", selected = list_contamination_types[2])
        updateNumericInput(session, "nmcin_crd05_p_occurrence", value = 100)
        updateNumericInput(session, "nmcin_crd05_parameter_1", value = 0)
        updateNumericInput(session, "nmcin_crd05_parameter_2", value = 0)
        updateNumericInput(session, "nmcin_crd05_parameter_3", value = 0)
        # Card 06
        updateCheckboxInput(session, "ckbxin_crd06", value = FALSE)
        updateTextInput(session, "txtin_crd06_model_step", value = "")
        updateSelectInput(session, "slctin_crd06_microbial_process", selected = list_processes[2])
        updateSelectInput(session, "slctin_crd06_distribution", selected = list_distributions[1])
        updateSelectInput(session, "slctin_crd06_contamination_basis", selected = list_contamination_types[2])
        updateNumericInput(session, "nmcin_crd06_p_occurrence", value = 100)
        updateNumericInput(session, "nmcin_crd06_parameter_1", value = 0)
        updateNumericInput(session, "nmcin_crd06_parameter_2", value = 0)
        updateNumericInput(session, "nmcin_crd06_parameter_3", value = 0)
        # Card 07
        updateCheckboxInput(session, "ckbxin_crd07", value = FALSE)
        updateTextInput(session, "txtin_crd07_model_step", value = "")
        updateSelectInput(session, "slctin_crd07_microbial_process", selected = list_processes[2])
        updateSelectInput(session, "slctin_crd07_distribution", selected = list_distributions[1])
        updateSelectInput(session, "slctin_crd07_contamination_basis", selected = list_contamination_types[2])
        updateNumericInput(session, "nmcin_crd07_p_occurrence", value = 100)
        updateNumericInput(session, "nmcin_crd07_parameter_1", value = 0)
        updateNumericInput(session, "nmcin_crd07_parameter_2", value = 0)
        updateNumericInput(session, "nmcin_crd07_parameter_3", value = 0)
        # Card 08
        updateCheckboxInput(session, "ckbxin_crd08", value = FALSE)
        updateTextInput(session, "txtin_crd08_model_step", value = "")
        updateSelectInput(session, "slctin_crd08_microbial_process", selected = list_processes[2])
        updateSelectInput(session, "slctin_crd08_distribution", selected = list_distributions[1])
        updateSelectInput(session, "slctin_crd08_contamination_basis", selected = list_contamination_types[2])
        updateNumericInput(session, "nmcin_crd08_p_occurrence", value = 100)
        updateNumericInput(session, "nmcin_crd08_parameter_1", value = 0)
        updateNumericInput(session, "nmcin_crd08_parameter_2", value = 0)
        updateNumericInput(session, "nmcin_crd08_parameter_3", value = 0)
        # Card 09
        updateCheckboxInput(session, "ckbxin_crd09", value = FALSE)
        updateTextInput(session, "txtin_crd09_model_step", value = "")
        updateSelectInput(session, "slctin_crd09_microbial_process", selected = list_processes[2])
        updateSelectInput(session, "slctin_crd09_distribution", selected = list_distributions[1])
        updateSelectInput(session, "slctin_crd09_contamination_basis", selected = list_contamination_types[2])
        updateNumericInput(session, "nmcin_crd09_p_occurrence", value = 100)
        updateNumericInput(session, "nmcin_crd09_parameter_1", value = 0)
        updateNumericInput(session, "nmcin_crd09_parameter_2", value = 0)
        updateNumericInput(session, "nmcin_crd09_parameter_3", value = 0)
    }
    
    updateSteps <- function() {
        steps <- data.frame(
            model_step = character(),
            microbial_process = character(),
            distribution = character(),
            contamination_basis = character(),
            freq = numeric(),
            parameter_1 = numeric(),
            parameter_2 = numeric(),
            parameter_3 = numeric(),
            stringsAsFactors = FALSE
        )
        step_names <- colnames(steps)
        if (rvs$crd01 == TRUE) {
            new_step <- data.frame(
                as.character(rvs$crd01_model_step),
                as.character(rvs$crd01_microbial_process),
                as.character(rvs$crd01_distribution),
                as.character(rvs$crd01_contamination_basis),
                as.numeric(rvs$crd01_p_occurrence)/100,
                as.numeric(rvs$crd01_parameter_1),
                as.numeric(rvs$crd01_parameter_2),
                as.numeric(rvs$crd01_parameter_3),
                stringsAsFactors = FALSE
            )
            colnames(new_step) <- step_names
            steps <- rbind(steps, new_step)
            colnames(steps) <- step_names
        }
        if (rvs$crd02 == TRUE) {
            new_step <- data.frame(
                as.character(rvs$crd02_model_step),
                as.character(rvs$crd02_microbial_process),
                as.character(rvs$crd02_distribution),
                as.character(rvs$crd02_contamination_basis),
                as.numeric(rvs$crd02_p_occurrence)/100,
                as.numeric(rvs$crd02_parameter_1),
                as.numeric(rvs$crd02_parameter_2),
                as.numeric(rvs$crd02_parameter_3),
                stringsAsFactors = FALSE
            )
            colnames(new_step) <- step_names
            steps <- rbind(steps, new_step)
            colnames(steps) <- step_names
        }
        if (rvs$crd03 == TRUE) {
            new_step <- data.frame(
                as.character(rvs$crd03_model_step),
                as.character(rvs$crd03_microbial_process),
                as.character(rvs$crd03_distribution),
                as.character(rvs$crd03_contamination_basis),
                as.numeric(rvs$crd03_p_occurrence)/100,
                as.numeric(rvs$crd03_parameter_1),
                as.numeric(rvs$crd03_parameter_2),
                as.numeric(rvs$crd03_parameter_3),
                stringsAsFactors = FALSE
            )
            colnames(new_step) <- step_names
            steps <- rbind(steps, new_step)
            colnames(steps) <- step_names
        }
        if (rvs$crd04 == TRUE) {
            new_step <- data.frame(
                as.character(rvs$crd04_model_step),
                as.character(rvs$crd04_microbial_process),
                as.character(rvs$crd04_distribution),
                as.character(rvs$crd04_contamination_basis),
                as.numeric(rvs$crd04_p_occurrence)/100,
                as.numeric(rvs$crd04_parameter_1),
                as.numeric(rvs$crd04_parameter_2),
                as.numeric(rvs$crd04_parameter_3),
                stringsAsFactors = FALSE
            )
            colnames(new_step) <- step_names
            steps <- rbind(steps, new_step)
            colnames(steps) <- step_names
        }
        if (rvs$crd05 == TRUE) {
            new_step <- data.frame(
                as.character(rvs$crd05_model_step),
                as.character(rvs$crd05_microbial_process),
                as.character(rvs$crd05_distribution),
                as.character(rvs$crd05_contamination_basis),
                as.numeric(rvs$crd05_p_occurrence)/100,
                as.numeric(rvs$crd05_parameter_1),
                as.numeric(rvs$crd05_parameter_2),
                as.numeric(rvs$crd05_parameter_3),
                stringsAsFactors = FALSE
            )
            colnames(new_step) <- step_names
            steps <- rbind(steps, new_step)
            colnames(steps) <- step_names
        }
        if (rvs$crd06 == TRUE) {
            new_step <- data.frame(
                as.character(rvs$crd06_model_step),
                as.character(rvs$crd06_microbial_process),
                as.character(rvs$crd06_distribution),
                as.character(rvs$crd06_contamination_basis),
                as.numeric(rvs$crd06_p_occurrence)/100,
                as.numeric(rvs$crd06_parameter_1),
                as.numeric(rvs$crd06_parameter_2),
                as.numeric(rvs$crd06_parameter_3),
                stringsAsFactors = FALSE
            )
            colnames(new_step) <- step_names
            steps <- rbind(steps, new_step)
            colnames(steps) <- step_names
        }
        if (rvs$crd07 == TRUE) {
            new_step <- data.frame(
                as.character(rvs$crd07_model_step),
                as.character(rvs$crd07_microbial_process),
                as.character(rvs$crd07_distribution),
                as.character(rvs$crd07_contamination_basis),
                as.numeric(rvs$crd07_p_occurrence)/100,
                as.numeric(rvs$crd07_parameter_1),
                as.numeric(rvs$crd07_parameter_2),
                as.numeric(rvs$crd07_parameter_3),
                stringsAsFactors = FALSE
            )
            colnames(new_step) <- step_names
            steps <- rbind(steps, new_step)
            colnames(steps) <- step_names
        }
        if (rvs$crd08 == TRUE) {
            new_step <- data.frame(
                as.character(rvs$crd08_model_step),
                as.character(rvs$crd08_microbial_process),
                as.character(rvs$crd08_distribution),
                as.character(rvs$crd08_contamination_basis),
                as.numeric(rvs$crd08_p_occurrence)/100,
                as.numeric(rvs$crd08_parameter_1),
                as.numeric(rvs$crd08_parameter_2),
                as.numeric(rvs$crd08_parameter_3),
                stringsAsFactors = FALSE
            )
            colnames(new_step) <- step_names
            steps <- rbind(steps, new_step)
            colnames(steps) <- step_names
        }
        if (rvs$crd09 == TRUE) {
            new_step <- data.frame(
                as.character(rvs$crd09_model_step),
                as.character(rvs$crd09_microbial_process),
                as.character(rvs$crd09_distribution),
                as.character(rvs$crd09_contamination_basis),
                as.numeric(rvs$crd09_p_occurrence)/100,
                as.numeric(rvs$crd09_parameter_1),
                as.numeric(rvs$crd09_parameter_2),
                as.numeric(rvs$crd09_parameter_3),
                stringsAsFactors = FALSE
            )
            colnames(new_step) <- step_names
            steps <- rbind(steps, new_step)
            colnames(steps) <- step_names
        }
        if (any(steps$microbial_process == "Increase/Reduction")) {
            steps$contamination_basis[which(steps$microbial_process == "Increase/Reduction")] <- NA
        }
        if (any(steps$microbial_process == "Risk Output Test" | steps$microbial_process == "Product Test")) {
            steps$distribution[which(steps$microbial_process == "Risk Output Test" | steps$microbial_process == "Product Test")] <- NA
            steps$contamination_basis[which(steps$microbial_process == "Risk Output Test" | steps$microbial_process == "Product Test")] <- NA
        }
        if (any(steps$distribution == "Normal (log10)" | steps$distribution == "Uniform (log10)" | steps$distribution == "Normal (linear)" | steps$distribution == "Uniform (linear)")) {
            steps$parameter_3[which(steps$distribution == "Normal (log10)" | steps$distribution == "Uniform (log10)" | steps$distribution == "Normal (linear)" | steps$distribution == "Uniform (linear)")] <- NA
        }
        return(steps)
    }
    
    observeEvent(input$btn_reset, {
        reset_steps()
    })
    
    observeEvent(input$btn_run, {
        n_iterations <- as.numeric(rvs$n_iterations)
        if (n_iterations < as.numeric(detectCores())) {
            cl <- NA
            parallel_mode <- FALSE
        } else {
            cl <- makeCluster(as.numeric(detectCores()))
            parallel_mode <- TRUE
        }
        # parallel_mode <- FALSE # Parallel Override for test purposes
        n_simulations <- 1
        n_cores <- as.numeric(rvs$n_cores)
        field_mass <- as.numeric(rvs$field_mass)
        serving_mass_g <- as.numeric(rvs$serving_mass)
        alpha <- as.numeric(rvs$alpha)
        beta <- as.numeric(rvs$beta)
        unit <- as.character(rvs$field_mass_unit)
        
        steps <- updateSteps()
        steps$carried_contamination <- NA
        steps$carried_contamination[1] <- 0
        lhs_width <- lhs_matrix_width(steps)
        
        if (length(which(steps$microbial_process == "Risk Output Test")) == 0) {
            showNotification(id = "ntfcn_error_regulatory_testing", "Error! Risk Output Test step is disabled or missing.", duration = 10, closeButton = TRUE, type = "error")
        } else if (length(which(steps$microbial_process == "Risk Output Test")) > 1) {
            showNotification(id = "ntfcn_error_regulatory_testing", "Error! Only 1 Risk Output Test is currently supported.", duration = 10, closeButton = TRUE, type = "error")
        } else {
            if ((steps$microbial_process[nrow(steps)] != "Risk Output Test")) {
                shinyalert(
                    title = "Alert",
                    text = "One or more of the steps you have entered are below the product test at retail.\nPlease confirm your entries.",
                    size = "s", 
                    closeOnEsc = TRUE,
                    closeOnClickOutside = FALSE,
                    html = FALSE,
                    type = "warning",
                    showConfirmButton = TRUE,
                    showCancelButton = TRUE,
                    confirmButtonText = "Continue",
                    confirmButtonCol = "#AEDEF4",
                    cancelButtonText = "Cancel",
                    timer = 0,
                    imageUrl = "",
                    animation = TRUE,
                    callbackR = function(x) {
                        if (x == TRUE) {
                            set.seed(42, "Mersenne-Twister")
                            lhs_matrix <- lhs::randomLHS(n_iterations, lhs_width)
                            time_started <- as.character(format(Sys.time(), "%X"))
                            removeNotification(id = "ntfcn_finished")
                            removeNotification(id = "ntfcn_warning_product_testing")
                            if (length(which(steps$microbial_process == "Product Test")) == 0) {
                                showNotification(id = "ntfcn_warning_product_testing", paste0("Warning! Product Test step is disabled or missing."), duration = NULL, closeButton = TRUE, type = "warning")
                            }
                            showNotification(id = "ntfcn_running", paste0("Running! Please Wait.\nStarted at ", time_started), duration = NULL, closeButton = FALSE, type = "message")
                            if (parallel_mode == TRUE) {
                                clusterExport(
                                    cl,
                                    varlist = c(
                                        "n_iterations",
                                        "n_simulations",
                                        "n_cores",
                                        "steps",
                                        "lhs_matrix",
                                        "field_mass",
                                        "serving_mass_g",
                                        "alpha",
                                        "beta",
                                        "unit"
                                    ),
                                    envir = environment()
                                )
                            }
                            rvs$model_results <- scrm_clustered(
                                n_iterations,
                                1,
                                n_cores,
                                steps,
                                lhs_matrix,
                                field_mass,
                                serving_mass_g,
                                alpha,
                                beta,
                                cl,
                                parallel_mode,
                                unit
                            )
                            rvs$model_results$detected_by_retail <-  as.numeric(rvs$model_results$p_positive_regulatory_test)>1 & as.numeric(rvs$model_results$testing_outcome_product_sum>0)
                            rvs$model_results$detected_by_fp_testing <-  as.numeric(rvs$model_results$p_positive_regulatory_test)>1 & as.numeric(rvs$model_results$testing_outcome_regulatory_sum>0)
                            rvs$p_positive_regulatory <- p_positive_regulatory_summarizer(rvs$model_results)
                            
                            rvs$n_positives <- sum(as.numeric(rvs$model_results$p_positive_regulatory_test))
                            rvs$risk_of_positive_bins <- risk_binner(rvs$model_results, rvs$n_iterations)
                            rvs$risk_of_positive_bins_new <- risk_binner_02(rvs$model_results)
                            
                            if (1/rvs$p_positive_regulatory$mean > 100) {
                                rvs$overall <- as.integer(round(1/rvs$p_positive_regulatory$mean, digits = -2))
                            } else {
                                rvs$overall <- as.integer(round(1/rvs$p_positive_regulatory$mean, digits = 0))
                            }
                            
                            summary_table_counts <- data.frame(
                                scenario = c(
                                    as.character(rvs$scenario_name)
                                ),
                                overall = c(
                                    rvs$overall
                                ),
                                x_bt_1_10 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$counts[2]))
                                )
                            )
                            colnames(summary_table_counts) <- c("Scenario Name", "Overall Risk of Positive Test (1 in ...)", "Count of Lots with Risk of >1 in 10 Positives")
                            output$summary_table_counts <- renderTable(summary_table_counts, rownames = FALSE)
                            
                            scenario_name <- rvs$scenario_name
                            overall <- rvs$overall
                            highest_risk_count <- as.integer(rvs$risk_of_positive_bins_new$counts[2])
                            risk_output_step <- steps$model_step[which(steps$microbial_process == "Risk Output Test")]
                            
                            output$summary_text <- renderUI(paste0("The overall risk of a ", risk_output_step, " sample testing positive under the ", scenario_name, " scenario was approximately 1 in ", overall, ", and resulted in ", highest_risk_count, " highest-risk lots, meaning those with greater than 1 in 10 change of being Positive."))
                            
                            lot_table_counts <- data.frame(
                                overall = c(
                                    rvs$overall
                                ),
                                x_bt_1_10 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$counts[2]))
                                ),
                                x_bt_10_100 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$counts[3]))
                                ),
                                x_bt_100_1000 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$counts[4]))
                                ),
                                x_bt_1000_10000 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$counts[5]))
                                ),
                                x_gt_10000 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$counts[6]))
                                ),
                                x_et_0 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$counts[7]))
                                ),
                                total = c(
                                    sum(
                                        as.integer(rvs$risk_of_positive_bins_new$counts[2]),
                                        as.integer(rvs$risk_of_positive_bins_new$counts[3]),
                                        as.integer(rvs$risk_of_positive_bins_new$counts[4]),
                                        as.integer(rvs$risk_of_positive_bins_new$counts[5]),
                                        as.integer(rvs$risk_of_positive_bins_new$counts[6]),
                                        as.integer(rvs$risk_of_positive_bins_new$counts[7])
                                    )
                                )
                            )
                            colnames(lot_table_counts) <- c("Overall Risk of Positive Test (1 in ...)", ">1 in 10", ">1 in 100", ">1 in 1,000", ">1 in 10,000", "< 1 in 10,000", "Non-contaminated/Rejected", "Total")
                            output$lot_table_counts <- renderTable(lot_table_counts, rownames = FALSE)
                            rvs$lot_table_counts <- lot_table_counts
                            
                            lot_table_regulatory_tests <- data.frame(
                                x_bt_1_10 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[2]))
                                ),
                                x_bt_10_100 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[3]))
                                ),
                                x_bt_100_1000 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[4]))
                                ),
                                x_bt_1000_10000 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[5]))
                                ),
                                x_gt_10000 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[6]))
                                ),
                                x_et_0 = c(
                                    (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[7]))
                                ),
                                total = c(
                                    sum(
                                        as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[2]),
                                        as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[3]),
                                        as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[4]),
                                        as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[5]),
                                        as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[6])
                                    )
                                )
                            )
                            colnames(lot_table_regulatory_tests) <- c(">1 in 10", ">1 in 100", ">1 in 1,000", ">1 in 10,000", "< 1 in 10,000", "Non-contaminated/Rejected", "Total Positive")
                            output$lot_table_regulatory_tests <- renderTable(lot_table_regulatory_tests, rownames = FALSE)
                            rvs$lot_table_regulatory_tests <- lot_table_regulatory_tests
                            
                            time_finished <- as.character(format(Sys.time(), "%X"))
                            rvs$timestamp <- as.character(paste0(format(Sys.Date()),"-", format(Sys.time(), "%H-%M-%S")))
                            removeNotification(id = "ntfcn_running")
                            showNotification(id = "ntfcn_finished", paste0("Finished at ", time_finished, " (Started at ", time_started, ")"), duration = NULL, closeButton = TRUE, type = "message")
                            if (rvs$beep == TRUE) {
                                beep()
                            }
                            # write.csv(rvs$summary_table_counts, "summary_table_counts.csv", row.names = FALSE)
                        }
                    }
                )
            } else {
                set.seed(42, "Mersenne-Twister")
                lhs_matrix <- lhs::randomLHS(n_iterations, lhs_width)
                time_started <- as.character(format(Sys.time(), "%X"))
                removeNotification(id = "ntfcn_warning_product_testing")
                if (length(which(steps$microbial_process == "Product Test")) == 0) {
                    showNotification(id = "ntfcn_warning_product_testing", paste0("Warning! Product Test step is disabled or missing."), duration = NULL, closeButton = TRUE, type = "warning")
                }
                showNotification(id = "ntfcn_running", paste0("Running! Please Wait.\nStarted at ", time_started), duration = NULL, closeButton = FALSE, type = "message")
                if (parallel_mode == TRUE) {
                    clusterExport(
                        cl,
                        varlist = c(
                            "n_iterations",
                            "n_simulations",
                            "n_cores",
                            "steps",
                            "lhs_matrix",
                            "field_mass",
                            "serving_mass_g",
                            "alpha",
                            "beta",
                            "unit"
                        ),
                        envir = environment()
                    )
                }
                rvs$model_results <- scrm_clustered(
                    n_iterations,
                    1,
                    n_cores,
                    steps,
                    lhs_matrix,
                    field_mass,
                    serving_mass_g,
                    alpha,
                    beta,
                    cl,
                    parallel_mode,
                    unit
                )
                
                rvs$model_results$detected_by_retail <-  as.numeric(rvs$model_results$p_positive_regulatory_test)>1 & as.numeric(rvs$model_results$testing_outcome_product_sum>0)
                rvs$model_results$detected_by_fp_testing <-  as.numeric(rvs$model_results$p_positive_regulatory_test)>1 & as.numeric(rvs$model_results$testing_outcome_regulatory_sum>0)
                rvs$p_positive_regulatory <- p_positive_regulatory_summarizer(rvs$model_results)
                
                rvs$n_positives <- sum(as.numeric(rvs$model_results$p_positive_regulatory_test))
                rvs$risk_of_positive_bins <- risk_binner(rvs$model_results, rvs$n_iterations)
                rvs$risk_of_positive_bins_new <- risk_binner_02(rvs$model_results)
                
                if (1/rvs$p_positive_regulatory$mean > 100) {
                    rvs$overall <- as.integer(round(1/rvs$p_positive_regulatory$mean, digits = -2))
                } else {
                    rvs$overall <- as.integer(round(1/rvs$p_positive_regulatory$mean, digits = 0))
                }
                
                summary_table_counts <- data.frame(
                    scenario = c(
                        as.character(rvs$scenario_name)
                    ),
                    overall = c(
                        rvs$overall
                    ),
                    x_bt_1_10 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$counts[2]))
                    )
                )
                colnames(summary_table_counts) <- c("Scenario Name", "Overall Risk of Positive Test (1 in ...)", "Count of Lots with Risk of >1 in 10 Positives")
                output$summary_table_counts <- renderTable(summary_table_counts, rownames = FALSE)
                
                scenario_name <- rvs$scenario_name
                overall <- rvs$overall
                highest_risk_count <- as.integer(rvs$risk_of_positive_bins_new$counts[2])
                
                output$summary_text <- renderUI(paste0("The overall risk under the ", scenario_name, " scenario was approximately 1 in ", overall, ", and resulted in ", highest_risk_count, " highest-risk lots."))
                
                lot_table_counts <- data.frame(
                    overall = c(
                        rvs$overall
                    ),
                    x_bt_1_10 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$counts[2]))
                    ),
                    x_bt_10_100 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$counts[3]))
                    ),
                    x_bt_100_1000 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$counts[4]))
                    ),
                    x_bt_1000_10000 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$counts[5]))
                    ),
                    x_gt_10000 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$counts[6]))
                    ),
                    x_et_0 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$counts[7]))
                    ),
                    total = c(
                        sum(
                            as.integer(rvs$risk_of_positive_bins_new$counts[2]),
                            as.integer(rvs$risk_of_positive_bins_new$counts[3]),
                            as.integer(rvs$risk_of_positive_bins_new$counts[4]),
                            as.integer(rvs$risk_of_positive_bins_new$counts[5]),
                            as.integer(rvs$risk_of_positive_bins_new$counts[6]),
                            as.integer(rvs$risk_of_positive_bins_new$counts[7])
                        )
                    )
                )
                colnames(lot_table_counts) <- c("Overall Risk of Positive Test (1 in ...)", ">1 in 10", ">1 in 100", ">1 in 1,000", ">1 in 10,000", "< 1 in 10,000", "Non-contaminated/Rejected", "Total")
                output$lot_table_counts <- renderTable(lot_table_counts, rownames = FALSE)
                rvs$lot_table_counts <- lot_table_counts
                
                lot_table_regulatory_tests <- data.frame(
                    x_bt_1_10 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[2]))
                    ),
                    x_bt_10_100 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[3]))
                    ),
                    x_bt_100_1000 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[4]))
                    ),
                    x_bt_1000_10000 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[5]))
                    ),
                    x_gt_10000 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[6]))
                    ),
                    x_et_0 = c(
                        (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[7]))
                    ),
                    total = c(
                        sum(
                            as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[2]),
                            as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[3]),
                            as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[4]),
                            as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[5]),
                            as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[6])
                        )
                    )
                )
                colnames(lot_table_regulatory_tests) <- c(">1 in 10", ">1 in 100", ">1 in 1,000", ">1 in 10,000", "< 1 in 10,000", "Non-contaminated/Rejected", "Total Positive")
                output$lot_table_regulatory_tests <- renderTable(lot_table_regulatory_tests, rownames = FALSE)
                rvs$lot_table_regulatory_tests <- lot_table_regulatory_tests
                
                time_finished <- as.character(format(Sys.time(), "%X"))
                rvs$timestamp <- as.character(paste0(format(Sys.Date()),"-", format(Sys.time(), "%H-%M-%S")))
                removeNotification(id = "ntfcn_running")
                showNotification(id = "ntfcn_finished", paste0("Finished at ", time_finished, " (Started at ", time_started, ")"), duration = NULL, closeButton = TRUE, type = "message")
                if (rvs$beep == TRUE) {
                    beep()   
                }
                # write.csv(rvs$summary_table_counts, "summary_table_counts.csv", row.names = FALSE)
            }   
        }
    })
    
    observeEvent(input$btn_save_to_slot, {
        if (length(rvs$overall)>0) {
            rvs$results_lot_table_counts[rvs$save_slot,1] <- rvs$scenario_name
            rvs$results_lot_table_counts[rvs$save_slot,2] <- rvs$overall
            rvs$results_lot_table_counts[rvs$save_slot,3] <- (as.integer(rvs$risk_of_positive_bins_new$counts[2]))
            rvs$results_lot_table_counts[rvs$save_slot,4] <- (as.integer(rvs$risk_of_positive_bins_new$counts[3]))
            rvs$results_lot_table_counts[rvs$save_slot,5] <- (as.integer(rvs$risk_of_positive_bins_new$counts[4]))
            rvs$results_lot_table_counts[rvs$save_slot,6] <- (as.integer(rvs$risk_of_positive_bins_new$counts[5]))
            rvs$results_lot_table_counts[rvs$save_slot,7] <- (as.integer(rvs$risk_of_positive_bins_new$counts[6]))
            rvs$results_lot_table_counts[rvs$save_slot,8] <- (as.integer(rvs$risk_of_positive_bins_new$counts[7]))
            rvs$results_lot_table_counts[rvs$save_slot,9] <- sum(
                as.integer(rvs$risk_of_positive_bins_new$counts[2]),
                as.integer(rvs$risk_of_positive_bins_new$counts[3]),
                as.integer(rvs$risk_of_positive_bins_new$counts[4]),
                as.integer(rvs$risk_of_positive_bins_new$counts[5]),
                as.integer(rvs$risk_of_positive_bins_new$counts[6]),
                as.integer(rvs$risk_of_positive_bins_new$counts[7])
            )
            
            output$results_lot_table_counts <- DT::renderDataTable(rvs$results_lot_table_counts, rownames = FALSE)
            
            rvs$results_lot_table_regulatory_tests[rvs$save_slot,1] <- rvs$scenario_name
            rvs$results_lot_table_regulatory_tests[rvs$save_slot,2] <- (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[2]))
            rvs$results_lot_table_regulatory_tests[rvs$save_slot,3] <- (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[3]))
            rvs$results_lot_table_regulatory_tests[rvs$save_slot,4] <- (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[4]))
            rvs$results_lot_table_regulatory_tests[rvs$save_slot,5] <- (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[5]))
            rvs$results_lot_table_regulatory_tests[rvs$save_slot,6] <- (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[6]))
            rvs$results_lot_table_regulatory_tests[rvs$save_slot,7] <- (as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[7]))
            rvs$results_lot_table_regulatory_tests[rvs$save_slot,8] <- sum(
                as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[2]),
                as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[3]),
                as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[4]),
                as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[5]),
                as.integer(rvs$risk_of_positive_bins_new$regulatory_tests[6])
            )
            
            output$results_lot_table_regulatory_tests <- DT::renderDataTable(rvs$results_lot_table_regulatory_tests, rownames = FALSE)   
        }
    })
    
    observeEvent(input$btn_load_scenario, {
        rvs$scenario <- as.character(input$slctin_scenario)
        if (rvs$scenario != "NA") {
            scenario_index <- which(list_scenario_names == rvs$scenario)
            target_scenario <- as.data.frame(list_scenarios_data[scenario_index])
            
            reset_steps()
            
            updateTextInput(session, "txtin_scenario_name", value = rvs$scenario)
            rvs$scenario_name <- input$txtin_scenario_name
            # Card 01
            if (nrow(target_scenario) >= 1) {
                updateCheckboxInput(session, "ckbxin_crd01", value = TRUE)
                updateTextInput(session, "txtin_crd01_model_step", value = target_scenario$model_step[1])
                updateSelectInput(session, "slctin_crd01_microbial_process", selected = target_scenario$microbial_process[1])
                updateNumericInput(session, "nmcin_crd01_p_occurrence", value = target_scenario$p_occurrence[1])
                updateSelectInput(session, "slctin_crd01_distribution", selected = target_scenario$distribution[1])
                updateNumericInput(session, "slctin_crd01_contamination_basis", value = target_scenario$contamination_basis[1])
                updateNumericInput(session, "nmcin_crd01_parameter_1", value = target_scenario$parameter_1[1])
                updateNumericInput(session, "nmcin_crd01_parameter_2", value = target_scenario$parameter_2[1])
                updateNumericInput(session, "nmcin_crd01_parameter_3", value = target_scenario$parameter_3[1])
                rvs$crd01 <- TRUE
                rvs$crd01_model_step <- target_scenario$model_step[1]
                rvs$crd01_microbial_process <- target_scenario$microbial_process[1]
                rvs$crd01_p_occurrence <- target_scenario$p_occurrence[1]
                rvs$crd01_distribution <- target_scenario$distribution[1]
                rvs$crd01_contamination_basis <- target_scenario$contamination_basis[1]
                rvs$crd01_parameter_1 <- target_scenario$parameter_1[1]
                rvs$crd01_parameter_2 <- target_scenario$parameter_2[1]
                rvs$crd01_parameter_3 <- target_scenario$parameter_3[1]
            }
            # Card 02
            if (nrow(target_scenario) >= 2) {
                updateCheckboxInput(session, "ckbxin_crd02", value = TRUE)
                updateTextInput(session, "txtin_crd02_model_step", value = target_scenario$model_step[2])
                updateSelectInput(session, "slctin_crd02_microbial_process", selected = target_scenario$microbial_process[2])
                updateNumericInput(session, "nmcin_crd02_p_occurrence", value = target_scenario$p_occurrence[2])
                updateSelectInput(session, "slctin_crd02_distribution", selected = target_scenario$distribution[2])
                updateNumericInput(session, "slctin_crd02_contamination_basis", value = target_scenario$contamination_basis[2])
                updateNumericInput(session, "nmcin_crd02_parameter_1", value = target_scenario$parameter_1[2])
                updateNumericInput(session, "nmcin_crd02_parameter_2", value = target_scenario$parameter_2[2])
                updateNumericInput(session, "nmcin_crd02_parameter_3", value = target_scenario$parameter_3[2])
                rvs$crd02 <- TRUE
                rvs$crd02_model_step <- target_scenario$model_step[2]
                rvs$crd02_microbial_process <- target_scenario$microbial_process[2]
                rvs$crd02_p_occurrence <- target_scenario$p_occurrence[2]
                rvs$crd02_distribution <- target_scenario$distribution[2]
                rvs$crd02_contamination_basis <- target_scenario$contamination_basis[2]
                rvs$crd02_parameter_1 <- target_scenario$parameter_1[2]
                rvs$crd02_parameter_2 <- target_scenario$parameter_2[2]
                rvs$crd02_parameter_3 <- target_scenario$parameter_3[2]
            }
            # Card 03
            if (nrow(target_scenario) >= 3) {
                updateCheckboxInput(session, "ckbxin_crd03", value = TRUE)
                updateTextInput(session, "txtin_crd03_model_step", value = target_scenario$model_step[3])
                updateSelectInput(session, "slctin_crd03_microbial_process", selected = target_scenario$microbial_process[3])
                updateNumericInput(session, "nmcin_crd03_p_occurrence", value = target_scenario$p_occurrence[3])
                updateSelectInput(session, "slctin_crd03_distribution", selected = target_scenario$distribution[3])
                updateNumericInput(session, "slctin_crd03_contamination_basis", value = target_scenario$contamination_basis[3])
                updateNumericInput(session, "nmcin_crd03_parameter_1", value = target_scenario$parameter_1[3])
                updateNumericInput(session, "nmcin_crd03_parameter_2", value = target_scenario$parameter_2[3])
                updateNumericInput(session, "nmcin_crd03_parameter_3", value = target_scenario$parameter_3[3])
                rvs$crd03 <- TRUE
                rvs$crd03_model_step <- target_scenario$model_step[3]
                rvs$crd03_microbial_process <- target_scenario$microbial_process[3]
                rvs$crd03_p_occurrence <- target_scenario$p_occurrence[3]
                rvs$crd03_distribution <- target_scenario$distribution[3]
                rvs$crd03_contamination_basis <- target_scenario$contamination_basis[3]
                rvs$crd03_parameter_1 <- target_scenario$parameter_1[3]
                rvs$crd03_parameter_2 <- target_scenario$parameter_2[3]
                rvs$crd03_parameter_3 <- target_scenario$parameter_3[3]
            }
            # Card 04
            if (nrow(target_scenario) >= 4) {
                updateCheckboxInput(session, "ckbxin_crd04", value = TRUE)
                updateTextInput(session, "txtin_crd04_model_step", value = target_scenario$model_step[4])
                updateSelectInput(session, "slctin_crd04_microbial_process", selected = target_scenario$microbial_process[4])
                updateNumericInput(session, "nmcin_crd04_p_occurrence", value = target_scenario$p_occurrence[4])
                updateSelectInput(session, "slctin_crd04_distribution", selected = target_scenario$distribution[4])
                updateNumericInput(session, "slctin_crd04_contamination_basis", value = target_scenario$contamination_basis[4])
                updateNumericInput(session, "nmcin_crd04_parameter_1", value = target_scenario$parameter_1[4])
                updateNumericInput(session, "nmcin_crd04_parameter_2", value = target_scenario$parameter_2[4])
                updateNumericInput(session, "nmcin_crd04_parameter_3", value = target_scenario$parameter_3[4])
                rvs$crd04 <- TRUE
                rvs$crd04_model_step <- target_scenario$model_step[4]
                rvs$crd04_microbial_process <- target_scenario$microbial_process[4]
                rvs$crd04_p_occurrence <- target_scenario$p_occurrence[4]
                rvs$crd04_distribution <- target_scenario$distribution[4]
                rvs$crd04_contamination_basis <- target_scenario$contamination_basis[4]
                rvs$crd04_parameter_1 <- target_scenario$parameter_1[4]
                rvs$crd04_parameter_2 <- target_scenario$parameter_2[4]
                rvs$crd04_parameter_3 <- target_scenario$parameter_3[4]
            }
            # Card 05
            if (nrow(target_scenario) >= 5) {
                updateCheckboxInput(session, "ckbxin_crd05", value = TRUE)
                updateTextInput(session, "txtin_crd05_model_step", value = target_scenario$model_step[5])
                updateSelectInput(session, "slctin_crd05_microbial_process", selected = target_scenario$microbial_process[5])
                updateNumericInput(session, "nmcin_crd05_p_occurrence", value = target_scenario$p_occurrence[5])
                updateSelectInput(session, "slctin_crd05_distribution", selected = target_scenario$distribution[5])
                updateNumericInput(session, "slctin_crd05_contamination_basis", value = target_scenario$contamination_basis[5])
                updateNumericInput(session, "nmcin_crd05_parameter_1", value = target_scenario$parameter_1[5])
                updateNumericInput(session, "nmcin_crd05_parameter_2", value = target_scenario$parameter_2[5])
                updateNumericInput(session, "nmcin_crd05_parameter_3", value = target_scenario$parameter_3[5])
                rvs$crd05 <- TRUE
                rvs$crd05_model_step <- target_scenario$model_step[5]
                rvs$crd05_microbial_process <- target_scenario$microbial_process[5]
                rvs$crd05_p_occurrence <- target_scenario$p_occurrence[5]
                rvs$crd05_distribution <- target_scenario$distribution[5]
                rvs$crd05_contamination_basis <- target_scenario$contamination_basis[5]
                rvs$crd05_parameter_1 <- target_scenario$parameter_1[5]
                rvs$crd05_parameter_2 <- target_scenario$parameter_2[5]
                rvs$crd05_parameter_3 <- target_scenario$parameter_3[5]
            }
            # Card 06
            if (nrow(target_scenario) >= 6) {
                updateCheckboxInput(session, "ckbxin_crd06", value = TRUE)
                updateTextInput(session, "txtin_crd06_model_step", value = target_scenario$model_step[6])
                updateSelectInput(session, "slctin_crd06_microbial_process", selected = target_scenario$microbial_process[6])
                updateNumericInput(session, "nmcin_crd06_p_occurrence", value = target_scenario$p_occurrence[6])
                updateSelectInput(session, "slctin_crd06_distribution", selected = target_scenario$distribution[6])
                updateNumericInput(session, "slctin_crd06_contamination_basis", value = target_scenario$contamination_basis[6])
                updateNumericInput(session, "nmcin_crd06_parameter_1", value = target_scenario$parameter_1[6])
                updateNumericInput(session, "nmcin_crd06_parameter_2", value = target_scenario$parameter_2[6])
                updateNumericInput(session, "nmcin_crd06_parameter_3", value = target_scenario$parameter_3[6])
                rvs$crd06 <- TRUE
                rvs$crd06_model_step <- target_scenario$model_step[6]
                rvs$crd06_microbial_process <- target_scenario$microbial_process[6]
                rvs$crd06_p_occurrence <- target_scenario$p_occurrence[6]
                rvs$crd06_distribution <- target_scenario$distribution[6]
                rvs$crd06_contamination_basis <- target_scenario$contamination_basis[6]
                rvs$crd06_parameter_1 <- target_scenario$parameter_1[6]
                rvs$crd06_parameter_2 <- target_scenario$parameter_2[6]
                rvs$crd06_parameter_3 <- target_scenario$parameter_3[6]
            }
            # Card 07
            if (nrow(target_scenario) >= 7) {
                updateCheckboxInput(session, "ckbxin_crd07", value = TRUE)
                updateTextInput(session, "txtin_crd07_model_step", value = target_scenario$model_step[7])
                updateSelectInput(session, "slctin_crd07_microbial_process", selected = target_scenario$microbial_process[7])
                updateNumericInput(session, "nmcin_crd07_p_occurrence", value = target_scenario$p_occurrence[7])
                updateSelectInput(session, "slctin_crd07_distribution", selected = target_scenario$distribution[7])
                updateNumericInput(session, "slctin_crd07_contamination_basis", value = target_scenario$contamination_basis[7])
                updateNumericInput(session, "nmcin_crd07_parameter_1", value = target_scenario$parameter_1[7])
                updateNumericInput(session, "nmcin_crd07_parameter_2", value = target_scenario$parameter_2[7])
                updateNumericInput(session, "nmcin_crd07_parameter_3", value = target_scenario$parameter_3[7])
                rvs$crd07 <- TRUE
                rvs$crd07_model_step <- target_scenario$model_step[7]
                rvs$crd07_microbial_process <- target_scenario$microbial_process[7]
                rvs$crd07_p_occurrence <- target_scenario$p_occurrence[7]
                rvs$crd07_distribution <- target_scenario$distribution[7]
                rvs$crd07_contamination_basis <- target_scenario$contamination_basis[7]
                rvs$crd07_parameter_1 <- target_scenario$parameter_1[7]
                rvs$crd07_parameter_2 <- target_scenario$parameter_2[7]
                rvs$crd07_parameter_3 <- target_scenario$parameter_3[7]
            }
            # Card 08
            if (nrow(target_scenario) >= 8) {
                updateCheckboxInput(session, "ckbxin_crd08", value = TRUE)
                updateTextInput(session, "txtin_crd08_model_step", value = target_scenario$model_step[8])
                updateSelectInput(session, "slctin_crd08_microbial_process", selected = target_scenario$microbial_process[8])
                updateNumericInput(session, "nmcin_crd08_p_occurrence", value = target_scenario$p_occurrence[8])
                updateSelectInput(session, "slctin_crd08_distribution", selected = target_scenario$distribution[8])
                updateNumericInput(session, "slctin_crd08_contamination_basis", value = target_scenario$contamination_basis[8])
                updateNumericInput(session, "nmcin_crd08_parameter_1", value = target_scenario$parameter_1[8])
                updateNumericInput(session, "nmcin_crd08_parameter_2", value = target_scenario$parameter_2[8])
                updateNumericInput(session, "nmcin_crd08_parameter_3", value = target_scenario$parameter_3[8])
                rvs$crd08 <- TRUE
                rvs$crd08_model_step <- target_scenario$model_step[8]
                rvs$crd08_microbial_process <- target_scenario$microbial_process[8]
                rvs$crd08_p_occurrence <- target_scenario$p_occurrence[8]
                rvs$crd08_distribution <- target_scenario$distribution[8]
                rvs$crd08_contamination_basis <- target_scenario$contamination_basis[8]
                rvs$crd08_parameter_1 <- target_scenario$parameter_1[8]
                rvs$crd08_parameter_2 <- target_scenario$parameter_2[8]
                rvs$crd08_parameter_3 <- target_scenario$parameter_3[8]
            }
            # Card 09
            if (nrow(target_scenario) == 9) {
                updateCheckboxInput(session, "ckbxin_crd09", value = TRUE)
                updateTextInput(session, "txtin_crd09_model_step", value = target_scenario$model_step[9])
                updateSelectInput(session, "slctin_crd09_microbial_process", selected = target_scenario$microbial_process[9])
                updateNumericInput(session, "nmcin_crd09_p_occurrence", value = target_scenario$p_occurrence[9])
                updateSelectInput(session, "slctin_crd09_distribution", selected = target_scenario$distribution[9])
                updateNumericInput(session, "slctin_crd09_contamination_basis", value = target_scenario$contamination_basis[9])
                updateNumericInput(session, "nmcin_crd09_parameter_1", value = target_scenario$parameter_1[9])
                updateNumericInput(session, "nmcin_crd09_parameter_2", value = target_scenario$parameter_2[9])
                updateNumericInput(session, "nmcin_crd09_parameter_3", value = target_scenario$parameter_3[9])
                rvs$crd09 <- TRUE
                rvs$crd09_model_step <- target_scenario$model_step[9]
                rvs$crd09_microbial_process <- target_scenario$microbial_process[9]
                rvs$crd09_p_occurrence <- target_scenario$p_occurrence[9]
                rvs$crd09_distribution <- target_scenario$distribution[9]
                rvs$crd09_contamination_basis <- target_scenario$contamination_basis[9]
                rvs$crd09_parameter_1 <- target_scenario$parameter_1[9]
                rvs$crd09_parameter_2 <- target_scenario$parameter_2[9]
                rvs$crd09_parameter_3 <- target_scenario$parameter_3[9]
            }
        }
    })
    
    observeEvent(toListenTextinScenarioName(), {
        rvs$scenario_name <- input$txtin_scenario_name
    })
    
    observeEvent(toListenNumericinFieldMassData(), {
        rvs$field_mass_unit <- input$slctin_field_mass_unit
        if (rvs$field_mass_unit == "lb") {
            rvs$field_mass <- input$nmcin_field_mass*454
        } else {
            rvs$field_mass <- input$nmcin_field_mass
        }
    })
    
    observeEvent(toListenNumericinServingMass(), {
        # rvs$serving_mass <- input$nmcin_serving_mass
    })
    
    observeEvent(toListenNumericinAlpha(), {
        # rvs$alpha <- input$nmcin_alpha
    })
    
    observeEvent(toListenNumericinBeta(), {
        # rvs$beta <- input$nmcin_beta
    })
    
    observeEvent(toListenSelectinSaveSlot(), {
        rvs$save_slot <- input$slctin_save_slot
    })
    
    observeEvent(toListenCheckboxBeep(), {
        rvs$beep <- input$ckbxin_beep
    })
    
    observeEvent(toListenNumericinIterations(), {
        rvs$n_iterations <- input$nmcin_n_iterations
    })
    
    observeEvent(toListenCheckboxCard01(), {
        rvs$crd01 <- input$ckbxin_crd01
    })
    observeEvent(toListenTextinCard01ModelStep(), {
        rvs$crd01_model_step <- input$txtin_crd01_model_step
    })
    observeEvent(toListenSelectinCard01MicrobialProcess(), {
        rvs$crd01_microbial_process <- input$slctin_crd01_microbial_process
        if (rvs$crd01_microbial_process == "Contamination/Removal") {
            update_tooltip("tltp_crd01_microbial_process", "Contamination Maths:\nCarried Contamination (CFU/g) = Carried Contamination (CFU/g) + 10^Modeled Value (log10 CFU/g)")
            enable("slctin_crd01_contamination_basis")
            show("slctin_crd01_contamination_basis")
        } else if (rvs$crd01_microbial_process == "Increase/Reduction") {
            update_tooltip("tltp_crd01_microbial_process", "Increase/Reduction Maths:\nCarried Contamination (CFU/g) = 10^(log10(Carried Contamination (CFU/g)) + Modeled Value (log10 CFU/g))")
            disable("slctin_crd01_contamination_basis")
            hide("slctin_crd01_contamination_basis")
        } else {
            update_tooltip("tltp_crd01_microbial_process", "")
            disable("slctin_crd01_contamination_basis")
            hide("slctin_crd01_contamination_basis")
        }
        if (rvs$crd01_microbial_process == "Risk Output Test" | rvs$crd01_microbial_process == "Product Test") {
            disable("slctin_crd01_distribution")
            hide("slctin_crd01_distribution")
            updateNumericInput(session, "nmcin_crd01_parameter_1", label = "Composite Mass (g)")
            updateNumericInput(session, "nmcin_crd01_parameter_2", label = "Grabs (N)")
            updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Tests (N)")
        } else if (rvs$crd01_microbial_process == "Contamination/Removal") {
            enable("slctin_crd01_distribution")
            show("slctin_crd01_distribution")
            updateSelectInput(session, "slctin_crd01_distribution", choices = list_contamination_distributions, selected = list_contamination_distributions[1])
        } else if (rvs$crd01_microbial_process == "Increase/Reduction") {
            enable("slctin_crd01_distribution")
            show("slctin_crd01_distribution")
            updateSelectInput(session, "slctin_crd01_distribution", choices = list_distributions, selected = list_distributions[1])
        }
    })
    observeEvent(toListenSelectinCard01ProcessAndDist(), {
        rvs$crd01_microbial_process <- input$slctin_crd01_microbial_process
        if (rvs$crd01_microbial_process == "Contamination/Removal") {
            rvs$crd01_distribution <- input$slctin_crd01_distribution
            if (rvs$crd01_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd01_parameter_3")
            } else if (rvs$crd01_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd01_parameter_3")
            } else if (rvs$crd01_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd01_parameter_3")
            } else if (rvs$crd01_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd01_parameter_3")
            } else if (rvs$crd01_distribution == "Normal (linear)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd01_parameter_3")
            } else if (rvs$crd01_distribution == "Uniform (linear)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd01_parameter_3")
            } else if (rvs$crd01_distribution == "Pert (linear)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd01_parameter_3")
            } else if (rvs$crd01_distribution == "Triangular (linear)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd01_parameter_3")
            }
        } else if (rvs$crd01_microbial_process == "Increase/Reduction") {
            rvs$crd01_distribution <- input$slctin_crd01_distribution
            if (rvs$crd01_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = "Mean (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = "SD (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd01_parameter_3")
            } else if (rvs$crd01_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = "Lower Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = "Upper Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd01_parameter_3")
            } else if (rvs$crd01_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd01_parameter_3")
            } else if (rvs$crd01_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd01_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd01_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd01_parameter_3")
            }
        } else {
            shinyjs::show(id = "nmcin_crd01_parameter_3")
        }
    })
    observeEvent(toListenNumberinCard01POccurrence(), {
        rvs$crd01_p_occurrence <- input$nmcin_crd01_p_occurrence
    })
    observeEvent(toListenSelectinCard01ContaminationBasis(), {
        rvs$crd01_contamination_basis <- input$slctin_crd01_contamination_basis
        if (rvs$crd01_contamination_basis == list_contamination_types[1]) {
            rvs$contamination_unit <- ""
        } else if (rvs$crd01_contamination_basis == list_contamination_types[2]) {
            rvs$contamination_unit <- "/g"
        } else if (rvs$crd01_contamination_basis == list_contamination_types[3]) {
            rvs$contamination_unit <- "/lb"
        }
        if (rvs$crd01_distribution == "Normal (log10)") {
            updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Not Used")
        } else if (rvs$crd01_distribution == "Uniform (log10)") {
            updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Not Used")
        } else if (rvs$crd01_distribution == "Pert (log10)") {
            updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
            shinyjs::show(id = "nmcin_crd01_parameter_3")
        } else if (rvs$crd01_distribution == "Triangular (log10)") {
            updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd01_distribution == "Normal (linear)") {
            updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Not Used")
        } else if (rvs$crd01_distribution == "Uniform (linear)") {
            updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_3", label = "Not Used")
        } else if (rvs$crd01_distribution == "Pert (linear)") {
            updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd01_distribution == "Triangular (linear)") {
            updateNumericInput(session, "nmcin_crd01_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd01_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        }
    })
    observeEvent(toListenNumberinCard01Parameter1(), {
        rvs$crd01_parameter_1 <- input$nmcin_crd01_parameter_1 
    })
    observeEvent(toListenNumberinCard01Parameter2(), {
        rvs$crd01_parameter_2 <- input$nmcin_crd01_parameter_2 
    })
    observeEvent(toListenNumberinCard01Parameter3(), {
        rvs$crd01_parameter_3 <- input$nmcin_crd01_parameter_3 
    })
    
    observeEvent(toListenCheckboxCard02(), {
        rvs$crd02 <- input$ckbxin_crd02
    })
    observeEvent(toListenTextinCard02ModelStep(), {
        rvs$crd02_model_step <- input$txtin_crd02_model_step
    })
    observeEvent(toListenSelectinCard02MicrobialProcess(), {
        rvs$crd02_microbial_process <- input$slctin_crd02_microbial_process
        if (rvs$crd02_microbial_process == "Contamination/Removal") {
            update_tooltip("tltp_crd02_microbial_process", "Contamination Maths:\nCarried Contamination (CFU/g) = Carried Contamination (CFU/g) + 10^Modeled Value (log10 CFU/g)")
            enable("slctin_crd02_contamination_basis")
            show("slctin_crd02_contamination_basis")
        } else if (rvs$crd02_microbial_process == "Increase/Reduction") {
            update_tooltip("tltp_crd02_microbial_process", "Increase/Reduction Maths:\nCarried Contamination (CFU/g) = 10^(log10(Carried Contamination (CFU/g)) + Modeled Value (log10 CFU/g))")
            disable("slctin_crd02_contamination_basis")
            hide("slctin_crd02_contamination_basis")
        } else {
            update_tooltip("tltp_crd02_microbial_process", "")
            disable("slctin_crd02_contamination_basis")
            hide("slctin_crd02_contamination_basis")
        }
        if (rvs$crd02_microbial_process == "Risk Output Test" | rvs$crd02_microbial_process == "Product Test") {
            disable("slctin_crd02_distribution")
            hide("slctin_crd02_distribution")
            updateNumericInput(session, "nmcin_crd02_parameter_1", label = "Composite Mass (g)")
            updateNumericInput(session, "nmcin_crd02_parameter_2", label = "Grabs (N)")
            updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Tests (N)")
        } else if (rvs$crd02_microbial_process == "Contamination/Removal") {
            enable("slctin_crd02_distribution")
            show("slctin_crd02_distribution")
            updateSelectInput(session, "slctin_crd02_distribution", choices = list_contamination_distributions, selected = list_contamination_distributions[1])
        } else if (rvs$crd02_microbial_process == "Increase/Reduction") {
            enable("slctin_crd02_distribution")
            show("slctin_crd02_distribution")
            updateSelectInput(session, "slctin_crd02_distribution", choices = list_distributions, selected = list_distributions[1])
        }
    })
    observeEvent(toListenSelectinCard02ProcessAndDist(), {
        rvs$crd02_microbial_process <- input$slctin_crd02_microbial_process
        if (rvs$crd02_microbial_process == "Contamination/Removal") {
            rvs$crd02_distribution <- input$slctin_crd02_distribution
            if (rvs$crd02_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd02_parameter_3")
            } else if (rvs$crd02_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd02_parameter_3")
            } else if (rvs$crd02_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd02_parameter_3")
            } else if (rvs$crd02_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd02_parameter_3")
            } else if (rvs$crd02_distribution == "Normal (linear)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd02_parameter_3")
            } else if (rvs$crd02_distribution == "Uniform (linear)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd02_parameter_3")
            } else if (rvs$crd02_distribution == "Pert (linear)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd02_parameter_3")
            } else if (rvs$crd02_distribution == "Triangular (linear)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd02_parameter_3")
            }
        } else if (rvs$crd02_microbial_process == "Increase/Reduction") {
            rvs$crd02_distribution <- input$slctin_crd02_distribution
            if (rvs$crd02_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = "Mean (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = "SD (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd02_parameter_3")
            } else if (rvs$crd02_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = "Lower Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = "Upper Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd02_parameter_3")
            } else if (rvs$crd02_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd02_parameter_3")
            } else if (rvs$crd02_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd02_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd02_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd02_parameter_3")
            }
        } else {
            shinyjs::show(id = "nmcin_crd02_parameter_3")
        }
    })
    observeEvent(toListenNumberinCard02POccurrence(), {
        rvs$crd02_p_occurrence <- input$nmcin_crd02_p_occurrence
    })
    observeEvent(toListenSelectinCard02ContaminationBasis(), {
        rvs$crd02_contamination_basis <- input$slctin_crd02_contamination_basis
        if (rvs$crd02_contamination_basis == list_contamination_types[1]) {
            rvs$contamination_unit <- ""
        } else if (rvs$crd02_contamination_basis == list_contamination_types[2]) {
            rvs$contamination_unit <- "/g"
        } else if (rvs$crd02_contamination_basis == list_contamination_types[3]) {
            rvs$contamination_unit <- "/lb"
        }
        if (rvs$crd02_distribution == "Normal (log10)") {
            updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Not Used")
        } else if (rvs$crd02_distribution == "Uniform (log10)") {
            updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Not Used")
        } else if (rvs$crd02_distribution == "Pert (log10)") {
            updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
            shinyjs::show(id = "nmcin_crd02_parameter_3")
        } else if (rvs$crd02_distribution == "Triangular (log10)") {
            updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd02_distribution == "Normal (linear)") {
            updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Not Used")
        } else if (rvs$crd02_distribution == "Uniform (linear)") {
            updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_3", label = "Not Used")
        } else if (rvs$crd02_distribution == "Pert (linear)") {
            updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd02_distribution == "Triangular (linear)") {
            updateNumericInput(session, "nmcin_crd02_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd02_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        }
    })
    observeEvent(toListenNumberinCard02Parameter1(), {
        rvs$crd02_parameter_1 <- input$nmcin_crd02_parameter_1 
    })
    observeEvent(toListenNumberinCard02Parameter2(), {
        rvs$crd02_parameter_2 <- input$nmcin_crd02_parameter_2 
    })
    observeEvent(toListenNumberinCard02Parameter3(), {
        rvs$crd02_parameter_3 <- input$nmcin_crd02_parameter_3
    })
    
    observeEvent(toListenCheckboxCard03(), {
        rvs$crd03 <- input$ckbxin_crd03
    })
    observeEvent(toListenTextinCard03ModelStep(), {
        rvs$crd03_model_step <- input$txtin_crd03_model_step
    })
    observeEvent(toListenSelectinCard03MicrobialProcess(), {
        rvs$crd03_microbial_process <- input$slctin_crd03_microbial_process
        if (rvs$crd03_microbial_process == "Contamination/Removal") {
            update_tooltip("tltp_crd03_microbial_process", "Contamination Maths:\nCarried Contamination (CFU/g) = Carried Contamination (CFU/g) + 10^Modeled Value (log10 CFU/g)")
            enable("slctin_crd03_contamination_basis")
            show("slctin_crd03_contamination_basis")
        } else if (rvs$crd03_microbial_process == "Increase/Reduction") {
            update_tooltip("tltp_crd03_microbial_process", "Increase/Reduction Maths:\nCarried Contamination (CFU/g) = 10^(log10(Carried Contamination (CFU/g)) + Modeled Value (log10 CFU/g))")
            disable("slctin_crd03_contamination_basis")
            hide("slctin_crd03_contamination_basis")
        } else {
            update_tooltip("tltp_crd03_microbial_process", "")
            disable("slctin_crd03_contamination_basis")
            hide("slctin_crd03_contamination_basis")
        }
        if (rvs$crd03_microbial_process == "Risk Output Test" | rvs$crd03_microbial_process == "Product Test") {
            disable("slctin_crd03_distribution")
            hide("slctin_crd03_distribution")
            updateNumericInput(session, "nmcin_crd03_parameter_1", label = "Composite Mass (g)")
            updateNumericInput(session, "nmcin_crd03_parameter_2", label = "Grabs (N)")
            updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Tests (N)")
        } else if (rvs$crd03_microbial_process == "Contamination/Removal") {
            enable("slctin_crd03_distribution")
            show("slctin_crd03_distribution")
            updateSelectInput(session, "slctin_crd03_distribution", choices = list_contamination_distributions, selected = list_contamination_distributions[1])
        } else if (rvs$crd03_microbial_process == "Increase/Reduction") {
            enable("slctin_crd03_distribution")
            show("slctin_crd03_distribution")
            updateSelectInput(session, "slctin_crd03_distribution", choices = list_distributions, selected = list_distributions[1])
        }
    })
    observeEvent(toListenSelectinCard03ProcessAndDist(), {
        rvs$crd03_microbial_process <- input$slctin_crd03_microbial_process
        if (rvs$crd03_microbial_process == "Contamination/Removal") {
            rvs$crd03_distribution <- input$slctin_crd03_distribution
            if (rvs$crd03_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd03_parameter_3")
            } else if (rvs$crd03_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd03_parameter_3")
            } else if (rvs$crd03_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd03_parameter_3")
            } else if (rvs$crd03_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd03_parameter_3")
            } else if (rvs$crd03_distribution == "Normal (linear)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd03_parameter_3")
            } else if (rvs$crd03_distribution == "Uniform (linear)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd03_parameter_3")
            } else if (rvs$crd03_distribution == "Pert (linear)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd03_parameter_3")
            } else if (rvs$crd03_distribution == "Triangular (linear)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd03_parameter_3")
            }
        } else if (rvs$crd03_microbial_process == "Increase/Reduction") {
            rvs$crd03_distribution <- input$slctin_crd03_distribution
            if (rvs$crd03_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = "Mean (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = "SD (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd03_parameter_3")
            } else if (rvs$crd03_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = "Lower Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = "Upper Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd03_parameter_3")
            } else if (rvs$crd03_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd03_parameter_3")
            } else if (rvs$crd03_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd03_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd03_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd03_parameter_3")
            }
        } else {
            shinyjs::show(id = "nmcin_crd03_parameter_3")
        }
    })
    observeEvent(toListenNumberinCard03POccurrence(), {
        rvs$crd03_p_occurrence <- input$nmcin_crd03_p_occurrence
    })
    observeEvent(toListenSelectinCard03ContaminationBasis(), {
        rvs$crd03_contamination_basis <- input$slctin_crd03_contamination_basis
        if (rvs$crd03_contamination_basis == list_contamination_types[1]) {
            rvs$contamination_unit <- ""
        } else if (rvs$crd03_contamination_basis == list_contamination_types[2]) {
            rvs$contamination_unit <- "/g"
        } else if (rvs$crd03_contamination_basis == list_contamination_types[3]) {
            rvs$contamination_unit <- "/lb"
        }
        if (rvs$crd03_distribution == "Normal (log10)") {
            updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Not Used")
        } else if (rvs$crd03_distribution == "Uniform (log10)") {
            updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Not Used")
        } else if (rvs$crd03_distribution == "Pert (log10)") {
            updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
            shinyjs::show(id = "nmcin_crd03_parameter_3")
        } else if (rvs$crd03_distribution == "Triangular (log10)") {
            updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd03_distribution == "Normal (linear)") {
            updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Not Used")
        } else if (rvs$crd03_distribution == "Uniform (linear)") {
            updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_3", label = "Not Used")
        } else if (rvs$crd03_distribution == "Pert (linear)") {
            updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd03_distribution == "Triangular (linear)") {
            updateNumericInput(session, "nmcin_crd03_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd03_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        }
    })
    observeEvent(toListenNumberinCard03Parameter1(), {
        rvs$crd03_parameter_1 <- input$nmcin_crd03_parameter_1 
    })
    observeEvent(toListenNumberinCard03Parameter2(), {
        rvs$crd03_parameter_2 <- input$nmcin_crd03_parameter_2 
    })
    observeEvent(toListenNumberinCard03Parameter3(), {
        rvs$crd03_parameter_3 <- input$nmcin_crd03_parameter_3 
    })
    
    observeEvent(toListenCheckboxCard04(), {
        rvs$crd04 <- input$ckbxin_crd04
    })
    observeEvent(toListenTextinCard04ModelStep(), {
        rvs$crd04_model_step <- input$txtin_crd04_model_step
    })
    observeEvent(toListenSelectinCard04MicrobialProcess(), {
        rvs$crd04_microbial_process <- input$slctin_crd04_microbial_process
        if (rvs$crd04_microbial_process == "Contamination/Removal") {
            update_tooltip("tltp_crd04_microbial_process", "Contamination Maths:\nCarried Contamination (CFU/g) = Carried Contamination (CFU/g) + 10^Modeled Value (log10 CFU/g)")
            enable("slctin_crd04_contamination_basis")
            show("slctin_crd04_contamination_basis")
        } else if (rvs$crd04_microbial_process == "Increase/Reduction") {
            update_tooltip("tltp_crd04_microbial_process", "Increase/Reduction Maths:\nCarried Contamination (CFU/g) = 10^(log10(Carried Contamination (CFU/g)) + Modeled Value (log10 CFU/g))")
            disable("slctin_crd04_contamination_basis")
            hide("slctin_crd04_contamination_basis")
        } else {
            update_tooltip("tltp_crd04_microbial_process", "")
            disable("slctin_crd04_contamination_basis")
            hide("slctin_crd04_contamination_basis")
        }
        if (rvs$crd04_microbial_process == "Risk Output Test" | rvs$crd04_microbial_process == "Product Test") {
            disable("slctin_crd04_distribution")
            hide("slctin_crd04_distribution")
            updateNumericInput(session, "nmcin_crd04_parameter_1", label = "Composite Mass (g)")
            updateNumericInput(session, "nmcin_crd04_parameter_2", label = "Grabs (N)")
            updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Tests (N)")
        } else if (rvs$crd04_microbial_process == "Contamination/Removal") {
            enable("slctin_crd04_distribution")
            show("slctin_crd04_distribution")
            updateSelectInput(session, "slctin_crd04_distribution", choices = list_contamination_distributions, selected = list_contamination_distributions[1])
        } else if (rvs$crd04_microbial_process == "Increase/Reduction") {
            enable("slctin_crd04_distribution")
            show("slctin_crd04_distribution")
            updateSelectInput(session, "slctin_crd04_distribution", choices = list_distributions, selected = list_distributions[1])
        }
    })
    observeEvent(toListenSelectinCard04ProcessAndDist(), {
        rvs$crd04_microbial_process <- input$slctin_crd04_microbial_process
        if (rvs$crd04_microbial_process == "Contamination/Removal") {
            rvs$crd04_distribution <- input$slctin_crd04_distribution
            if (rvs$crd04_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd04_parameter_3")
            } else if (rvs$crd04_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd04_parameter_3")
            } else if (rvs$crd04_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd04_parameter_3")
            } else if (rvs$crd04_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd04_parameter_3")
            } else if (rvs$crd04_distribution == "Normal (linear)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd04_parameter_3")
            } else if (rvs$crd04_distribution == "Uniform (linear)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd04_parameter_3")
            } else if (rvs$crd04_distribution == "Pert (linear)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd04_parameter_3")
            } else if (rvs$crd04_distribution == "Triangular (linear)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd04_parameter_3")
            }
        } else if (rvs$crd04_microbial_process == "Increase/Reduction") {
            rvs$crd04_distribution <- input$slctin_crd04_distribution
            if (rvs$crd04_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = "Mean (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = "SD (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd04_parameter_3")
            } else if (rvs$crd04_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = "Lower Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = "Upper Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd04_parameter_3")
            } else if (rvs$crd04_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd04_parameter_3")
            } else if (rvs$crd04_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd04_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd04_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd04_parameter_3")
            }
        } else {
            shinyjs::show(id = "nmcin_crd04_parameter_3")
        }
    })
    observeEvent(toListenNumberinCard04POccurrence(), {
        rvs$crd04_p_occurrence <- input$nmcin_crd04_p_occurrence
    })
    observeEvent(toListenSelectinCard04ContaminationBasis(), {
        rvs$crd04_contamination_basis <- input$slctin_crd04_contamination_basis
        if (rvs$crd04_contamination_basis == list_contamination_types[1]) {
            rvs$contamination_unit <- ""
        } else if (rvs$crd04_contamination_basis == list_contamination_types[2]) {
            rvs$contamination_unit <- "/g"
        } else if (rvs$crd04_contamination_basis == list_contamination_types[3]) {
            rvs$contamination_unit <- "/lb"
        }
        if (rvs$crd04_distribution == "Normal (log10)") {
            updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Not Used")
        } else if (rvs$crd04_distribution == "Uniform (log10)") {
            updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Not Used")
        } else if (rvs$crd04_distribution == "Pert (log10)") {
            updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
            shinyjs::show(id = "nmcin_crd04_parameter_3")
        } else if (rvs$crd04_distribution == "Triangular (log10)") {
            updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd04_distribution == "Normal (linear)") {
            updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Not Used")
        } else if (rvs$crd04_distribution == "Uniform (linear)") {
            updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_3", label = "Not Used")
        } else if (rvs$crd04_distribution == "Pert (linear)") {
            updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd04_distribution == "Triangular (linear)") {
            updateNumericInput(session, "nmcin_crd04_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd04_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        }
    })
    observeEvent(toListenNumberinCard04Parameter1(), {
        rvs$crd04_parameter_1 <- input$nmcin_crd04_parameter_1 
    })
    observeEvent(toListenNumberinCard04Parameter2(), {
        rvs$crd04_parameter_2 <- input$nmcin_crd04_parameter_2 
    })
    observeEvent(toListenNumberinCard04Parameter3(), {
        rvs$crd04_parameter_3 <- input$nmcin_crd04_parameter_3 
    })
    
    observeEvent(toListenCheckboxCard05(), {
        rvs$crd05 <- input$ckbxin_crd05
    })
    observeEvent(toListenTextinCard05ModelStep(), {
        rvs$crd05_model_step <- input$txtin_crd05_model_step
    })
    observeEvent(toListenSelectinCard05MicrobialProcess(), {
        rvs$crd05_microbial_process <- input$slctin_crd05_microbial_process
        if (rvs$crd05_microbial_process == "Contamination/Removal") {
            update_tooltip("tltp_crd05_microbial_process", "Contamination Maths:\nCarried Contamination (CFU/g) = Carried Contamination (CFU/g) + 10^Modeled Value (log10 CFU/g)")
            enable("slctin_crd05_contamination_basis")
            show("slctin_crd05_contamination_basis")
        } else if (rvs$crd05_microbial_process == "Increase/Reduction") {
            update_tooltip("tltp_crd05_microbial_process", "Increase/Reduction Maths:\nCarried Contamination (CFU/g) = 10^(log10(Carried Contamination (CFU/g)) + Modeled Value (log10 CFU/g))")
            disable("slctin_crd05_contamination_basis")
            hide("slctin_crd05_contamination_basis")
        } else {
            update_tooltip("tltp_crd05_microbial_process", "")
            disable("slctin_crd05_contamination_basis")
            hide("slctin_crd05_contamination_basis")
        }
        if (rvs$crd05_microbial_process == "Risk Output Test" | rvs$crd05_microbial_process == "Product Test") {
            disable("slctin_crd05_distribution")
            hide("slctin_crd05_distribution")
            updateNumericInput(session, "nmcin_crd05_parameter_1", label = "Composite Mass (g)")
            updateNumericInput(session, "nmcin_crd05_parameter_2", label = "Grabs (N)")
            updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Tests (N)")
        } else if (rvs$crd05_microbial_process == "Contamination/Removal") {
            enable("slctin_crd05_distribution")
            show("slctin_crd05_distribution")
            updateSelectInput(session, "slctin_crd05_distribution", choices = list_contamination_distributions, selected = list_contamination_distributions[1])
        } else if (rvs$crd05_microbial_process == "Increase/Reduction") {
            enable("slctin_crd05_distribution")
            show("slctin_crd05_distribution")
            updateSelectInput(session, "slctin_crd05_distribution", choices = list_distributions, selected = list_distributions[1])
        }
    })
    observeEvent(toListenSelectinCard05ProcessAndDist(), {
        rvs$crd05_microbial_process <- input$slctin_crd05_microbial_process
        if (rvs$crd05_microbial_process == "Contamination/Removal") {
            rvs$crd05_distribution <- input$slctin_crd05_distribution
            if (rvs$crd05_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd05_parameter_3")
            } else if (rvs$crd05_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd05_parameter_3")
            } else if (rvs$crd05_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd05_parameter_3")
            } else if (rvs$crd05_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd05_parameter_3")
            } else if (rvs$crd05_distribution == "Normal (linear)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd05_parameter_3")
            } else if (rvs$crd05_distribution == "Uniform (linear)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd05_parameter_3")
            } else if (rvs$crd05_distribution == "Pert (linear)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd05_parameter_3")
            } else if (rvs$crd05_distribution == "Triangular (linear)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd05_parameter_3")
            }
        } else if (rvs$crd05_microbial_process == "Increase/Reduction") {
            rvs$crd05_distribution <- input$slctin_crd05_distribution
            if (rvs$crd05_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = "Mean (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = "SD (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd05_parameter_3")
            } else if (rvs$crd05_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = "Lower Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = "Upper Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd05_parameter_3")
            } else if (rvs$crd05_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd05_parameter_3")
            } else if (rvs$crd05_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd05_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd05_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd05_parameter_3")
            }
        } else {
            shinyjs::show(id = "nmcin_crd05_parameter_3")
        }
    })
    observeEvent(toListenNumberinCard05POccurrence(), {
        rvs$crd05_p_occurrence <- input$nmcin_crd05_p_occurrence
    })
    observeEvent(toListenSelectinCard05ContaminationBasis(), {
        rvs$crd05_contamination_basis <- input$slctin_crd05_contamination_basis
        if (rvs$crd05_contamination_basis == list_contamination_types[1]) {
            rvs$contamination_unit <- ""
        } else if (rvs$crd05_contamination_basis == list_contamination_types[2]) {
            rvs$contamination_unit <- "/g"
        } else if (rvs$crd05_contamination_basis == list_contamination_types[3]) {
            rvs$contamination_unit <- "/lb"
        }
        if (rvs$crd05_distribution == "Normal (log10)") {
            updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Not Used")
        } else if (rvs$crd05_distribution == "Uniform (log10)") {
            updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Not Used")
        } else if (rvs$crd05_distribution == "Pert (log10)") {
            updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
            shinyjs::show(id = "nmcin_crd05_parameter_3")
        } else if (rvs$crd05_distribution == "Triangular (log10)") {
            updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd05_distribution == "Normal (linear)") {
            updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Not Used")
        } else if (rvs$crd05_distribution == "Uniform (linear)") {
            updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_3", label = "Not Used")
        } else if (rvs$crd05_distribution == "Pert (linear)") {
            updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd05_distribution == "Triangular (linear)") {
            updateNumericInput(session, "nmcin_crd05_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd05_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        }
    })
    observeEvent(toListenNumberinCard05Parameter1(), {
        rvs$crd05_parameter_1 <- input$nmcin_crd05_parameter_1 
    })
    observeEvent(toListenNumberinCard05Parameter2(), {
        rvs$crd05_parameter_2 <- input$nmcin_crd05_parameter_2 
    })
    observeEvent(toListenNumberinCard05Parameter3(), {
        rvs$crd05_parameter_3 <- input$nmcin_crd05_parameter_3 
    })
    
    observeEvent(toListenCheckboxCard06(), {
        rvs$crd06 <- input$ckbxin_crd06
    })
    observeEvent(toListenTextinCard06ModelStep(), {
        rvs$crd06_model_step <- input$txtin_crd06_model_step
    })
    observeEvent(toListenSelectinCard06MicrobialProcess(), {
        rvs$crd06_microbial_process <- input$slctin_crd06_microbial_process
        if (rvs$crd06_microbial_process == "Contamination/Removal") {
            update_tooltip("tltp_crd06_microbial_process", "Contamination Maths:\nCarried Contamination (CFU/g) = Carried Contamination (CFU/g) + 10^Modeled Value (log10 CFU/g)")
            enable("slctin_crd06_contamination_basis")
            show("slctin_crd06_contamination_basis")
        } else if (rvs$crd06_microbial_process == "Increase/Reduction") {
            update_tooltip("tltp_crd06_microbial_process", "Increase/Reduction Maths:\nCarried Contamination (CFU/g) = 10^(log10(Carried Contamination (CFU/g)) + Modeled Value (log10 CFU/g))")
            disable("slctin_crd06_contamination_basis")
            hide("slctin_crd06_contamination_basis")
        } else {
            update_tooltip("tltp_crd06_microbial_process", "")
            disable("slctin_crd06_contamination_basis")
            hide("slctin_crd06_contamination_basis")
        }
        if (rvs$crd06_microbial_process == "Risk Output Test" | rvs$crd06_microbial_process == "Product Test") {
            disable("slctin_crd06_distribution")
            hide("slctin_crd06_distribution")
            updateNumericInput(session, "nmcin_crd06_parameter_1", label = "Composite Mass (g)")
            updateNumericInput(session, "nmcin_crd06_parameter_2", label = "Grabs (N)")
            updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Tests (N)")
        } else if (rvs$crd06_microbial_process == "Contamination/Removal") {
            enable("slctin_crd06_distribution")
            show("slctin_crd06_distribution")
            updateSelectInput(session, "slctin_crd06_distribution", choices = list_contamination_distributions, selected = list_contamination_distributions[1])
        } else if (rvs$crd06_microbial_process == "Increase/Reduction") {
            enable("slctin_crd06_distribution")
            show("slctin_crd06_distribution")
            updateSelectInput(session, "slctin_crd06_distribution", choices = list_distributions, selected = list_distributions[1])
        }
    })
    observeEvent(toListenSelectinCard06ProcessAndDist(), {
        rvs$crd06_microbial_process <- input$slctin_crd06_microbial_process
        if (rvs$crd06_microbial_process == "Contamination/Removal") {
            rvs$crd06_distribution <- input$slctin_crd06_distribution
            if (rvs$crd06_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd06_parameter_3")
            } else if (rvs$crd06_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd06_parameter_3")
            } else if (rvs$crd06_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd06_parameter_3")
            } else if (rvs$crd06_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd06_parameter_3")
            } else if (rvs$crd06_distribution == "Normal (linear)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd06_parameter_3")
            } else if (rvs$crd06_distribution == "Uniform (linear)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd06_parameter_3")
            } else if (rvs$crd06_distribution == "Pert (linear)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd06_parameter_3")
            } else if (rvs$crd06_distribution == "Triangular (linear)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd06_parameter_3")
            }
        } else if (rvs$crd06_microbial_process == "Increase/Reduction") {
            rvs$crd06_distribution <- input$slctin_crd06_distribution
            if (rvs$crd06_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = "Mean (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = "SD (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd06_parameter_3")
            } else if (rvs$crd06_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = "Lower Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = "Upper Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd06_parameter_3")
            } else if (rvs$crd06_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd06_parameter_3")
            } else if (rvs$crd06_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd06_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd06_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd06_parameter_3")
            }
        } else {
            shinyjs::show(id = "nmcin_crd06_parameter_3")
        }
    })
    observeEvent(toListenNumberinCard06POccurrence(), {
        rvs$crd06_p_occurrence <- input$nmcin_crd06_p_occurrence
    })
    observeEvent(toListenSelectinCard06ContaminationBasis(), {
        rvs$crd06_contamination_basis <- input$slctin_crd06_contamination_basis
        if (rvs$crd06_contamination_basis == list_contamination_types[1]) {
            rvs$contamination_unit <- ""
        } else if (rvs$crd06_contamination_basis == list_contamination_types[2]) {
            rvs$contamination_unit <- "/g"
        } else if (rvs$crd06_contamination_basis == list_contamination_types[3]) {
            rvs$contamination_unit <- "/lb"
        }
        if (rvs$crd06_distribution == "Normal (log10)") {
            updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Not Used")
        } else if (rvs$crd06_distribution == "Uniform (log10)") {
            updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Not Used")
        } else if (rvs$crd06_distribution == "Pert (log10)") {
            updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
            shinyjs::show(id = "nmcin_crd06_parameter_3")
        } else if (rvs$crd06_distribution == "Triangular (log10)") {
            updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd06_distribution == "Normal (linear)") {
            updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Not Used")
        } else if (rvs$crd06_distribution == "Uniform (linear)") {
            updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_3", label = "Not Used")
        } else if (rvs$crd06_distribution == "Pert (linear)") {
            updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd06_distribution == "Triangular (linear)") {
            updateNumericInput(session, "nmcin_crd06_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd06_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        }
    })
    observeEvent(toListenNumberinCard06Parameter1(), {
        rvs$crd06_parameter_1 <- input$nmcin_crd06_parameter_1 
    })
    observeEvent(toListenNumberinCard06Parameter2(), {
        rvs$crd06_parameter_2 <- input$nmcin_crd06_parameter_2 
    })
    observeEvent(toListenNumberinCard06Parameter3(), {
        rvs$crd06_parameter_3 <- input$nmcin_crd06_parameter_3 
    })
    
    observeEvent(toListenCheckboxCard07(), {
        rvs$crd07 <- input$ckbxin_crd07
    })
    observeEvent(toListenTextinCard07ModelStep(), {
        rvs$crd07_model_step <- input$txtin_crd07_model_step
    })
    observeEvent(toListenSelectinCard07MicrobialProcess(), {
        rvs$crd07_microbial_process <- input$slctin_crd07_microbial_process
        if (rvs$crd07_microbial_process == "Contamination/Removal") {
            update_tooltip("tltp_crd07_microbial_process", "Contamination Maths:\nCarried Contamination (CFU/g) = Carried Contamination (CFU/g) + 10^Modeled Value (log10 CFU/g)")
            enable("slctin_crd07_contamination_basis")
            show("slctin_crd07_contamination_basis")
        } else if (rvs$crd07_microbial_process == "Increase/Reduction") {
            update_tooltip("tltp_crd07_microbial_process", "Increase/Reduction Maths:\nCarried Contamination (CFU/g) = 10^(log10(Carried Contamination (CFU/g)) + Modeled Value (log10 CFU/g))")
            disable("slctin_crd07_contamination_basis")
            hide("slctin_crd07_contamination_basis")
        } else {
            update_tooltip("tltp_crd07_microbial_process", "")
            disable("slctin_crd07_contamination_basis")
            hide("slctin_crd07_contamination_basis")
        }
        if (rvs$crd07_microbial_process == "Risk Output Test" | rvs$crd07_microbial_process == "Product Test") {
            disable("slctin_crd07_distribution")
            hide("slctin_crd07_distribution")
            updateNumericInput(session, "nmcin_crd07_parameter_1", label = "Composite Mass (g)")
            updateNumericInput(session, "nmcin_crd07_parameter_2", label = "Grabs (N)")
            updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Tests (N)")
        } else if (rvs$crd07_microbial_process == "Contamination/Removal") {
            enable("slctin_crd07_distribution")
            show("slctin_crd07_distribution")
            updateSelectInput(session, "slctin_crd07_distribution", choices = list_contamination_distributions, selected = list_contamination_distributions[1])
        } else if (rvs$crd07_microbial_process == "Increase/Reduction") {
            enable("slctin_crd07_distribution")
            show("slctin_crd07_distribution")
            updateSelectInput(session, "slctin_crd07_distribution", choices = list_distributions, selected = list_distributions[1])
        }
    })
    observeEvent(toListenSelectinCard07ProcessAndDist(), {
        rvs$crd07_microbial_process <- input$slctin_crd07_microbial_process
        if (rvs$crd07_microbial_process == "Contamination/Removal") {
            rvs$crd07_distribution <- input$slctin_crd07_distribution
            if (rvs$crd07_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd07_parameter_3")
            } else if (rvs$crd07_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd07_parameter_3")
            } else if (rvs$crd07_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd07_parameter_3")
            } else if (rvs$crd07_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd07_parameter_3")
            } else if (rvs$crd07_distribution == "Normal (linear)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd07_parameter_3")
            } else if (rvs$crd07_distribution == "Uniform (linear)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd07_parameter_3")
            } else if (rvs$crd07_distribution == "Pert (linear)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd07_parameter_3")
            } else if (rvs$crd07_distribution == "Triangular (linear)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd07_parameter_3")
            }
        } else if (rvs$crd07_microbial_process == "Increase/Reduction") {
            rvs$crd07_distribution <- input$slctin_crd07_distribution
            if (rvs$crd07_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = "Mean (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = "SD (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd07_parameter_3")
            } else if (rvs$crd07_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = "Lower Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = "Upper Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd07_parameter_3")
            } else if (rvs$crd07_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd07_parameter_3")
            } else if (rvs$crd07_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd07_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd07_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd07_parameter_3")
            }
        } else {
            shinyjs::show(id = "nmcin_crd07_parameter_3")
        }
    })
    observeEvent(toListenNumberinCard07POccurrence(), {
        rvs$crd07_p_occurrence <- input$nmcin_crd07_p_occurrence
    })
    observeEvent(toListenSelectinCard07ContaminationBasis(), {
        rvs$crd07_contamination_basis <- input$slctin_crd07_contamination_basis
        if (rvs$crd07_contamination_basis == list_contamination_types[1]) {
            rvs$contamination_unit <- ""
        } else if (rvs$crd07_contamination_basis == list_contamination_types[2]) {
            rvs$contamination_unit <- "/g"
        } else if (rvs$crd07_contamination_basis == list_contamination_types[3]) {
            rvs$contamination_unit <- "/lb"
        }
        if (rvs$crd07_distribution == "Normal (log10)") {
            updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Not Used")
        } else if (rvs$crd07_distribution == "Uniform (log10)") {
            updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Not Used")
        } else if (rvs$crd07_distribution == "Pert (log10)") {
            updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
            shinyjs::show(id = "nmcin_crd07_parameter_3")
        } else if (rvs$crd07_distribution == "Triangular (log10)") {
            updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd07_distribution == "Normal (linear)") {
            updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Not Used")
        } else if (rvs$crd07_distribution == "Uniform (linear)") {
            updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_3", label = "Not Used")
        } else if (rvs$crd07_distribution == "Pert (linear)") {
            updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd07_distribution == "Triangular (linear)") {
            updateNumericInput(session, "nmcin_crd07_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd07_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        }
    })
    observeEvent(toListenNumberinCard07Parameter1(), {
        rvs$crd07_parameter_1 <- input$nmcin_crd07_parameter_1 
    })
    observeEvent(toListenNumberinCard07Parameter2(), {
        rvs$crd07_parameter_2 <- input$nmcin_crd07_parameter_2 
    })
    observeEvent(toListenNumberinCard07Parameter3(), {
        rvs$crd07_parameter_3 <- input$nmcin_crd07_parameter_3 
    })
    
    observeEvent(toListenCheckboxCard08(), {
        rvs$crd08 <- input$ckbxin_crd08
    })
    observeEvent(toListenTextinCard08ModelStep(), {
        rvs$crd08_model_step <- input$txtin_crd08_model_step
    })
    observeEvent(toListenSelectinCard08MicrobialProcess(), {
        rvs$crd08_microbial_process <- input$slctin_crd08_microbial_process
        if (rvs$crd08_microbial_process == "Contamination/Removal") {
            update_tooltip("tltp_crd08_microbial_process", "Contamination Maths:\nCarried Contamination (CFU/g) = Carried Contamination (CFU/g) + 10^Modeled Value (log10 CFU/g)")
            enable("slctin_crd08_contamination_basis")
            show("slctin_crd08_contamination_basis")
        } else if (rvs$crd08_microbial_process == "Increase/Reduction") {
            update_tooltip("tltp_crd08_microbial_process", "Increase/Reduction Maths:\nCarried Contamination (CFU/g) = 10^(log10(Carried Contamination (CFU/g)) + Modeled Value (log10 CFU/g))")
            disable("slctin_crd08_contamination_basis")
            hide("slctin_crd08_contamination_basis")
        } else {
            update_tooltip("tltp_crd08_microbial_process", "")
            disable("slctin_crd08_contamination_basis")
            hide("slctin_crd08_contamination_basis")
        }
        if (rvs$crd08_microbial_process == "Risk Output Test" | rvs$crd08_microbial_process == "Product Test") {
            disable("slctin_crd08_distribution")
            hide("slctin_crd08_distribution")
            updateNumericInput(session, "nmcin_crd08_parameter_1", label = "Composite Mass (g)")
            updateNumericInput(session, "nmcin_crd08_parameter_2", label = "Grabs (N)")
            updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Tests (N)")
        } else if (rvs$crd08_microbial_process == "Contamination/Removal") {
            enable("slctin_crd08_distribution")
            show("slctin_crd08_distribution")
            updateSelectInput(session, "slctin_crd08_distribution", choices = list_contamination_distributions, selected = list_contamination_distributions[1])
        } else if (rvs$crd08_microbial_process == "Increase/Reduction") {
            enable("slctin_crd08_distribution")
            show("slctin_crd08_distribution")
            updateSelectInput(session, "slctin_crd08_distribution", choices = list_distributions, selected = list_distributions[1])
        }
    })
    observeEvent(toListenSelectinCard08ProcessAndDist(), {
        rvs$crd08_microbial_process <- input$slctin_crd08_microbial_process
        if (rvs$crd08_microbial_process == "Contamination/Removal") {
            rvs$crd08_distribution <- input$slctin_crd08_distribution
            if (rvs$crd08_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd08_parameter_3")
            } else if (rvs$crd08_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd08_parameter_3")
            } else if (rvs$crd08_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd08_parameter_3")
            } else if (rvs$crd08_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd08_parameter_3")
            } else if (rvs$crd08_distribution == "Normal (linear)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd08_parameter_3")
            } else if (rvs$crd08_distribution == "Uniform (linear)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd08_parameter_3")
            } else if (rvs$crd08_distribution == "Pert (linear)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd08_parameter_3")
            } else if (rvs$crd08_distribution == "Triangular (linear)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd08_parameter_3")
            }
        } else if (rvs$crd08_microbial_process == "Increase/Reduction") {
            rvs$crd08_distribution <- input$slctin_crd08_distribution
            if (rvs$crd08_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = "Mean (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = "SD (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd08_parameter_3")
            } else if (rvs$crd08_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = "Lower Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = "Upper Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd08_parameter_3")
            } else if (rvs$crd08_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd08_parameter_3")
            } else if (rvs$crd08_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd08_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd08_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd08_parameter_3")
            }
        } else {
            shinyjs::show(id = "nmcin_crd08_parameter_3")
        }
    })
    observeEvent(toListenNumberinCard08POccurrence(), {
        rvs$crd08_p_occurrence <- input$nmcin_crd08_p_occurrence
    })
    observeEvent(toListenSelectinCard08ContaminationBasis(), {
        rvs$crd08_contamination_basis <- input$slctin_crd08_contamination_basis
        if (rvs$crd08_contamination_basis == list_contamination_types[1]) {
            rvs$contamination_unit <- ""
        } else if (rvs$crd08_contamination_basis == list_contamination_types[2]) {
            rvs$contamination_unit <- "/g"
        } else if (rvs$crd08_contamination_basis == list_contamination_types[3]) {
            rvs$contamination_unit <- "/lb"
        }
        if (rvs$crd08_distribution == "Normal (log10)") {
            updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Not Used")
        } else if (rvs$crd08_distribution == "Uniform (log10)") {
            updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Not Used")
        } else if (rvs$crd08_distribution == "Pert (log10)") {
            updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
            shinyjs::show(id = "nmcin_crd08_parameter_3")
        } else if (rvs$crd08_distribution == "Triangular (log10)") {
            updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd08_distribution == "Normal (linear)") {
            updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Not Used")
        } else if (rvs$crd08_distribution == "Uniform (linear)") {
            updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_3", label = "Not Used")
        } else if (rvs$crd08_distribution == "Pert (linear)") {
            updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd08_distribution == "Triangular (linear)") {
            updateNumericInput(session, "nmcin_crd08_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd08_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        }
    })
    observeEvent(toListenNumberinCard08Parameter1(), {
        rvs$crd08_parameter_1 <- input$nmcin_crd08_parameter_1 
    })
    observeEvent(toListenNumberinCard08Parameter2(), {
        rvs$crd08_parameter_2 <- input$nmcin_crd08_parameter_2 
    })
    observeEvent(toListenNumberinCard08Parameter3(), {
        rvs$crd08_parameter_3 <- input$nmcin_crd08_parameter_3 
    })
    
    observeEvent(toListenCheckboxCard09(), {
        rvs$crd09 <- input$ckbxin_crd09
    })
    observeEvent(toListenTextinCard09ModelStep(), {
        rvs$crd09_model_step <- input$txtin_crd09_model_step
    })
    observeEvent(toListenSelectinCard09MicrobialProcess(), {
        rvs$crd09_microbial_process <- input$slctin_crd09_microbial_process
        if (rvs$crd09_microbial_process == "Contamination/Removal") {
            update_tooltip("tltp_crd09_microbial_process", "Contamination Maths:\nCarried Contamination (CFU/g) = Carried Contamination (CFU/g) + 10^Modeled Value (log10 CFU/g)")
            enable("slctin_crd09_contamination_basis")
            show("slctin_crd09_contamination_basis")
        } else if (rvs$crd09_microbial_process == "Increase/Reduction") {
            update_tooltip("tltp_crd09_microbial_process", "Increase/Reduction Maths:\nCarried Contamination (CFU/g) = 10^(log10(Carried Contamination (CFU/g)) + Modeled Value (log10 CFU/g))")
            disable("slctin_crd09_contamination_basis")
            hide("slctin_crd09_contamination_basis")
        } else {
            update_tooltip("tltp_crd09_microbial_process", "")
            disable("slctin_crd09_contamination_basis")
            hide("slctin_crd09_contamination_basis")
        }
        if (rvs$crd09_microbial_process == "Risk Output Test" | rvs$crd09_microbial_process == "Product Test") {
            disable("slctin_crd09_distribution")
            hide("slctin_crd09_distribution")
            updateNumericInput(session, "nmcin_crd09_parameter_1", label = "Composite Mass (g)")
            updateNumericInput(session, "nmcin_crd09_parameter_2", label = "Grabs (N)")
            updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Tests (N)")
        } else if (rvs$crd09_microbial_process == "Contamination/Removal") {
            enable("slctin_crd09_distribution")
            show("slctin_crd09_distribution")
            updateSelectInput(session, "slctin_crd09_distribution", choices = list_contamination_distributions, selected = list_contamination_distributions[1])
        } else if (rvs$crd09_microbial_process == "Increase/Reduction") {
            enable("slctin_crd09_distribution")
            show("slctin_crd09_distribution")
            updateSelectInput(session, "slctin_crd09_distribution", choices = list_distributions, selected = list_distributions[1])
        }
    })
    observeEvent(toListenSelectinCard09ProcessAndDist(), {
        rvs$crd09_microbial_process <- input$slctin_crd09_microbial_process
        if (rvs$crd09_microbial_process == "Contamination/Removal") {
            rvs$crd09_distribution <- input$slctin_crd09_distribution
            if (rvs$crd09_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd09_parameter_3")
            } else if (rvs$crd09_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd09_parameter_3")
            } else if (rvs$crd09_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd09_parameter_3")
            } else if (rvs$crd09_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd09_parameter_3")
            } else if (rvs$crd09_distribution == "Normal (linear)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd09_parameter_3")
            } else if (rvs$crd09_distribution == "Uniform (linear)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd09_parameter_3")
            } else if (rvs$crd09_distribution == "Pert (linear)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd09_parameter_3")
            } else if (rvs$crd09_distribution == "Triangular (linear)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
                shinyjs::show(id = "nmcin_crd09_parameter_3")
            }
        } else if (rvs$crd09_microbial_process == "Increase/Reduction") {
            rvs$crd09_distribution <- input$slctin_crd09_distribution
            if (rvs$crd09_distribution == "Normal (log10)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = "Mean (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = "SD (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd09_parameter_3")
            } else if (rvs$crd09_distribution == "Uniform (log10)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = "Lower Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = "Upper Bound (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Not Used")
                shinyjs::hide(id = "nmcin_crd09_parameter_3")
            } else if (rvs$crd09_distribution == "Pert (log10)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd09_parameter_3")
            } else if (rvs$crd09_distribution == "Triangular (log10)") {
                updateNumericInput(session, "nmcin_crd09_parameter_1", label = "Minimum (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd09_parameter_2", label = "Mode (log10 CFU/g)")
                updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Maximum (log10 CFU/g)")
                shinyjs::show(id = "nmcin_crd09_parameter_3")
            }
        } else {
            shinyjs::show(id = "nmcin_crd09_parameter_3")
        }
    })
    observeEvent(toListenNumberinCard09POccurrence(), {
        rvs$crd09_p_occurrence <- input$nmcin_crd09_p_occurrence
    })
    observeEvent(toListenSelectinCard09ContaminationBasis(), {
        rvs$crd09_contamination_basis <- input$slctin_crd09_contamination_basis
        if (rvs$crd09_contamination_basis == list_contamination_types[1]) {
            rvs$contamination_unit <- ""
        } else if (rvs$crd09_contamination_basis == list_contamination_types[2]) {
            rvs$contamination_unit <- "/g"
        } else if (rvs$crd09_contamination_basis == list_contamination_types[3]) {
            rvs$contamination_unit <- "/lb"
        }
        if (rvs$crd09_distribution == "Normal (log10)") {
            updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Mean (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("SD (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Not Used")
        } else if (rvs$crd09_distribution == "Uniform (log10)") {
            updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Lower Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Upper Bound (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Not Used")
        } else if (rvs$crd09_distribution == "Pert (log10)") {
            updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
            shinyjs::show(id = "nmcin_crd09_parameter_3")
        } else if (rvs$crd09_distribution == "Triangular (log10)") {
            updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Minimum (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Mode (log10 CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_3", label = paste0("Maximum (log10 CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd09_distribution == "Normal (linear)") {
            updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Mean (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("SD (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Not Used")
        } else if (rvs$crd09_distribution == "Uniform (linear)") {
            updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Lower Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Upper Bound (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_3", label = "Not Used")
        } else if (rvs$crd09_distribution == "Pert (linear)") {
            updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        } else if (rvs$crd09_distribution == "Triangular (linear)") {
            updateNumericInput(session, "nmcin_crd09_parameter_1", label = paste0("Minimum (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_2", label = paste0("Mode (CFU", rvs$contamination_unit, ")"))
            updateNumericInput(session, "nmcin_crd09_parameter_3", label = paste0("Maximum (CFU", rvs$contamination_unit, ")"))
        }
    })
    observeEvent(toListenNumberinCard09Parameter1(), {
        rvs$crd09_parameter_1 <- input$nmcin_crd09_parameter_1 
    })
    observeEvent(toListenNumberinCard09Parameter2(), {
        rvs$crd09_parameter_2 <- input$nmcin_crd09_parameter_2 
    })
    observeEvent(toListenNumberinCard09Parameter3(), {
        rvs$crd09_parameter_3 <- input$nmcin_crd09_parameter_3 
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
