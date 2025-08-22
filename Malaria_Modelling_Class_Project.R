###### The following exercise 
# Network Analysis and SIS Epidemic Model
# Load required libraries
library(igraph)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Create the dataset
data <- data.frame(
  Participant = c("KAM", "SIN", "RAH", "DAW", "LUC", "TUH", "BAS", "MBA", "MUS", "DAV", "CHA", "EMM", "LYD", "VAL", "GOD", "AGN", "ROS", "FRA", "JAM", "DGB", "ELI", "WAN", "MUT", "GEO", "GEF", "SAH", "PAT", "TAB", "BRI", "STA", "OUM", "STE", "THU", "EMI", "CAM", "ZEN", "BIL", "JOH", "ALI", "PUR", "WAI", "MIL", "BRL", "TRI"), 
  Country = c("Cameroon", "Ethiopia", "Tanzania", "Ethiopia", "Zambia", "Uganda", "Nigeria", "Cameroon", "Zimbabwe", "Zambia", "Kenya", "Kenya", "Kenya", "Kenya", "Malawi", "Malawi", "Tanzania", "Tanzania", "Kenya", "Ethiopia", "Tanzania", "Kenya", "Kenya", "Kenya", "Kenya", "Tunisia", "Kenya", "Kenya", "Kenya", "Kenya", "Tunisia", "Kenya", "Kenya", "Switzerland", "Switzerland", "Ghana", "Switzerland", "Kenya", "Kenya", "Kenya", "Kenya", "Kenya", "Kenya", "Kenya")
)

#####PART 1. NETWORK CONSTRUCTION

##Remember Networks are based on matrices - adjacency matrix defines who is conneted to whom.
# Create adjacency matrix based on connection rules: a) People from same country are connected; b) Non-Kenyans are also connected to Kenyans. 
create_network <- function(data) {
  n <- nrow(data)
  adj_matrix <- matrix(0, n, n)
  rownames(adj_matrix) <- colnames(adj_matrix) <- data$Participant
  
  for(i in 1:n) {
    for(j in 1:n) {
      if(i != j) {
        ##Same country connect
        if(data$Country[i] == data$Country[j]) {
          adj_matrix[i,j] <- 1
        }
        ##Non-Kenyans connect to Kenyans
        else if(data$Country[i] != "Kenya" & data$Country[j] == "Kenya") {
          adj_matrix[i,j] <- 1
        }
        else if(data$Country[i] == "Kenya" & data$Country[j] != "Kenya") {
          adj_matrix[i,j] <- 1
        }
      }
    }
  }
  return(adj_matrix)
}

# Create the network and add attributes (CORRECTED)
adj_matrix <- create_network(data) ### Explore the adjacency 
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

# Add attributes - CORRECTED VERSION
V(g)$country <- data$Country
V(g)$name <- data$Participant

##Create and visualize the network
plot_network <- function(g, data) {
  ##Color nodes by country
  countries <- unique(data$Country)
  colors <- RColorBrewer::brewer.pal(min(length(countries), 11), "Spectral")
  country_colors <- setNames(colors[1:length(countries)], countries)
  
  ##Set vertex colors using the country attribute
  V(g)$color <- country_colors[V(g)$country]
  
  ##Set the graph layout - there are many layouts - please explore them
  set.seed(123)
  layout <- layout_with_fr(g)
  
  # Plot
  plot(g, 
       layout = layout,
       vertex.size = 8,
       vertex.label = V(g)$name,  # Explicitly use vertex names
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       vertex.frame.color = "white",
       edge.color = "gray70",
       edge.width = 0.8,
       main = "Participant Interaction Network")
  
  ##Add legend - you can change position or text size
  legend("topright", legend = countries, fill = country_colors, cex = 0.8, title = "Country")
}

###Plot the network
plot_network(g, data)

### Now you have a network - 

# 2. SIS EPIDEMIC MODdata = # 2. SIS EPIDEMIC MODEL FUNCTION
sis_model <- function(adj_matrix, beta, gamma, initial_infected, days, superspreaders = NULL, mask_day = NULL, mask_efficacy = 0.5) {
  n <- nrow(adj_matrix)
  participants <- rownames(adj_matrix)
  
  # Initialize states: 0 = Susceptible, 1 = Infected
  states <- matrix(0, nrow = days + 1, ncol = n)
  colnames(states) <- participants
  
  # Set initial infected
  states[1, initial_infected] <- 1
  
  # Track daily infections for incidence calculation
  daily_new_infections <- rep(0, days + 1)
  
  # Superspreader multiplier
  ss_multiplier <- rep(1, n)
  if(!is.null(superspreaders)) {
    ss_idx <- which(participants %in% superspreaders)
    ss_multiplier[ss_idx] <- 3  # 3x transmission rate
  }
  
  for(day in 1:days) {
    current_state <- states[day, ]
    new_state <- current_state
    new_infections <- 0
    
    # Apply mask effect if applicable
    current_beta <- beta
    if(!is.null(mask_day) && day >= mask_day) {
      current_beta <- beta * (1 - mask_efficacy)
    }
    
    for(i in 1:n) {
      if(current_state[i] == 0) {  # If susceptible
        # Calculate infection probability from infected neighbors
        infected_neighbors <- which(current_state == 1 & adj_matrix[i, ] == 1)
        if(length(infected_neighbors) > 0) {
          # Calculate total force of infection
          total_force <- 0
          for(neighbor in infected_neighbors) {
            total_force <- total_force + current_beta * ss_multiplier[neighbor]
          }
          
          # Probability of infection (1 - probability of avoiding all infections)
          prob_infection <- 1 - exp(-total_force)
          
          if(runif(1) < prob_infection) {
            new_state[i] <- 1
            new_infections <- new_infections + 1
          }
        }
      } else {  # If infected
        # Recovery
        if(runif(1) < gamma) {
          new_state[i] <- 0
        }
      }
    }
    
    states[day + 1, ] <- new_state
    daily_new_infections[day + 1] <- new_infections
  }
  
  return(list(states = states, daily_new_infections = daily_new_infections))
}

# 3. RUN SIMULATIONS

# Set parameters
beta <- 0.4
gamma <- 0.1
days <- 8

# Find participant indices
zeinabu_idx <- which(data$Participant == "ZEN")
agnes_idx <- which(data$Participant == "AGN")
initial_infected <- c("ZEN", "AGN")

# Scenario 1: Basic SIS model
set.seed(42)
result1 <- sis_model(adj_matrix, beta, gamma, initial_infected, days)

# Scenario 2: With masks from day 4
set.seed(42)
result2 <- sis_model(adj_matrix, beta, gamma, initial_infected, days, mask_day = 4)

# Scenario 3: With superspreaders
superspreaders <- c("BAS", "BIL", "MUT", "PUR")
set.seed(42)
result3 <- sis_model(adj_matrix, beta, gamma, initial_infected, days, superspreaders = superspreaders)

# 4. ANALYSIS AND RESULTS

cat("=== EPIDEMIC ANALYSIS RESULTS ===\n\n")

# a) Number of infected on days 2, 4, 8
infected_day2 <- sum(result1$states[3, ])  # Day 2 (row 3, 0-indexed)
infected_day4 <- sum(result1$states[5, ])  # Day 4
infected_day8 <- sum(result1$states[9, ])  # Day 8

cat("a) Number of infected individuals:\n")
cat("   Day 2:", infected_day2, "people\n")
cat("   Day 4:", infected_day4, "people\n")
cat("   Day 8:", infected_day8, "people\n\n")

# b) Incidence at day 8
incidence_day8 <- result1$daily_new_infections[9]
cat("b) Incidence on Day 8:", incidence_day8, "new infections\n\n")

# c) Infections averted with masks
infected_day8_masks <- sum(result2$states[9, ])
infections_averted <- infected_day8 - infected_day8_masks
cat("c) With masks from Day 4:\n")
cat("   Infections on Day 8 with masks:", infected_day8_masks, "\n")
cat("   Infections averted:", infections_averted, "\n\n")

# d) Network metrics
cat("d) Network Metrics:\n")

# i) Density
density_val <- edge_density(g)
cat("   i) Network density:", round(density_val, 3), "\n")

# ii) Superspreader analysis
cat("   ii) Superspreader analysis:\n")
infected_normal <- sum(result1$states[9, ])
infected_superspreaders <- sum(result3$states[9, ])
cat("       Without superspreaders: ", infected_normal, " infected\n")
cat("       With superspreaders: ", infected_superspreaders, " infected\n")
cat("       Additional infections from superspreaders: ", infected_superspreaders - infected_normal, "\n")

# iii) Effect of interaction rate changes
cat("   iii) Effect of changing interaction rates:\n")
betas <- c(0.2, 0.3, 0.4, 0.5, 0.6)
final_infected <- numeric(length(betas))

for(i in seq_along(betas)) {
  set.seed(42)
  temp_result <- sis_model(adj_matrix, betas[i], gamma, initial_infected, days)
  final_infected[i] <- sum(temp_result$states[9, ])
}

beta_results <- data.frame(Beta = betas, Final_Infected = final_infected)
print(beta_results)

# 5. VISUALIZATION OF EPIDEMIC PROGRESSION
plot_epidemic <- function(states) {
  days <- nrow(states) - 1
  infected_counts <- rowSums(states)
  
  df <- data.frame(
    Day = 0:days,
    Infected = infected_counts
  )
  
  ggplot(df, aes(x = Day, y = Infected)) +
    geom_line(color = "red", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(title = "Epidemic Progression",
         x = "Day",
         y = "Number of Infected Individuals") +
    theme_minimal() +
    scale_x_continuous(breaks = 0:days)
}

# Plot epidemic progression
print(plot_epidemic(result1$states))

# Comparison plot
compare_df <- data.frame(
  Day = rep(0:days, 3),
  Infected = c(rowSums(result1$states), rowSums(result2$states), rowSums(result3$states)),
  Scenario = rep(c("Basic", "With Masks", "With Superspreaders"), each = days + 1)
)

comparison_plot <- ggplot(compare_df, aes(x = Day, y = Infected, color = Scenario)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Comparison of Epidemic Scenarios",
       x = "Day",
       y = "Number of Infected Individuals") +
  theme_minimal() +
  scale_x_continuous(breaks = 0:days) +
  scale_color_manual(values = c("red", "blue", "orange"))

print(comparison_plot)

cat("\n=== SUMMARY ===\n")
cat("Network has", vcount(g), "participants and", ecount(g), "connections\n")
cat("Density:", round(density_val, 3), "- indicates a", 
    ifelse(density_val > 0.5, "highly", "moderately"), "connected network\n")
cat("Superspreaders significantly increase transmission\n")
cat("Masks can effectively reduce transmission when implemented\n")