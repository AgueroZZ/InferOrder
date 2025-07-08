# Function to recover the ordering based on OTSP
recover_ordering_otsp <- function(data, method = "nn") {
  # Compute the distance matrix
  dist_mat <- as.matrix(dist(data))

  # Identify the two most distant points
  far_points <- which(dist_mat == max(dist_mat), arr.ind = TRUE)[1, ]
  start_node <- far_points[1]
  end_node <- far_points[2]

  # Add a dummy node
  n <- nrow(dist_mat)
  dummy_row <- rep(Inf, n)
  dummy_col <- rep(Inf, n + 1)
  dist_mat <- rbind(cbind(dist_mat, dummy_row), dummy_col)

  # Connect the dummy node to the start and end points with zero cost
  dist_mat[n + 1, start_node] <- 0
  dist_mat[start_node, n + 1] <- 0
  dist_mat[n + 1, end_node] <- 0
  dist_mat[end_node, n + 1] <- 0

  # Ensure infinite cost for dummy-to-dummy transitions
  dist_mat[n + 1, n + 1] <- Inf

  # Solve the TSP
  tsp <- TSP::ATSP(dist_mat)
  tour <- TSP::solve_TSP(tsp, start = (n+1), method = method, rep = 10)

  # Wrap the order around the dummy node
  tour_order <- as.integer(tour)
  dummy_index <- which(tour_order == (n + 1))
  wrapped_order <- c(tour_order[dummy_index:length(tour_order)],
                     tour_order[1:(dummy_index - 1)])

  # Remove the dummy node from the order
  open_path <- wrapped_order[wrapped_order != (n + 1)]

  return(open_path)
}

# Function to recover the ordering based on OTSP by group
# 'groups' should be a vector of group assignments for each observation
# Output: A list of orderings, one for each group
recover_ordering_otsp_grouped <- function(data, method = "nn", groups) {
  if (length(groups) != nrow(data)) {
    stop("'groups' must be the same length as the number of rows in 'data'")
  }

  unique_groups <- unique(groups)
  group_orderings <- list()

  for (group in unique_groups) {
    # Subset the data for the current group
    group_data <- data[groups == group, , drop = FALSE]

    # Recover the ordering for the current group
    group_ordering <- recover_ordering_otsp(group_data, method = method)

    # Store the group-specific ordering (using the group-specific row indices)
    group_orderings[[group]] <- group_ordering
  }

  return(group_orderings)
}

# Function to create a structure plot
plot_loadings <- function(loadings, order, which_to_plot = NULL) {
  # Ensure the order is valid
  if (length(order) != nrow(loadings)) {
    stop("The length of 'order' must match the number of rows in the 'loadings' matrix.")
  }

  # Reorder the loadings based on the specified order
  loadings_ordered <- loadings[order, , drop = FALSE]

  # Convert the loadings matrix into a tidy data frame
  loadings_df <- as.data.frame(loadings_ordered)

  # rename the columns
  colnames(loadings_df) <- paste0("loading", 1:ncol(loadings_df))

  # check which factor to plot
  if (is.null(which_to_plot)) {
    loadings_df <- loadings_df[,1:ncol(loadings_df)]
  } else{
    loadings_df <- loadings_df[,which_to_plot]
  }

  loadings_df$x <- seq_len(nrow(loadings_df))

  loadings_long <- pivot_longer(loadings_df, cols = starts_with("loading"),
                                names_to = "Component", values_to = "Value")



  # Plot the structure plot
  ggplot(loadings_long, aes(x = x, y = Value, fill = Component)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_x_continuous("Observation", breaks = seq(0, nrow(loadings), by = 10)) +
    scale_y_continuous("Loading Value") +
    theme_minimal() +
    ggtitle("Structure Plot of Loadings")
}

# Function to plot loadings by group
plot_loadings_grouped <- function(loadings, order_list, groups, which_to_plot = NULL) {
  # Check that inputs are valid
  if (length(groups) != nrow(loadings)) {
    stop("'groups' must have the same length as the number of rows in 'loadings'")
  }

  if (length(order_list) != length(unique(groups))) {
    stop("'order_list' must have one ordering for each unique group in 'groups'")
  }

  # Initialize a list to store reordered loadings
  reordered_loadings_list <- list()

  # Reorder loadings for each group
  unique_groups <- unique(groups)
  for (group in unique_groups) {
    # Subset loadings and order based on the group's ordering
    group_indices <- which(groups == group)
    group_order <- order_list[[group]]
    reordered_loadings <- loadings[group_indices[group_order], , drop = FALSE]

    # Store reordered loadings in the list with group label
    reordered_loadings_list[[group]] <- reordered_loadings
  }

  # Combine reordered loadings into one data frame for plotting
  loadings_combined <- do.call(rbind, reordered_loadings_list)
  group_labels <- unlist(mapply(rep, unique_groups, sapply(reordered_loadings_list, nrow)))
  loadings_combined <- as.data.frame(loadings_combined)
  colnames(loadings_combined) <- paste0("loading", 1:ncol(loadings_combined))

  # check which factor to plot
  if (is.null(which_to_plot)) {
    loadings_combined <- loadings_combined[,1:ncol(loadings_combined)]
  } else{
    loadings_combined <- loadings_combined[,which_to_plot]
  }

  loadings_combined$Group <- factor(group_labels)
  loadings_combined$x <- seq_len(nrow(loadings_combined))

  # Convert to long format for ggplot
  loadings_long <- pivot_longer(loadings_combined, cols = starts_with("loading"),
                                names_to = "Component", values_to = "Value")

  # Plot the structure plot grouped by Group
  ggplot(loadings_long, aes(x = x, y = Value, fill = Component)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_x_continuous("Observation", breaks = seq(0, nrow(loadings), by = 10)) +
    scale_y_continuous("Loading Value") +
    facet_wrap(~Group, scales = "free_x", nrow = 1) +  # One facet per group
    theme_minimal() +
    ggtitle("Structure Plot of Loadings by Group")
}

# Function to compute the loss function for a given order with a customizable distance metric
compute_loss <- function(loadings, order, distance_func = function(x, y) sqrt(sum((x - y)^2)), normalized = FALSE) {
  # Ensure the order is valid
  if (length(order) != nrow(loadings)) {
    stop("The length of 'order' must match the number of rows in the 'loadings' matrix.")
  }

  # Reorder the loadings based on the specified order
  loadings_ordered <- loadings[order, , drop = FALSE]

  # Normalize the loadings so each row summing to 1 if requested
  if (normalized) {
    if(ncol(loadings_ordered) == 1) {
      return(0)
    }
    row_sums <- rowSums(loadings_ordered)
    if(any(row_sums == 0)) {
      stop("Please rerun the function with normalized = FALSE")
    }
    loadings_ordered <- sweep(loadings_ordered, 1, row_sums, FUN = "/")
  }

  # Compute pairwise distances using the provided distance function
  total_distance <- 0
  for (i in seq_len(nrow(loadings_ordered) - 1)) {
    dist <- distance_func(loadings_ordered[i + 1, ], loadings_ordered[i, ])
    total_distance <- total_distance + dist
  }

  return(total_distance)
}

# Function to compute the loss function for a given order_list with optional group weighting
compute_loss_grouped <- function(loadings, order_list, groups,
                                 distance_func = function(x, y) sqrt(sum((x - y)^2)),
                                 normalized = FALSE,
                                 group_weight = FALSE) {
  # Check inputs
  if (!is.matrix(loadings)) stop("'loadings' must be a matrix")
  if (length(groups) != nrow(loadings)) stop("'groups' must have the same length as the number of rows in 'loadings'")
  if (!is.list(order_list)) stop("'order_list' must be a list")

  # Initialize variables
  group_losses <- numeric(0)
  group_sizes <- numeric(0)

  # Compute loss for each group
  for (group in names(order_list)) {
    # Get group-specific data and ordering
    group_data <- loadings[groups == group, , drop = FALSE]
    group_order <- order_list[[group]]

    # Ensure the group ordering matches the group data
    if (length(group_order) != nrow(group_data)) {
      stop(paste("Order length does not match the number of rows for group", group))
    }

    # Compute the loss for this group
    group_loss <- compute_loss(group_data, group_order, distance_func = distance_func, normalized = normalized)
    group_losses <- c(group_losses, group_loss)
    group_sizes <- c(group_sizes, nrow(group_data))
  }

  # Compute weighted or unweighted average loss
  if (group_weight) {
    total_loss <- sum(group_losses * group_sizes) / sum(group_sizes)  # Weighted average
  } else {
    total_loss <- mean(group_losses)  # Unweighted mean
  }

  return(total_loss)
}


# Function to generate an elbow plot for factor selection
elbow_plot_selection <- function(loadings, start_index = 1, plot = TRUE, normalized = FALSE, method = "farthest_insertion", fast = FALSE) {
  # Check inputs
  if (!is.matrix(loadings)) stop("'loadings' must be a matrix")

  n_factors <- ncol(loadings)
  remaining_indices <- seq_len(n_factors)
  selected_indices <- numeric(0)
  loss_values <- numeric(0)

  # if fast = TRUE, then random sample 3000 rows to compute the loss when the number of rows is larger than 3000
  if (fast & nrow(loadings) > 3000) {
    loadings <- loadings[sample(nrow(loadings), 3000),]
  }

  # Helper function to compute loss for a given subset of loadings
  compute_subset_loss <- function(indices) {
    subset_loadings <- loadings[, indices, drop = FALSE]
    ordering <- recover_ordering_otsp(subset_loadings, method = method)
    compute_loss(subset_loadings, ordering, normalized = normalized)
  }

  # Add the start_index as the first selected factor
  selected_indices <- c(start_index)
  remaining_indices <- setdiff(remaining_indices, start_index)
  initial_loss <- compute_subset_loss(selected_indices)
  loss_values <- c(loss_values, initial_loss)

  # Start the selection process
  for (step in seq_len(n_factors - 1)) {  # Loop through remaining factors
    best_loss <- Inf
    best_index <- NULL

    for (idx in remaining_indices) {
      current_indices <- c(selected_indices, idx)
      current_loss <- compute_subset_loss(current_indices)

      if (current_loss < best_loss) {
        best_loss <- current_loss
        best_index <- idx
      }
    }

    # Update selected and remaining indices
    selected_indices <- c(selected_indices, best_index)
    remaining_indices <- setdiff(remaining_indices, best_index)
    loss_values <- c(loss_values, best_loss)
  }

  # Optionally plot the elbow plot
  if (plot) {
    plot(seq_along(loss_values), loss_values, type = "b", pch = 19,
         xlab = "Number of Factors Included", ylab = "Loss Function Value",
         main = "Elbow Plot for Factor Selection")
    # Add text labels for the selected indices
    text(seq_along(loss_values), loss_values, labels = selected_indices, pos = 3, cex = 0.8)
  }

  # Return results
  list(selected_order = selected_indices, loss_values = loss_values)
}


# Function to generate an elbow plot for factor selection (grouped version)
elbow_plot_selection_grouped <- function(loadings, start_index = 1, plot = TRUE,
                                         normalized = FALSE, method = "nn",
                                         groups, group_weight = FALSE,
                                         fast = FALSE, fast.num = 3000) {
  # Check inputs
  if (!is.matrix(loadings)) stop("'loadings' must be a matrix")
  if (length(groups) != nrow(loadings)) stop("'groups' must match the number of rows in 'loadings'")

  # if fast = TRUE, then for each group, random sample fast.num rows for each group (or all rows if the number of rows is less than fast.num in this group)
  if (fast) {
    new_loadings <- matrix(0, 0, ncol(loadings))
    new_groups <- NULL
    unique_groups <- unique(groups)
    for (group in unique_groups) {
      group_indices <- which(groups == group)
      if (length(group_indices) > fast.num) {
        sample_indices <- sample(group_indices, fast.num)
        new_loadings <- rbind(new_loadings, loadings[sample_indices,])
        new_groups <- c(new_groups, as.character(groups[sample_indices]))
      } else {
        new_loadings <- rbind(new_loadings, loadings[group_indices,])
        new_groups <- c(new_groups, as.character(rep(group, length(group_indices))))
      }
    }
    loadings <- new_loadings
    groups <- new_groups
  }

  n_factors <- ncol(loadings)
  remaining_indices <- seq_len(n_factors)
  selected_indices <- numeric(0)
  loss_values <- numeric(0)

  # Helper function to compute grouped loss for a given subset of loadings
  compute_subset_loss_grouped <- function(indices) {
    subset_loadings <- loadings[, indices, drop = FALSE]
    order_list <- recover_ordering_otsp_grouped(subset_loadings, groups = groups, method = method)
    compute_loss_grouped(subset_loadings, order_list, groups,
                         distance_func = function(x, y) sqrt(sum((x - y)^2)),
                         normalized = normalized, group_weight = group_weight)
  }

  # Add the start_index as the first selected factor
  selected_indices <- c(start_index)
  remaining_indices <- setdiff(remaining_indices, start_index)
  initial_loss <- compute_subset_loss_grouped(selected_indices)
  loss_values <- c(loss_values, initial_loss)

  # Start the selection process
  for (step in seq_len(n_factors - 1)) {  # Loop through remaining factors
    best_loss <- Inf
    best_index <- NULL

    for (idx in remaining_indices) {
      current_indices <- c(selected_indices, idx)
      current_loss <- compute_subset_loss_grouped(current_indices)

      if (current_loss < best_loss) {
        best_loss <- current_loss
        best_index <- idx
      }
    }

    # Update selected and remaining indices
    selected_indices <- c(selected_indices, best_index)
    remaining_indices <- setdiff(remaining_indices, best_index)
    loss_values <- c(loss_values, best_loss)
  }

  # Optionally plot the elbow plot
  if (plot) {
    plot(seq_along(loss_values), loss_values, type = "b", pch = 19,
         xlab = "Number of Factors Included", ylab = "Loss Function Value",
         main = "Elbow Plot for Factor Selection (Grouped)")
    # Add text labels for the selected indices
    text(seq_along(loss_values), loss_values, labels = selected_indices, pos = 3, cex = 0.8)
  }

  # Return results
  list(selected_order = selected_indices, loss_values = loss_values)
}

