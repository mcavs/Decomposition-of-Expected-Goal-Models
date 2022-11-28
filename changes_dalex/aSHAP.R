library(iBreakDown)
library(ggplot2)

shap_aggregated <- function(explainer, new_observations, order = NULL, B = 25, ...) {
  ret_raw <- data.frame(contribution = c(), variable_name = c(), label = c())

  new_observations <- new_observations[,colnames(explainer$data)]

  for(i in 1:nrow(new_observations)){
    new_obs <- new_observations[i,]
    shap_vals <- iBreakDown::shap(explainer, new_observation = new_obs, B = B, ...)
    shap_vals <- shap_vals[shap_vals$B != 0, c('contribution', 'variable_name', 'label')]
    ret_raw <- rbind(ret_raw, shap_vals)
  }
  
  data_preds <- predict(explainer, explainer$data)
  mean_prediction <- mean(data_preds)
  
  subset_preds <- predict(explainer, new_observations)
  mean_subset <- mean(subset_preds)
  
  if(is.null(order)) {
    order <- calculate_order(explainer, mean_prediction, new_observations, predict)
  }
  
  ret <- raw_to_aggregated(ret_raw, mean_prediction, mean_subset, order, explainer$label)

  predictions_new <- data.frame(contribution = subset_preds, variable_name='prediction', label=ret$label[1])
  predictions_old <- data.frame(contribution = data_preds, variable_name='intercept', label=ret$label[1])
  ret_raw <- rbind(ret_raw, predictions_new, predictions_old)

  out <- list(aggregated = ret, raw = ret_raw)
  class(out) <- c('shap_aggregated', class(out))

  out
}

raw_to_aggregated <- function(ret_raw, mean_prediction, mean_subset, order, label){
  ret <- aggregate(ret_raw$contribution, list(ret_raw$variable_name, ret_raw$label), FUN=mean)
  colnames(ret) <- c('variable', 'label', 'contribution')
  ret$variable <- as.character(ret$variable)
  rownames(ret) <- ret$variable
  
  ret <- ret[order,]
  
  ret$position <- (nrow(ret) + 1):2
  ret$sign <- ifelse(ret$contribution >= 0, "1", "-1")
  
  ret <- rbind(ret, data.frame(variable = "intercept",
                               label = label,
                               contribution = mean_prediction,
                               position = max(ret$position) + 1,
                               sign = "X"),
               make.row.names=FALSE)
  
  ret <- rbind(ret, data.frame(variable = "prediction",
                               label = label,
                               contribution = mean_subset,
                               position = 1,
                               sign = "X"),
               make.row.names=FALSE)
  
  ret <- ret[call_order_func(ret$position, decreasing = TRUE), ]
  
  ret$cumulative <- cumsum(ret$contribution)
  ret$cumulative[nrow(ret)] <- ret$contribution[nrow(ret)]
  
  # some numerical errors cause little translation of the last bars
  ret$cumulative[nrow(ret) - 1] <- ret$cumulative[nrow(ret)]
  
  ret$variable_name <- ret$variable
  ret$variable_name <- factor(ret$variable_name, levels=c(ret$variable_name, ''))
  ret$variable_name[nrow(ret)] <- ''
  
  ret$variable_value <- '' # column for consistency
  
  ret
}

predict_parts_shap_aggregated <- function(explainer, new_observations, ...) {
  #test_explainer(explainer, has_data = TRUE, function_name = "predict_parts_shap_aggregated")

  res <- shap_aggregated(explainer,
                         new_observations = new_observations,
                         ...)

  class(res) <- c('predict_parts', class(res))

  res
}

call_order_func <- function(...) {
  order(...)
}

calculate_1d_changes <- function(model, new_observation, data, predict_function) {
  average_yhats <- list()
  j <- 1
  for (i in colnames(new_observation)) {
    current_data <- data
    current_data[,i] <- new_observation[,i]
    yhats <- predict_function(model, current_data)
    average_yhats[[j]] <- colMeans(as.data.frame(yhats))
    j <- j + 1
  }
  names(average_yhats) <- colnames(new_observation)
  average_yhats
}

generate_average_observation <- function(subset) {
  # takes average
  numeric_cols <- unlist(lapply(subset, is.numeric))
  numeric_cols <- names(numeric_cols[numeric_cols == TRUE])
  df_numeric <- t(as.data.frame(colMeans(subset[,numeric_cols])))

  # takes most frequent one
  factor_cols <- unlist(lapply(subset, is.factor))
  factor_cols <- names(factor_cols[factor_cols == TRUE])
  df_factory <- as.data.frame(lapply(factor_cols, function(col) {
    factor(names(which.max(table(subset[,col]))), levels = levels(subset[,col]))
  }))
  colnames(df_factory) <- factor_cols

  # also takes most frequent one
  other_cols <- unlist(lapply(subset, function(x){(!is.numeric(x)) & (!is.factor(x))}))
  other_cols <- names(other_cols[other_cols == TRUE])
  df_others <- as.data.frame(lapply(other_cols, function(col) {
    tab <- table(subset[,col])
    tab_names <- attr(tab, 'dimnames')[[1]]
    tab_names[which.max(tab)]
  }), stringsAsFactors = FALSE)
  colnames(df_others) <- other_cols

  outs <- list()
  if(!ncol(df_numeric) == 0){outs <- append(list(df_numeric), outs)}
  if(!ncol(df_factory) == 0){outs <- append(list(df_factory), outs)}
  if(!ncol(df_others) == 0){outs <- append(list(df_others), outs)}

  do.call("cbind", outs)[,colnames(subset)]
}

calculate_order <- function(x, mean_prediction, new_data, predict_function) {
  baseline_yhat <- mean_prediction

  generated_obs <- generate_average_observation(new_data)

  average_yhats <- calculate_1d_changes(x, generated_obs, x$data, predict_function)
  diffs_1d <- sapply(seq_along(average_yhats), function(i) {
    sqrt(mean((average_yhats[[i]] - baseline_yhat)^2))
  })

  order(diffs_1d, decreasing = TRUE)
}

plot.shap_aggregated <- function(x, ..., shift_contributions = 0.05, add_contributions = TRUE, add_boxplots = TRUE, max_features = 10, title = "Aggregated SHAP") {
  x <- select_only_k_features(x, k = max_features, ...)
  aggregate <- x[[1]]
  raw <- x[[2]]
  class(aggregate) <- c('break_down', class(aggregate))

  # ret has at least 3 columns: first and last are intercept and prediction
  aggregate$mean_boxplot <- c(0, aggregate$cumulative[1:(nrow(aggregate)-2)], 0)
  raw <- merge(x = as.data.frame(aggregate[,c('variable', 'position', 'mean_boxplot')]), y = raw, by.x = "variable", by.y = "variable_name", all.y = TRUE)

  # max_features = max_features + 1 because we have one more class already - "+ all other features"
  p <- plot(aggregate, ..., add_contributions = FALSE, max_features = max_features + 1, title = title)

  if(add_boxplots){
    p <- p + geom_boxplot(data = raw,
                          aes(y = contribution + mean_boxplot,
                              x = position + 0.5,
                              group = position,
                              fill = "#371ea3",
                              xmin = min(contribution) - 0.85,
                              xmax = max(contribution) + 0.85),
                          color = "#371ea3",
                          fill = "#371ea3",
                          width = 0.15)
  }

  if (add_contributions) {
    aggregate$right_side <- pmax(aggregate$cumulative,  aggregate$cumulative - aggregate$contribution)
    drange <- diff(range(aggregate$cumulative))

    p <- p + geom_text(aes(y = right_side),
                       vjust = -1,
                       nudge_y = drange*shift_contributions,
                       hjust = -0.2,
                       color = "#371ea3")
  }

  p
}

select_only_k_features <- function(input, k = 10, use_default_filter = TRUE, ...) {
  x <- input[[1]]
  y <- input[[2]]
  
  if(use_default_filter){
    # filter-out redundant rows
    contribution_sum <- tapply(x$contribution, x$variable_name, function(contribution) sum(abs(contribution), na.rm = TRUE))
    contribution_ordered_vars <- names(sort(contribution_sum[!(names(contribution_sum) %in% c("", "intercept"))]))
    variables_keep <- tail(contribution_ordered_vars, k)
    variables_remove <- setdiff(contribution_ordered_vars, variables_keep)
  } else {
    vars <- x$variable_name[!(x$variable_name %in% c("", "intercept", "prediction"))]
    variables_keep <- head(vars, k)
    variables_remove <- setdiff(vars, variables_keep)
  }

  if (length(variables_remove) > 0) {
    x_remove   <- x[x$variable_name %in% variables_remove,]
    x_keep     <- x[!(x$variable_name %in% c(variables_remove, "")),]
    x_prediction <- x[x$variable == "prediction",]
    row.names(x_prediction) <- x_prediction$label
    remainings <- tapply(x_remove$contribution, x_remove$label, sum, na.rm=TRUE)
    # fix position and cumulative in x_keep
    x_keep$position <- as.numeric(as.factor(x_keep$position)) + 2
    for (i in 1:nrow(x_keep)) {
      if (x_keep[i,"variable_name"] == "intercept") {
        x_keep[i,"cumulative"] <- x_keep[i,"contribution"]
      } else {
        x_keep[i,"cumulative"] <- x_keep[i - 1,"cumulative"] + x_keep[i,"contribution"]
      }
    }
    # for each model we shall calculate the others statistic
    x_others <- data.frame(variable = "+ all other factors",
                           contribution = remainings,
                           variable_name = "+ all other factors",
                           variable_value = "",
                           cumulative = x_prediction[names(remainings),"cumulative"],
                           sign = sign(remainings),
                           position = 2,
                           label = names(remainings))
    #
    x <- rbind(x_keep, x_others, x_prediction)
    y$variable_name <- factor(ifelse(y$variable_name %in% variables_remove, "+ all other factors", as.character(y$variable_name)), levels = levels(x$variable_name))
  }

  list(aggregated = x, raw = y)
}

predict_parts_shap <- function(explainer, new_observation, ...) {
  # run checks against the explainer objects
  test_explainer(explainer, has_data = TRUE, function_name = "predict_parts_shap")

  # call the shap from iBreakDown
  res <- iBreakDown::shap(explainer,
                          new_observation = new_observation,
                          ...)
  class(res) <- c('predict_parts', class(res))
  res
}
