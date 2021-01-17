# implementation of RANSAC-algorithm
# (parallelization using package "doRNG": register parallel backends to use)
# input: formula: as for lm-function,
#        data: data.frame with all relevant variables
#        error_threshold: smallest deviation between observation and model
#                         for inliers
#        inlier_threshold: minimal number of observations in consensus set
#        iterations: number of subsamples
#        minimum_observations: minimal number of observations in a subsample
#        lm_arguments: additional arguments for lm-function
#        seed: seed for random numbers
# output: list of 'best' model and the input data with additional logical column
#         (".consensus_set") specifying the consensus set
#         (if satisfying model is not found, list of NULL and original data)
#         (if input data contains NAs in one of the variables in the formula,
#          omitted data is in the output)
ransaclm <- function(formula, data, error_threshold, inlier_threshold,
                     iterations = 100, minimum_observations = inlier_threshold,
                     lm_arguments = list(), seed = 123456) {
  require("doRNG")
  ### input check ##
  checkmate::assert_formula(formula)
  if (!is.data.frame(data)) {
    data <- try(as.data.frame(data))
    if (inherits(data, "try-error")) {
      stop("<data> not convertable to a data.frame")
    } else {
      message("<data> converted to data.frame.")
    }
  }
  checkmate::assert_data_frame(data, min.cols = 2, all.missing = FALSE)
  original_data <- data
  # check for NAs in data
  data <- try(model.frame(formula, original_data, na.action = na.fail), silent = TRUE)
  if (inherits(data, "try-error")) {
    data <- model.frame(formula, original_data, na.action = na.omit)
    data_omitted <- TRUE
    warning("omitted rows in <data>, where NA in one of the variables in the formula")
  }
  checkmate::assert_number(error_threshold, lower = 0, finite = TRUE)
  checkmate::assert_integerish(inlier_threshold, lower = 1, upper = nrow(data))
  checkmate::assert_integerish(iterations, lower = 1)
  checkmate::assert_integerish(minimum_observations, lower = 1, upper = nrow(data))
  checkmate::assert_list(lm_arguments)
  checkmate::assert_integerish(seed)

  number_of_fitted_coefs <- ncol(model.matrix(formula, data)) - 1
  if (minimum_observations < number_of_fitted_coefs) {
    stop("<minimum_observations> is smaller than the number of coefficients to fit")
  }
  ##################
  set.seed(seed)

  subsample_models <- foreach(seq_len(iterations),
    .multicombine = TRUE,
    .export = c(
      "fit_subsample_model",
      "find_best_consensus_set",
      "get_error", "get_error_not_subsample",
      "find_inlier", "fit_lm"
    )
  ) %dorng%
    fit_subsample_model(
      formula = formula,
      data = data,
      minimum_observations = minimum_observations,
      lm_arguments = lm_arguments,
      error_threshold = error_threshold
    )
  best_consensus_set <- find_best_consensus_set(
    subsample_models = subsample_models,
    inlier_threshold = inlier_threshold
  )
  best_model <- fit_lm(
    formula = formula,
    data = data[best_consensus_set, ],
    lm_arguments = lm_arguments
  )
  if (sum(best_consensus_set) < 0.1 * number_of_fitted_coefs) {
    warning("number of fitted coefficients is greater than 10% of the number of
            inliers in the final model")
  }
  best_data <- dplyr::mutate(data, ".consensus_set" = best_consensus_set)
  list(
    model = best_model,
    data = best_data
  )
}

# fit linear model on subsampled data and determine error vector, mean squared
# error, and inliers
# input: formula: as for lm-function,
#        data: full data,
#        minumum_observations: minimal number of observations in a subsample,
#        lm_arguments: additional arguments for lm-function
#        error_threshold: smallest deviation between observation and model
#                         for inliers
# output: list of
#           $subsampled_rows: subsampled rows,
#           $model: fitted model,
#           $errors: error vector,
#           $is_inlier: logical vector indicating inliers of model,
#           $consensus_set_size: number of inliers,
#           $mse: mean-squared-error
fit_subsample_model <- function(formula, data, minimum_observations, lm_arguments,
                                error_threshold) {
  subsampled_rows <- sample(1:nrow(data), minimum_observations)
  model <- fit_lm(
    formula = formula, data = data[subsampled_rows, ], lm_arguments = lm_arguments
  )
  errors <- get_error(
    data = data, subsampled_rows = subsampled_rows, model = model
  )
  mse <- sum(errors^2) / length(errors)
  is_inlier <- find_inlier(
    errors = errors, error_threshold = error_threshold
  )
  consensus_set_size <- sum(is_inlier)
  list(
    subsampled_rows = subsampled_rows,
    model = model,
    errors = errors,
    is_inlier = is_inlier,
    consensus_set_size = consensus_set_size,
    mse = mse
  )
}

# determine model that holds the inlier threshold and has the smallest
# mean-squared-error ('best' model)
# input: subsample_models: list of outputs of function fit_subsample_model,
#        inlier_threshold: minimum number of inliers in final model
# output: logical vector indicating whether observation is inlier of the 'best'
#         model (NULL is returned if consensus sets of all subsample models are
#         too small)
find_best_consensus_set <- function(subsample_models, inlier_threshold) {
  consensus_set_sizes <- sapply(subsample_models, "[[", "consensus_set_size")
  is_big_enough <- consensus_set_sizes > inlier_threshold

  if (sum(is_big_enough) == 0) {
    return(NULL)
  }
  position_of_smallest_mse <- which.min(sapply(
    subsample_models[is_big_enough], "[[", "mse"
  ))
  subsample_models[[position_of_smallest_mse]][["is_inlier"]]
}

# determine the deviation between observed target values and predicted ones
# input: data: full data,
#        subsampled_rows: row numbers of subsample,
#        model: subsample lm
# output: vector of deviations (error)
get_error <- function(data, subsampled_rows, model) {
  error <- vector(length = nrow(data))
  error[subsampled_rows] <- model[["residuals"]]

  not_subsampled_rows <- -subsampled_rows
  error[not_subsampled_rows] <- get_error_not_subsample(
    data = data, model = model, not_subsampled_rows = not_subsampled_rows
  )
  error
}

# determine the deviation between observed target values and predicted ones
# (when residuals are not given)
# input: data: full data,
#        model: subsample lm,
#        not_subsampled_rows: row numbers of observations not in the sample
# output: vector of deviations between observed rows and predictions (error)
get_error_not_subsample <- function(data, model, not_subsampled_rows) {
  predictions_not_subsampled <- unname(predict.lm(
    model, data[not_subsampled_rows, ]
  ))
  target <- names(model.frame(model))[1]
  target_not_subsampled <- data[not_subsampled_rows, target]
  predictions_not_subsampled - target_not_subsampled
}

# determine inliers
# input: errors: vector of deviation of observed target and prediction,
#        error_threshold: smallest deviation between observation and model
#                         for inliers
# output: logical vector whether error is smaller than threshold
find_inlier <- function(errors, error_threshold) {
  abs(errors) < error_threshold
}

# fit lm (additional lm-arguments might be given)
# input: formula: as for lm-function,
#        data: data.frame with relevant variables,
#        lm_arguments: additional arguments for lm-function
# output: fitted lm-model
fit_lm <- function(formula, data, lm_arguments) {
  lm_arguments[["formula"]] <- formula
  lm_arguments[["data"]] <- data
  do.call(lm, lm_arguments)
}

