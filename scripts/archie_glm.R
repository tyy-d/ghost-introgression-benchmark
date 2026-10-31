# ArchIE's published unpenalized binomial GLM with sample-size-dependent columns.
# Standardization improves numerical conditioning without changing the model family.
args <- commandArgs(trailingOnly=TRUE)
mode <- args[1]
data <- read.delim(args[2], check.names=FALSE)
if (mode == "train") {
    features <- setdiff(names(data), "label")
    means <- vapply(data[features], mean, numeric(1))
    scales <- vapply(data[features], sd, numeric(1))
    features <- features[is.finite(scales) & scales > 0]
    means <- means[features]
    scales <- scales[features]
    normalized <- as.data.frame(scale(data[features], center=means, scale=scales))
    normalized$label <- data$label
    if (length(unique(normalized$label)) != 2) stop("Training requires both ancestry classes")
    model <- glm(label ~ ., family=binomial(link="logit"), data=normalized,
                 control=glm.control(maxit=100))
    if (!model$converged) stop("ArchIE GLM did not converge")
    saveRDS(list(model=model, features=features, means=means, scales=scales), args[3])
} else if (mode == "predict") {
    fitted <- readRDS(args[3])
    normalized <- as.data.frame(scale(data[fitted$features],
                                      center=fitted$means, scale=fitted$scales))
    probabilities <- predict(fitted$model, newdata=normalized, type="response")
    if (any(!is.finite(probabilities))) stop("Nonfinite ArchIE probabilities")
    write.table(data.frame(probability=probabilities), args[4],
                sep="\t", quote=FALSE, row.names=FALSE)
} else stop("Expected train or predict")
