## Functions needed for evaluation of 2d experiments

# check if x has correct format:
test_dmv = function (x, ncols = NULL, len = NULL) {
  if (is.null(ncols) & is.null(len)) {
    tdx = test_data_frame(x = x, types = "numeric",
                          any.missing = FALSE, min.rows = 1, min.cols = 1,
                          null.ok = FALSE)
    tmx = test_matrix(x = x, mode = "numeric", min.rows = 1,
                      min.cols = 1, any.missing = FALSE, null.ok = FALSE)
    tvx = test_vector(x = x, strict = TRUE, any.missing = FALSE,
                      null.ok = FALSE)
    assert_true(tdx | tmx | tvx)
    if (tvx) {
      assert_numeric(x = x, finite = TRUE)
      x = matrix(x, ncol = 1)
    }
  } else if (is.null(len)) {
    tdx = test_data_frame(x = x, types = "numeric",
                          any.missing = FALSE, min.rows = 1, ncols = ncols,
                          null.ok = FALSE)
    tmx = test_matrix(x = x, mode = "numeric", min.rows = 1,
                      ncols = ncols, any.missing = FALSE, null.ok = FALSE)
    tvx = test_vector(x = x, strict = TRUE, any.missing = FALSE,
                      null.ok = FALSE)
    assert_true(tdx | tmx | tvx)
    if (tvx) {
      assert_numeric(x = x, finite = TRUE)
      x = matrix(x, ncol = ncols)
    }
  } else {
    tdx = test_data_frame(x = x, types = "numeric",
                          any.missing = FALSE, min.rows = 1, ncols = ncols,
                          null.ok = FALSE)
    tmx = test_matrix(x = x, mode = "numeric", min.rows = 1,
                      ncols = ncols, any.missing = FALSE, null.ok = FALSE)
    tvx = test_vector(x = x, strict = TRUE, any.missing = FALSE,
                      len = len, null.ok = FALSE)
    assert_true(tdx | tmx | tvx)
    if (tvx) {
      assert_numeric(x = x, finite = TRUE)
      x = matrix(x, nrow = 1)
    }
  }
  return(x)
}

# draw random numbers from multivariate normal distribution:
rmvnorm_par = function(n, par) {
  return(mixtools::rmvnorm(n = n, mu = par$mean, sigma = par$sigma))
}

# sampling calculating E and V:
sampling_ev = function(f, x, rn) {
  if (is.character(x = f)) f = get(f)
  assert_function(x = f, args = c("x", "c"), null.ok = FALSE)
  if (test_vector(x = x, strict = TRUE, any.missing = FALSE,
                  null.ok = FALSE)) {
    x = test_dmv(x = x, ncol = 1)
  } else {
    x = test_dmv(x = x, ncols = ncol(x))
  }
  colnames(x) = paste0("x", 1:ncol(x))
  if (test_vector(x = rn, strict = TRUE, any.missing = FALSE,
                  null.ok = FALSE)) {
    rn = test_dmv(x = rn, ncol = 1)
  } else {
    rn = test_dmv(x = rn, ncols = ncol(rn))
  }
  colnames(rn) = paste0("c", 1:ncol(rn))
  ev = apply(X = x, MARGIN = 1, FUN = function(y) {
    y = matrix(y, nrow = 1)
    y = y[rep(1, nrow(rn)), ]
    fv = f(x = y, c = rn)
    e = mean(fv)
    v = var(fv)
    return(c(e, v))
  })
  return(data.frame(x, e = ev[1, ], v = ev[2, ]))
}

# transformation of data:
transform_data = function (data, maximize = TRUE, sd = TRUE, n1 = "e", n2 = "v") {
  data = as.data.table(data)
  if (maximize) {
    ind = which(colnames(data) == n1)
    data[[ind]] = -data[[ind]]
    colnames(data)[ind] = paste0("m", n1)
  }
  if (sd) {
    ind = which(colnames(data) == n2)
    data[[ind]] = sqrt(data[[ind]])
    colnames(data)[ind] = "s"
  }
  return(data)
}

# p-norm:
pnorm = function(x, y, p = 2) {
  if(is.finite(p)) {
    return(as.numeric(dist(x = rbind(x, y), method = "minkowski", p = p)))
  } else {
    return(as.numeric(dist(x = rbind(x, y), method = "maximum")))
  }
}

# distance between x and Y:
dxY = function(x, Y, p = 2) {
  return(min(apply(Y, 1, function(z) pnorm(x = x, y = z, p = p))))
}

# GD:
GDp = function(X, Y, p = 2) {
  return((1 / nrow(X) * sum(apply(X, 1, function(z) dxY(z, Y, p = p)^p)))^(1 / p))
}

# IGD:
IGDp = function(X, Y, p = 2) {
  return((1 / nrow(Y) * sum(apply(Y, 1, function(z) dxY(z, X, p = p)^p)))^(1 / p))
}

# Delta_p:
deltap = function(X, Y, p = 2) {
  return(max(GDp(X, Y, p), IGDp(X, Y, p)))
}

# quality measures:
calculate_quality_measures = function(est, true, n1 = "me", n2 = "s") {
  colnames(est)[which(colnames(est) == n1)] = "n1"
  colnames(est)[which(colnames(est) == n2)] = "n2"
  colnames(true)[which(colnames(true) == n1)] = "n1"
  colnames(true)[which(colnames(true) == n2)] = "n2"
  est = data.table::as.data.table(est)
  est = est[, list(n1, n2)]
  est = apply(est, 2, as.numeric)
  true = as.data.table(true)
  true = true[, list(n1, n2)]
  true = apply(true, 2, as.numeric)
  indtrue = ecr::nondominated(t(true))
  indest = ecr::nondominated(t(est))
  deltap2 = deltap(as.matrix(est[indest, , drop = FALSE]),
                   as.matrix(true[indtrue, , drop = FALSE]), p = 2)
  distall = mean(vapply(1:nrow(est), function(x)
    pnorm(as.numeric(est[x, ]), as.numeric(true[x, ])), numeric(1)))
  return(c(deltap2 = deltap2, distall = distall))
}

# plot quality:
plot_quality = function (quality, measure, group, method, f, title = "",
                         xlab = method, ylab = measure, rawdata = FALSE) {
  colnames(quality)[which(colnames(quality) == measure)] = "measure"
  colnames(quality)[which(colnames(quality) == group)] = "group"
  colnames(quality)[which(colnames(quality) == method)] = "method"
  colnames(quality)[which(colnames(quality) == f)] = "f"
  if (rawdata) return(quality)
  g = ggplot(data = quality, mapping = aes(x = method, y = measure, group = method)) +
    geom_boxplot() + labs(title = title, x = xlab, y = ylab) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                       legend.key = element_rect(fill = "transparent", colour = NA),
                       legend.title = element_blank(), legend.position = "right",
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  if (!is.null(f) & !is.null(group)) {
    g = g + facet_grid(group ~ f)
  }
  else if (is.null(f)) {
    g = g + facet_wrap(~group, nrow = 1)
  }
  else if (is.null(group)) {
    g = g + facet_wrap(~f, nrow = 1)
  }
  return(g)
}
