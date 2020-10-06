## Functions needed for evaluation of 1d experiments

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
  }
  else if (is.null(len)) {
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
  }
  else {
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

# draw random numbers from normal distribution:
rnorm_par = function(n, par) {
  return(rnorm(n = n, mean = par$mean, sd = par$sd))
}

# draw random numbers from Gaussian mixture model:
rmix2norm_par = function(n, par) {
  ind = sample(x = 1:2, size = n, replace = TRUE, prob = par$lambdas)
  return(rnorm(n, mean = par$means[ind], sd = par$sigmas[ind]))
}

# get function for random number generation:
get_random_number_par = function (distribution) {
  assert_choice(x = distribution, choices = c("norm", "mix2norm"),
                null.ok = FALSE)
  random_number = switch(distribution, norm = rnorm_par, mix2norm = rmix2norm_par)
  return(random_number)
}

# 1d objective function:
f_artificial = function (x, c) {
  x = test_dmv(x = x, ncols = 1)
  c = test_dmv(x = c, ncols = 1)
  return(x[, 1] * cos(c[, 1] - x[, 1]^2))
}

# expectation and variance for 1d objective function:
f_artificial_ev = function (x, par, distribution) {
  assert_numeric(x = x, finite = TRUE, min.len = 1, any.missing = FALSE,
                 null.ok = FALSE)
  assert_list(par, min.len = 1, any.missing = FALSE, null.ok = FALSE)
  assert_choice(x = distribution, choices = c("norm", "mix2norm"),
                null.ok = FALSE)
  rnd = get_random_number_par(distribution = distribution)
  rn = rnd(n = 1e+05, par = par)
  e = vapply(X = x, FUN = function(y) mean(f_artificial(y, rn)),
             FUN.VALUE = numeric(1))
  v = vapply(X = x, FUN = function(y) var(f_artificial(y, rn)),
             FUN.VALUE = numeric(1))
  return(data.frame(x = x, e = e, v = v))
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

# plot domination:
plot_domination = function (data, dom = NULL, n1 = "me", n2 = "s", maximize = TRUE,
                            sd = TRUE, title = "", color = c("red", "black"),
                            size_nd = 1, xlab = NULL, ylab = NULL, rawdata = FALSE) {
  colnames(data)[which(colnames(data) == n1)] = "n1"
  colnames(data)[which(colnames(data) == n2)] = "n2"
  if (is.null(dom))
    dom = ecr::nondominated(t(data[, c("n1", "n2")]))
  data = cbind(data, dom = as.factor(ifelse(dom, "non dominated",
                                            "dominated")))
  if (rawdata)
    return(data)
  if (is.null(xlab)) {
    xlab = ifelse(maximize, "-E(f(x, C))", "E(f(x, C))")
  }
  if (is.null(ylab)) {
    ylab = ifelse(sd, "S(f(x, C))", "V(f(x, C))")
  }
  g = ggplot(data = data[data$dom == "non dominated",
                         ], mapping = aes(x = n1, y = n2, col = dom, size = dom)) +
    geom_point() + geom_point(data = data[data$dom == "dominated",
                                          ], mapping = aes(x = n1, y = n2, col = dom, size = dom)) +
    labs(title = title, x = xlab, y = ylab, size = dom) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                       legend.key = element_rect(fill = "transparent", colour = NA),
                       legend.title = element_blank(), legend.position = "right") +
    scale_color_manual("", values = color, labels = c("non dominated", "dominated"),
                       breaks = c("non dominated", "dominated")) +
    scale_size_manual("", values = c(size_nd, 1), labels = c("non dominated", "dominated"),
                      breaks = c("non dominated", "dominated"))
  return(g)
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

# integration limits:
int_lim = function (values, f = 3) {
  values = test_dmv(x = values)
  H = ks::hpi(x = values)
  lower = apply(values, 2, min) - f * H
  upper = apply(values, 2, max) + f * H
  return(list(lower = as.numeric(lower), upper = as.numeric(upper)))
}

# kernel density estimation:
edk = function (values, ...) {
  tm = test_matrix(x = values, min.rows = 1, min.cols = 1,
                   max.cols = 6, any.missing = FALSE, null.ok = FALSE)
  td = test_data_frame(x = values, min.rows = 1, min.cols = 1,
                       max.cols = 6, any.missing = FALSE, null.ok = FALSE)
  tv = test_vector(x = values, strict = TRUE, any.missing = FALSE,
                   null.ok = FALSE)
  assert_true(td | tm | tv)
  if (tv) values = matrix(values, ncol = 1)
  bounds = int_lim(values = values)
  par = list(bounds = bounds, data = values)
  e = function(x, par, ...) {
    dens = ks::kde(x = par$data, eval.points = x, ...)$estimate
    return(dens)
  }
  return(list(pdf = e, par = par))
}
