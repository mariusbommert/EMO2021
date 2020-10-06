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

# check names of design: 
check_design_names = function(design) {
  cn = colnames(design)
  sf = substr(cn, 1, 1)
  sl = as.numeric(substring(cn, 2))
  dsl = diff(sl)
  if (length(dsl) > 1) {
    cdn = test_names(x = unique(sf), identical.to = c("x", "c")) & 
      (length(which(dsl == -1)) == 1) & all(abs(dsl) == 1)
    nx = which(dsl == -1)
  }
  else {
    cdn = (test_names(x = sf, identical.to = c("x", "c")) & (dsl == 0))
    nx = 1
  }
  return(list(cdn = cdn, nx = nx))
}

# generate design using x and c:
design_xc = function (x, c) {
  x = test_dmv(x = x)
  c = test_dmv(x = c)
  assert_true(nrow(x) == nrow(c))
  design = cbind(x, c)
  colnames(design) = c(paste0("x", 1:ncol(x)), paste0("c", 1:ncol(c)))
  return(as.data.frame(design))
}

# density for normal distribution: 
dnorm_par = function(x, par) {
  x = as.vector(test_dmv(x = x, ncols = 1))
  return(dnorm(x = x, mean = par$mean, sd = par$sd))
}

# draw random numbers from normal distribution: 
rnorm_par = function(n, par) {
  return(rnorm(n = n, mean = par$mean, sd = par$sd))
}

# draw random numbers from gaussian mixture model: 
rmix2norm_par = function(n, par) {
  ind = sample(x = 1:2, size = n, replace = TRUE, prob = par$lambdas)
  return(rnorm(n, mean = par$means[ind], sd = par$sigmas[ind]))
}

# get function for drawing random numbers:
get_random_number_par = function (distribution) {
  assert_choice(x = distribution, choices = c("norm", "mix2norm"), 
                null.ok = FALSE)
  random_number = switch(distribution, norm = rnorm_par, mix2norm = rmix2norm_par)
  return(random_number)
}

# evaluate design: 
evaluate_design = function(design, f) {
  tmp = check_design_names(design)
  assert_true(x = tmp$cdn, na.ok = FALSE)
  assert_function(x = f, args = c("x", "c"), null.ok = FALSE)
  indx = 1:tmp$nx
  indc = (tmp$nx + 1):ncol(design)
  data = data.frame(design, value = apply(design, 1, function(x) 
    f(x[indx], matrix(x[indc], nrow = 1))))
  return(data)
}

# build a kriging model: 
build_kriging_model = function(data, covtype = "matern3_2", multistart = 100, 
                               factr = 1e+07, formula = ~1, pop.size = 20) {
  indxc = which(colnames(data) != "value")
  tmp = check_design_names(design = data[, indxc])
  model = DiceKriging::km(formula = formula, design = data[, indxc], 
                          response = data[, "value"], covtype = covtype, 
                          nugget.estim = TRUE, multistart = multistart, 
                          control = list(trace = 0, maxit = 1000, factr = factr, 
                                         pop.size = pop.size))
  return(model)
}

# prediction for kriging model: 
predict_model = function(x, c, model) {
  tmp = check_design_names(model@X)
  m = ncol(model@X)
  x = test_dmv(x = x, ncols = tmp$nx)
  c = test_dmv(x = c, ncols = m - tmp$nx)
  newdata = design_xc(x = x, c = c)
  return(DiceKriging::predict(object = model, newdata = newdata, 
                              type = "UK")$mean)
}

# convert model into function: 
model_to_f = function(model) {
  f = function(x, c) predict_model(x = x, c = c, model = model)
  return(f)
}

# caclulate integration limits: 
int_lim = function (values, f = 3) {
  values = test_dmv(x = values)
  H = ks::hpi(x = values)
  lower = apply(values, 2, min) - f * H
  upper = apply(values, 2, max) + f * H
  return(list(lower = as.numeric(lower), upper = as.numeric(upper)))
}

# numerical integration: 
mvQuad_ev = function(f, x, a, b, pdf, par, type = "cNC4", level = 10) {
  x = test_dmv(x = x)
  colnames(x) = paste0("x", 1:ncol(x))
  grid = createNIGrid(dim = length(a), type = type, level = level)
  rescale(object = grid, domain = cbind(a, b))
  nodes = getNodes(grid)
  colnames(nodes) = paste0("c", 1:length(a))
  weights = getWeights(grid)
  ev = apply(X = x, MARGIN = 1, FUN = function(y) {
    y = matrix(as.numeric(y), nrow = 1)
    y = y[rep(1, nrow(nodes)), ]
    ev = f(y, nodes) * pdf(x = nodes, par = par)
    e = sum(ev * weights)
    vv = (f(y, nodes) - e)^2 * pdf(x = nodes, par = par)
    v = sum(vv * weights)
    return(c(e, v))
  })
  return(data.frame(x, e = ev[1, ], v = ev[2, ]))
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

# 1d objective function: 
f_artificial = function (x, c) {
  x = test_dmv(x = x, ncols = 1)
  c = test_dmv(x = c, ncols = 1)
  return(x[, 1] * cos(c[, 1] - x[, 1]^2))
}