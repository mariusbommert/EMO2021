## Code for 1d experiments

library(batchtools)
library(checkmate)
library(data.table)

set.seed(1)

source("functions_experiments_1d.R")

reg = makeExperimentRegistry(file.dir = "experiments_emo21_1d_final")

generate_data = function(data, job, n, distribution, par) {
  ## generate design:
  xv1d = seq(from = -pi/2, to = pi/2, length.out = n)
  rn_par = get_random_number_par(distribution = distribution)
  cv1d = rn_par(n = n, par = par)
  design1d = design_xc(x = xv1d, c = cv1d)
  
  ## evaluate design:
  data1d = evaluate_design(design = design1d, f = f_artificial)
  
  ## build model:
  model1d = build_kriging_model(data = data1d)
  
  ## convert model to f:
  f1d = model_to_f(model = model1d)
  
  ## x-values for evaluation:
  xv1d = seq(from = -pi/2, to = pi/2, length.out = 101)
  
  return(list(xv1d = xv1d, cv1d = cv1d, f1d = f1d))
}

integration_1d_known = function(job, data, instance, level = 10) {
  library(mvQuad)
  
  ## estimate parameters:
  mean1d = mean(instance$cv1d)
  sd1d = sd(instance$cv1d)
  
  ## mvQuad:
  il1dk = int_lim(values = instance$cv1d)
  i1dk = mvQuad_ev(f = instance$f1d, a = il1dk$lower, b = il1dk$upper,
                   x = instance$xv1d, pdf = dnorm_par,
                   par = list(mean = mean1d, sd = sd1d), level = level,
                   type = "cNC4")
  return(i1dk)
}

integration_1d_unknown = function(job, data, instance, level = 10) {
  library(mvQuad)
  
  ## estimated density:
  edk1d = edk(values = instance$cv1d)
  
  ## mvQuad:
  i1du = mvQuad_ev(f = instance$f1d, a = edk1d$par$bounds$lower,
                   b = edk1d$par$bounds$upper, x = instance$xv1d, pdf = edk1d$pdf,
                   par = edk1d$par, level = level, type = "cNC4")
  return(i1du)
}

addProblem(name = "generate_data", fun = generate_data, seed = 1)
addAlgorithm(name = "integration_1d_known", fun = integration_1d_known)
addAlgorithm(name = "integration_1d_unknown", fun = integration_1d_unknown)

pd1 = data.table(n = c(10, 20, 50, 100, 200, 500))
pd2 = data.table(distribution = c("norm", "norm", "norm",
                                  "mix2norm", "mix2norm", "mix2norm", "mix2norm"),
                 par = list(list(mean = 1, sd = 0.5),
                            list(mean = 1, sd = 1),
                            list(mean = 5, sd = 1),
                            list(lambdas = c(0.95, 0.05), means = c(1, 10), sigmas = c(1, 1)),
                            list(lambdas = c(0.90, 0.10), means = c(1, 10), sigmas = c(1, 1)),
                            list(lambdas = c(0.80, 0.20), means = c(1, 10), sigmas = c(1, 1)),
                            list(lambdas = c(0.95, 0.05), means = c(1, 15), sigmas = c(1, 1))))
setkeyv(pd1[, k := 1], c(key(pd1), "k"))
setkeyv(pd2[, k := 1], c(key(pd2), "k"))
pdes1d = merge(pd1, pd2, by = "k", allow.cartesian = TRUE)
pdes1d[, k := NULL]

addExperiments(prob.designs = list(generate_data = pdes1d),
               algo.designs = list(integration_1d_known = data.frame()), repls = 100)
addExperiments(prob.designs = list(generate_data = pdes1d),
               algo.designs = list(integration_1d_unknown = data.frame()), repls = 100)

## to submit the Jobs to SLURM-Cluster:
# submitJobs(findNotSubmitted())

## to reduce the results after all jobs are finished:
# ujt1d = getJobTable()
# save(ujt1d, file = "ujt1d_emo21_final.RData")
# res1d = reduceResultsDataTable()
# save(res1d, file = "res1d_emo21_final.RData")
