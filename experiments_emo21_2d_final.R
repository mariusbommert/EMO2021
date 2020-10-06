## Code for 2d experiments

library(batchtools)
library(checkmate)

################################################################################

## calculate improved latin hypercube designs for the experiments

# set.seed(1)
# ns = c(50, 100, 200, 500, 1000, 2000, 5000, 20000)
# ilhds = vector("list", 8)
# for(i in 1:8) {
#   ilhds[[i]] = generate_lhd(n = ns[i], lower = c(-pi, -pi), upper = c(pi, pi),
#                              type = "improved")
# }
# names(ilhds) = ns
# save(ilhds, file = "ilhds.RData")

################################################################################

set.seed(1)

reg = makeExperimentRegistry(file.dir = "experiments_emo21_2d_final")

generate_data = function(data, job, n) {
  source("functions_experiments_2d.R")
  library(checkmate)
  
  ## generate design:
  load(file = "ilhds.RData")
  xv2d = ilhds[[which(names(ilhds) == n)]]
  cv2d = rmvnorm_par(n = n, par = list(mean = c(50, 70),
                                       sigma = matrix(c(9, -6, -6, 9), nrow = 2)))
  design2d = design_xc(x = xv2d, c = cv2d)

  ## evaluate design:
  data2d = evaluate_design(design = design2d, f = f_artificial_2dx2dc)

  ## build model:
  model2d = build_kriging_model(data = data2d)

  ## convert model to f:
  f2d = model_to_f(model = model2d)

  ## x-values for evaluation:
  xv2d = ilhds[[8]]

  return(list(xv2d = xv2d, cv2d = cv2d, f2d = f2d))
}

integration_2d_known = function(job, data, instance, level = 10, ncpus = 10) {
  source("functions_experiments_2d.R")
  library(checkmate)
  library(mvQuad)
  library(parallelMap)

  ## estimate parameters:
  mean2d = colMeans(instance$cv2d)
  sigma2d = cov(instance$cv2d)

  ## split x-values for different cores:
  nx = nrow(instance$xv2d)
  sx = split(instance$xv2d, rep(1:ncpus, each = ceiling(nx / ncpus))[1:nx])

  ## integration limits:
  il2dk = int_lim(values = instance$cv2d)

  ## export:
  parallelLibrary("checkmate", "mvQuad")
  parallelExport("instance", "level", "sx", "mean2d", "sigma2d", "il2dk", "test_dmv",
                 "dmvnorm_par", "rmvnorm_par", "design_xc", "evaluate_design", 
                 "f_artificial_2dx2dc", "build_kriging_model", "model_to_f")

  ## parallel mvQuad:
  i2dk = parallelLapply(sx, function(y) {
    r = apply(y, 1, function(z) {
      tmp = mvQuad_ev(f = instance$f2d, a = il2dk$lower, b = il2dk$upper,
                      x = matrix(z, nrow = 1), pdf = dmvnorm_par,
                      par = list(mean = mean2d, sigma = sigma2d), level = level,
                      type = "cNC4")
      return(tmp)
    })
    do.call(rbind, r)
  })
  i2dk = do.call(rbind, i2dk)
  return(i2dk)
}

integration_2d_unknown = function(job, data, instance, level = 10, ncpus = 10) {
  source("functions_experiments_2d.R")
  library(checkmate)
  library(mvQuad)
  library(parallelMap)

  ## estimated density:
  edk2d = edk(values = instance$cv2d)

  ## split x-values for different cores:
  nx = nrow(instance$xv2d)
  sx = split(instance$xv2d, rep(1:ncpus, each = ceiling(nx / ncpus))[1:nx])

  ## export:
  parallelLibrary("checkmate", "mvQuad")
  parallelExport("instance", "level", "sx", "edk2d", "rmvnorm_par", "design_xc", 
                 "test_dmv", "evaluate_design", "f_artificial_2dx2dc", 
                 "build_kriging_model", "model_to_f")

  ## parallel mvQuad:
  i2du = parallelLapply(sx, function(y) {
    r = apply(y, 1, function(z) {
      tmp = mvQuad_ev(f = instance$f2d, a = edk2d$par$bounds$lower,
                      b = edk2d$par$bounds$upper, x = matrix(z, nrow = 1),
                      pdf = edk2d$pdf, par = edk2d$par, level = level,
                      type = "cNC4")
      return(tmp)
    })
    do.call(rbind, r)
  })
  i2du = do.call(rbind, i2du)
  return(i2du)
}

addProblem(name = "generate_data", fun = generate_data, seed = 1)
addAlgorithm(name = "integration_2d_known", fun = integration_2d_known)
addAlgorithm(name = "integration_2d_unknown", fun = integration_2d_unknown)

pdes2d = data.frame(n = c(50, 100, 200, 500, 1000, 2000, 5000))

addExperiments(prob.designs = list(generate_data = pdes2d),
               algo.designs = list(integration_2d_known = data.frame()), repls = 100)
addExperiments(prob.designs = list(generate_data = pdes2d),
               algo.designs = list(integration_2d_unknown = data.frame()), repls = 100)

## to submit the Jobs to SLURM-Cluster:
# submitJobs(findNotSubmitted(), reso = list(memory = 6000, walltime = 48 * 3600, ncpus = 10, pm.backend = "multicore")))

## to reduce the results after all jobs are finished:
# ujt = unwrap(getJobTable())
# save(ujt, file = "ujt2d_emo21_final.RData")
# res = reduceResultsDataTable()
# save(res, file = "res2d_emo21_final.RData")