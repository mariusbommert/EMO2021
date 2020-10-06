## Code for evaluation of 1d experiments

library(batchtools)
library(checkmate)
library(data.table)
library(ggplot2)

set.seed(1)

## Some code is commented out. The final results are loaded instead. 
## The code which is commented out can be used for generating the final results.

source("functions_evaluation_1d.R")
# load("res1d_emo21_final.RData")
# load("ujt1d_emo21_final.RData")
# 
# ujt1d = unwrap(ujt1d, cols = c("algo.pars", "resources"))
# tmp = data.table(n = unlist(lapply(ujt1d$prob.pars, function(x) x[[1]])),
#                  distribution = unlist(lapply(ujt1d$prob.pars, function(x) x[[2]])),
#                  par = lapply(ujt1d$prob.pars, function(x) x[[3]]))
# ujt = cbind(ujt1d[, list(time.running, algorithm)], tmp)
# result_1d = cbind(ujt, res1d)

################################################################################

## calculate truth for all settings (only b, c and d are used for the paper)

xv1d = seq(from = -pi/2, to = pi/2, length.out = 101)
true1d_f1a = f_artificial_ev(x = xv1d, par = list(mean = 1, sd = 0.5),
                             distribution = "norm")
true1d_f1a = transform_data(true1d_f1a)
true1d_f1b = f_artificial_ev(x = xv1d, par = list(mean = 1, sd = 1),
                             distribution = "norm")
true1d_f1b = transform_data(true1d_f1b)
true1d_f1c = f_artificial_ev(x = xv1d, par = list(mean = 5, sd = 1),
                             distribution = "norm")
true1d_f1c = transform_data(true1d_f1c)

true1d_f1d = f_artificial_ev(x = xv1d, par = list(lambdas = c(0.95, 0.05), means = c(1, 10), sigmas = c(1, 1)),
                             distribution = "mix2norm")
true1d_f1d = transform_data(true1d_f1d)
true1d_f1e = f_artificial_ev(x = xv1d, par = list(lambdas = c(0.90, 0.10), means = c(1, 10), sigmas = c(1, 1)),
                             distribution = "mix2norm")
true1d_f1e = transform_data(true1d_f1e)
true1d_f1f = f_artificial_ev(x = xv1d, par = list(lambdas = c(0.80, 0.20), means = c(1, 10), sigmas = c(1, 1)),
                             distribution = "mix2norm")
true1d_f1f = transform_data(true1d_f1f)
true1d_f1g = f_artificial_ev(x = xv1d, par = list(lambdas = c(0.95, 0.05), means = c(1, 15), sigmas = c(1, 1)),
                             distribution = "mix2norm")
true1d_f1g = transform_data(true1d_f1g)

# plot truth for 1d and 2d:
f1a = plot_domination(data = true1d_f1a,
                      title = expression(paste(f[1](x, C), ", C ~ N(1, 0.25)")),
                      size_nd = 2) +
  theme(legend.position = "bottom")
legend = cowplot::get_legend(f1a)
f1a = f1a + theme(legend.position = "none")
f1b = plot_domination(data = true1d_f1b,
                      title = expression(paste(f[1](x, C), ", C ~ N(1, 1)")),
                      size_nd = 2) +
  theme(legend.position = "none")
f1c = plot_domination(data = true1d_f1c,
                      title = expression(paste(f[1](x, C), ", C ~ N(5, 1)")),
                      size_nd = 2) +
  theme(legend.position = "none")
f1d = plot_domination(data = true1d_f1d,
                      title = expression(paste(f[1](x, C), ", C ~ 0.95 x N(1, 1) + 0.05 x N(10, 1)")),
                      size_nd = 2) +
  theme(legend.position = "none")
f1e = plot_domination(data = true1d_f1e,
                      title = expression(paste(f[1](x, C), ", C ~ 0.9 x N(1, 1) + 0.1 x N(10, 1)")),
                      size_nd = 2) +
  theme(legend.position = "none")
f1f = plot_domination(data = true1d_f1f,
                      title = expression(paste(f[1](x, C), ", C ~ 0.8 x N(1, 1) + 0.2 x N(10, 1)")),
                      size_nd = 2) +
  theme(legend.position = "none")
f1g = plot_domination(data = true1d_f1g,
                      title = expression(paste(f[1](x, C), ", C ~ 0.95 x N(1, 1) + 0.05 x N(15, 1)")),
                      size_nd = 2) +
  theme(legend.position = "none")

################################################################################

## calculate quality measures

# qual = matrix(nrow = 8400, ncol = 2)
# inda = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(mean = 1, sd = 0.5))), is.logical)))
# indb = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(mean = 1, sd = 1))), is.logical)))
# indc = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(mean = 5, sd = 1))), is.logical)))
# indd = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(lambdas = c(0.95, 0.05), means = c(1, 10), sigmas = c(1, 1)))), is.logical)))
# inde = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(lambdas = c(0.90, 0.10), means = c(1, 10), sigmas = c(1, 1)))), is.logical)))
# indf = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(lambdas = c(0.80, 0.20), means = c(1, 10), sigmas = c(1, 1)))), is.logical)))
# indg = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(lambdas = c(0.95, 0.05), means = c(1, 15), sigmas = c(1, 1)))), is.logical)))
# for(i in 1:1200)  {
#   qual[inda[i], ] = calculate_quality_measures(est = transform_data(result_1d$result[[inda[i]]][, 2:3]),
#                                                true = true1d_f1a[, 2:3], n1 = "me", n2 = "s")
#   qual[indb[i], ] = calculate_quality_measures(est = transform_data(result_1d$result[[indb[i]]][, 2:3]),
#                                                true = true1d_f1b[, 2:3], n1 = "me", n2 = "s")
#   qual[indc[i], ] = calculate_quality_measures(est = transform_data(result_1d$result[[indc[i]]][, 2:3]),
#                                                true = true1d_f1c[, 2:3], n1 = "me", n2 = "s")
#   qual[indd[i], ] = calculate_quality_measures(est = transform_data(result_1d$result[[indd[i]]][, 2:3]),
#                                                true = true1d_f1d[, 2:3], n1 = "me", n2 = "s")
#   qual[inde[i], ] = calculate_quality_measures(est = transform_data(result_1d$result[[inde[i]]][, 2:3]),
#                                                true = true1d_f1e[, 2:3], n1 = "me", n2 = "s")
#   qual[indf[i], ] = calculate_quality_measures(est = transform_data(result_1d$result[[indf[i]]][, 2:3]),
#                                                true = true1d_f1f[, 2:3], n1 = "me", n2 = "s")
#   qual[indg[i], ] = calculate_quality_measures(est = transform_data(result_1d$result[[indg[i]]][, 2:3]),
#                                                true = true1d_f1g[, 2:3], n1 = "me", n2 = "s")
# }
# colnames(qual) = c("deltap2", "distall")
# result_1d = cbind(result_1d, qual)
# 
# save(result_1d, file = "result_1d_emo21_final.RData")

################################################################################

load("result_1d_emo21_final.RData")

## rename some levels of variables

result_1d$n = factor(result_1d$n, levels = c(10, 20, 50, 100, 200, 500))

result_1d$algorithm[which(result_1d$algorithm == "integration_1d_known")] = "k"
result_1d$algorithm[which(result_1d$algorithm == "integration_1d_unknown")] = "u"

inda = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(mean = 1, sd = 0.5))), is.logical)))
indb = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(mean = 1, sd = 1))), is.logical)))
indc = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(mean = 5, sd = 1))), is.logical)))
indd = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(lambdas = c(0.95, 0.05), means = c(1, 10), sigmas = c(1, 1)))), is.logical)))
inde = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(lambdas = c(0.90, 0.10), means = c(1, 10), sigmas = c(1, 1)))), is.logical)))
indf = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(lambdas = c(0.80, 0.20), means = c(1, 10), sigmas = c(1, 1)))), is.logical)))
indg = which(unlist(lapply(lapply(result_1d$par, function(x) all.equal(x, list(lambdas = c(0.95, 0.05), means = c(1, 15), sigmas = c(1, 1)))), is.logical)))

result_1d[, ms := algorithm]
result_1d$ms[inda] = "a"
result_1d$ms[indb] = "b"
result_1d$ms[indc] = "c"
result_1d$ms[indd] = "d"
result_1d$ms[inde] = "e"
result_1d$ms[indf] = "f"
result_1d$ms[indg] = "g"

################################################################################
################################################################################

## Code for Figure 1:

load("true2d_emo21_final.RData")

lab = function(labels) {
  list(list(expression(paste(f[1], " (a)")),
            expression(paste(f[1], " (b)")),
            expression(paste(f[1], " (c)"))))
}
f1b = plot_domination(data = true1d_f1b, rawdata = TRUE)
f1c = plot_domination(data = true1d_f1c, rawdata = TRUE)
f1d = plot_domination(data = true1d_f1d, rawdata = TRUE)
f2 = plot_domination(data = true2d, rawdata = TRUE)
data1d = rbind(cbind(f1c[, 2:4], f = "a"),
               cbind(f1b[, 2:4], f = "b"),
               cbind(f1d[, 2:4], f = "c"))
g1 = ggplot(data = data1d[data1d$dom == "non dominated", ],
            mapping = aes(x = n1, y = n2, col = dom, size = dom, shape = dom)) +
  geom_point() + geom_point(data = data1d[data1d$dom == "dominated", ],
                            mapping = aes(x = n1, y = n2, col = dom, size = dom)) +
  labs(title = "", x = "-E(f(x, C))", y = "S(f(x, C))", size = data1d$dom) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     legend.key = element_rect(fill = "transparent", colour = NA),
                     legend.title = element_blank(), legend.position = "right",
                     legend.text = element_text(size = 16),
                     axis.text = element_text(size = 12),
                     axis.title = element_text(size = 16),
                     strip.text = element_text(size = 14)) +
  scale_color_manual("", values = c("red", "black"), labels = c("non dominated", "dominated"),
                     breaks = c("non dominated", "dominated")) +
  scale_size_manual("", values = c(2, 1.5), labels = c("non dominated", "dominated"),
                    breaks = c("non dominated", "dominated")) +
  scale_shape_manual("", values = c(17, 16), labels = c("non dominated", "dominated"),
                     breaks = c("non dominated", "dominated")) +
  theme(legend.position = "bottom")
legend = cowplot::get_legend(g1)
p1 = g1 + facet_wrap(~f, labeller = lab) + theme(legend.position = "none")

lab = function(labels) {
  list(list(expression(paste(f[2]))))
}
data2d = cbind(f2[, 3:5], f = "f2(x, C)")
g2 = ggplot(data = data2d[data2d$dom == "non dominated", ],
            mapping = aes(x = n1, y = n2, col = dom, size = dom, shape = dom)) +
  geom_point() + geom_point(data = data2d[data2d$dom == "dominated", ],
                            mapping = aes(x = n1, y = n2, col = dom, size = dom)) +
  labs(title = "", x = "-E(f(x, C))", y = "S(f(x, C))", size = data2d$dom) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     legend.key = element_rect(fill = "transparent", colour = NA),
                     legend.title = element_blank(), legend.position = "right",
                     legend.text = element_text(size = 16),
                     axis.text = element_text(size = 12),
                     axis.title = element_text(size = 16),
                     strip.text = element_text(size = 14)) +
  scale_color_manual("", values = c("red", "black"), labels = c("non dominated", "dominated"),
                     breaks = c("non dominated", "dominated")) +
  scale_size_manual("", values = c(2, 1.5), labels = c("non dominated", "dominated"),
                    breaks = c("non dominated", "dominated")) +
  scale_shape_manual("", values = c(17, 16), labels = c("non dominated", "dominated"),
                     breaks = c("non dominated", "dominated"))
p2 = g2 + facet_wrap(~f, labeller = lab) + theme(legend.position = "none")
p = cowplot::plot_grid(p1, p2, rel_widths = c(3/4, 1/4))

pdf("true.pdf", width = 12, height = 5)
cowplot::plot_grid(p, legend, nrow = 2, rel_heights = c(0.9, 0.1))
dev.off()

################################################################################

## Code for Figure 2a

labr = function(labels) {
  list('paste(Delta[p])',
       'paste("MED")')
}

labc = function(labels) {
  list('paste(f[1], " (a), MIOS-PDE")',
       'paste(f[1], " (a), MIOS-KDE")',
       'paste(f[1], " (b), MIOS-PDE")',
       'paste(f[1], " (b), MIOS-KDE")',
       'paste(f[1], " (c), MIOS-PDE")',
       'paste(f[1], " (c), MIOS-KDE")')
}

q1d = result_1d
q1d$ms[indc] = "b"
q1d$ms[indb] = "c"
q1d$ms[indd] = "d"

q1ddp = plot_quality(quality = q1d[c(indb, indd, indc), list(deltap2, n, algorithm, ms)],
                     measure = "deltap2", group = "algorithm", method = "n",
                     f = "ms", rawdata = TRUE)
q1ddp = cbind(q1ddp, m = "Delta[p]")
q1dda = plot_quality(quality = q1d[c(indb, indd, indc), list(distall, n, algorithm, ms)],
                     measure = "distall", group = "algorithm", method = "n",
                     f = "ms", rawdata = TRUE)
q1dda = cbind(q1dda, m = "MED")
qual = rbind(q1ddp, q1dda)
qual[, q := paste(f, group)]

pdf("qual_1d.pdf", width = 12, height = 6) 
ggplot(data = qual, mapping = aes(x = method, y = measure, group = method)) +
  geom_boxplot() + labs(title = "", x = "", y = "") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.title = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)) +
  facet_grid(m ~ q, labeller = labeller(m = as_labeller(labr, label_parsed),
                                        q = as_labeller(labc, label_parsed))) +
  labs(x = "Number of observations for building the metamodel (n)")
dev.off()

################################################################################

## Code for Figure 4a

q1d = result_1d
q1d[, ddistall := distall]
q1d[, ddeltap2 := deltap2]
q1d$ddistall[which(q1d$algorithm == "k")] = q1d$distall[which(q1d$algorithm == "k")] -
  q1d$distall[which(q1d$algorithm == "u")]
q1d$ddeltap2[which(q1d$algorithm == "k")] = q1d$ddeltap2[which(q1d$algorithm == "k")] -
  q1d$ddeltap2[which(q1d$algorithm == "u")]
q1d = q1d[1:4200, ]
q1d$ms[indc[1:600]] = "b"
q1d$ms[indb[1:600]] = "c"
q1d$ms[indd[1:600]] = "d"

d1 = plot_quality(quality = q1d[c(indc[1:600], indb[1:600], indd[1:600]),
                                list(ddeltap2, n, algorithm, ms)], measure = "ddeltap2",
                  group = "ms", method = "n", f = NULL, rawdata = TRUE)
d1 = cbind(d1, m = "Difference in Delta[p], (known - unknown)")
d2 = plot_quality(quality = q1d[c(indc[1:600], indb[1:600], indd[1:600]),
                                list(ddistall, n, algorithm, ms)], measure = "ddistall",
                  group = "ms", method = "n", f = NULL, rawdata = TRUE)
d2 = cbind(d2, m = "Difference in MED (known - unknown)")
d = rbind(d1, d2)

labr = function(labels) {
  list('paste("Difference in ", Delta[p], " (MIOS-PDE - MIOS-KDE)")',
       'paste("Difference in MED (MIOS-PDE - MIOS-KDE)")')
}

labc = function(labels) {
  list('paste(f[1], " (a)")',
       'paste(f[1], " (b)")',
       'paste(f[1], " (c)")')
}

pdf("diff_qual_1d.pdf", width = 12, height = 9) 
ggplot(data = d, mapping = aes(x = method, y = measure, group = method)) +
  geom_boxplot() + labs(title = "", x = "", y = NULL) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.title = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 13)) +
  facet_grid(m ~ group, labeller = labeller(m = as_labeller(labr, label_parsed),
                                            group = as_labeller(labc, label_parsed))) +
  labs(x = "Number of observations for building the metamodel (n)") +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed")
dev.off()

################################################################################
################################################################################

## Code for Figure 3

set.seed(2)

data_dens = function(xv1d, cv1d, i) {
  dens = edk(values = cv1d)
  tmp = dens$pdf(xv1d, par = dens$par)
  return(data.frame(x = xv1d,
                    densknown = dnorm(x = xv1d, mean = mean(cv1d), sd = sd(cv1d)),
                    densunknown = tmp, group = i))
}

xv1d = seq(-5, 15, length.out = 200)
cv1d = rmix2norm_par(n = 1e5, par = list(lambdas = c(0.95, 0.05), means = c(1, 10), sigmas = c(1, 1)))
dens = edk(values = cv1d)
tdens = dens$pdf(x = seq(-5, 15, length.out = 200), par = dens$par)

values_dens_mix = function(n, xv1d, prop = 0.95, reps = 100) {
  res = vector("list", 100)
  for(i in 1:reps) {
    cv1d = rmix2norm_par(n = n, par = list(lambdas = c(prop, 1 - prop), means = c(1, 10), sigmas = c(1, 1)))
    res[[i]] = data_dens(xv1d = xv1d, cv1d = cv1d, i = i)
  }
  return(res)
}

dens_med = function(n, xv1d, tdens, prop = 0.95, reps = 100) {
  tmp = values_dens_mix(n = n, xv1d = xv1d, prop = prop, reps = reps)
  original = sapply(1:reps, function(x) mean(abs(tdens - tmp[[x]][, 2])))
  extended = sapply(1:reps, function(x) mean(abs(tdens - tmp[[x]][, 3])))
  return(data.frame(med = c(original, extended),
                    approach = c(rep("MIOS-PDE", reps), rep("MIOS-KDE", reps)), n = n))
}

med = rbind(dens_med(n = 10, xv1d = xv1d, tdens = tdens),
            dens_med(n = 20, xv1d = xv1d, tdens = tdens),
            dens_med(n = 50, xv1d = xv1d, tdens = tdens),
            dens_med(n = 100, xv1d = xv1d, tdens = tdens),
            dens_med(n = 200, xv1d = xv1d, tdens = tdens),
            dens_med(n = 500, xv1d = xv1d, tdens = tdens))
med$n = as.factor(med$n)
med$approach = factor(med$approach, levels = c("MIOS-PDE", "MIOS-KDE"))

pdf("diff_density_estimation.pdf", width = 12, height = 4)
ggplot(data = med, aes(x = n, y = med)) + geom_boxplot() + facet_wrap(~ approach, nrow = 1) +
  labs(title = "", x = "Number of observations for density estimation (n)",
       y = "MED") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.title = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))
dev.off()
