## Code for evaluation of 2d experiments

library(checkmate)
library(data.table)
library(ggplot2)

set.seed(1)

## Some code is commented out. The final results are loaded instead. 
## The code which is commented out can be used for generating the final results.

source("functions_evaluation_2d.R")
# load("ilhds.RData")
# load("ujt2d_emo21_final.RData")
# load("res2d_emo21_final.RData")
# load("true2d_emo21_final.RData")

## Calculate truth

# truth:
# xv2d = ilhds[[8]]
# cv2d = rmvnorm_par(n = 1e5, par = list(mean = c(50, 70),
#                                        sigma = matrix(c(9, -6, -6, 9), nrow = 2)))
# true2d = sampling_ev(f = f_artificial_2dx2dc, x = xv2d, rn = cv2d)
# true2d = transform_data(true2d)
# save(true2d, file = "true2d_emo21_final.RData")
#
# result:
# result_2d = cbind(ujt, res)

## calculate quality measures

# qual = matrix(nrow = 1400, ncol = 2)
# for(i in 1:1400)  {
#   qual[i, ] = calculate_quality_measures(est = transform_data(result_2d$result[[i]][, 3:4]),
#                                          true = true2d[, 3:4], n1 = "me", n2 = "s")
# }
# 
# colnames(qual) = c("deltap2", "distall")
# 
# result_2d = cbind(result_2d, qual)
# save(result_2d, file = "result_2d_emo21_final.RData")
# load("result_2d_emo21_final.RData")
# result_2d = result_2d[, list(time.running, algorithm, n, deltap2, distall)]
# save(result_2d, file = "result_2d_emo21_final_qual.RData")

################################################################################

load("result_2d_emo21_final_qual.RData")

## Code for Figure 2b

qf = result_2d
qf = qf[, list(time.running, algorithm, n, deltap2, distall)]
qf$n = factor(qf$n, levels = c(50, 100, 200, 500, 1000, 2000, 5000))
qf$algorithm[which(qf$algorithm == "integration_2d_known")] = "k"
qf$algorithm[which(qf$algorithm == "integration_2d_unknown")] = "u"

q2ddp = plot_quality(quality = qf[, list(deltap2, n, algorithm)],
                     measure = "deltap2", group = "algorithm", method = "n",
                     f = NULL, rawdata = TRUE)
q2ddp = cbind(q2ddp, m = "Delta[p]")
q2dda = plot_quality(quality = qf[, list(distall, n, algorithm)],
                     measure = "distall", group = "algorithm", method = "n",
                     f = NULL, rawdata = TRUE)
q2dda = cbind(q2dda, m = "MED")
qualf = rbind(q2ddp, q2dda)

lab = function(labels) {
  list(list(expression(paste(Delta[p], ", MIOS-PDE")),
            expression(paste(Delta[p], ", MIOS-KDE")),
            expression(paste("MED, MIOS-PDE")),
            expression(paste("MED, MIOS-KDE"))))
}

pdf("qual_2d.pdf", width = 12, height = 4)
ggplot(data = qualf, mapping = aes(x = method, y = measure, group = method)) +
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
  facet_grid(~ m + group, labeller = lab) +
  labs(x = "Number of observations for building the metamodel (n)")
dev.off()

################################################################################

## Code for Figure 4b

q2d = qf
q2d[, ddistall := distall]
q2d[, ddeltap2 := deltap2]
q2d$ddistall[which(q2d$algorithm == "k")] = q2d$distall[which(q2d$algorithm == "k")] -
  q2d$distall[which(q2d$algorithm == "u")]
q2d$ddeltap2[which(q2d$algorithm == "k")] = q2d$deltap2[which(q2d$algorithm == "k")] -
  q2d$deltap2[which(q2d$algorithm == "u")]
q2d = q2d[c(1:700), ]

qual1 = plot_quality(quality = q2d[, list(ddeltap2, n, algorithm)], measure = "ddeltap2",
                    group = "algorithm", method = "algorithm", f = NULL, rawdata = TRUE)
qual1 = cbind(qual1, m = "Difference in Delta[p], (known - unknown)")
qual2 = plot_quality(quality = q2d[, list(ddistall, n, algorithm)], measure = "ddistall",
                     group = "algorithm", method = "algorithm", f = NULL, rawdata = TRUE)
qual2 = cbind(qual2, m = "Difference in MED (known - unknown)")
quali = rbind(qual1, qual2)

lab = function(labels) {
  list(list(expression(paste("Difference in ", Delta[p], " (MIOS-PDE - MIOS-KDE)")),
            expression(paste("Difference in MED (MIOS-PDE - MIOS-KDE)"))))
}

pdf("diff_qual_2d.pdf", width = 12, height = 4) 
ggplot(data = quali, mapping = aes(x = n, y = measure)) +
  geom_boxplot() + theme_bw() +
  labs(x = "Number of observations for building the metamodel (n)", y = "") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.title = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)) +
  facet_wrap(~m, nrow = 1, labeller = lab) +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed")
dev.off()
