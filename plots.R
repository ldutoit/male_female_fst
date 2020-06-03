load("vectors_reanalysis.RData")

#error function
erf <- function(x){2*pnorm(x * sqrt(2)) - 1}

README()

#cumulative distribution plots (SuppMat Figs.)
n.H = 2/(1/94 + 1/94)
plot(sort(fly_fst_observed), erf(sqrt(n.H*sort(fly_fst_observed))), type = "l", xlim = c(0, 3/n.H), lwd = 3, col = "gray", xlab = "F.ST", ylab = "CDF")
lines(sort(fly_fst_permuted), 1:length(fly_fst_permuted)/length(fly_fst_permuted))
lines(sort(fly_fst_observed), 1:length(fly_fst_observed)/length(fly_fst_observed), col = "RED")
pdf("fish_4lines.pdf")
n.H = 2/(1/114 + 1/334)
plot(sort(fish_fst_permuted), erf(sqrt(n.H*sort(fish_fst_permuted))), type = "l", xlim = c(0, 3/n.H), lwd = 3, col = "gray", xlab = "F.ST", ylab = "CDF") # For humans and flycatcher, the null is the same for every SNP, but not for the pipefish as the number of missing data varies position by position, we created one vector that draws one point of the null for every SNP, creating a correct null
lines(sort(fish_theoretical_null), 1:length(fish_fst_permuted)/length(fish_fst_permuted), type = "l", xlim = c(0, 3/n.H), lwd = 1, col = "gray", xlab = "F.ST", ylab = "CDF",lty=2)
lines(sort(fish_fst_permuted), 1:length(fish_fst_permuted)/length(fish_fst_permuted))
lines(sort(fish_fst_observed), 1:length(fish_fst_observed)/length(fish_fst_observed), col = "RED")
dev.off()

n.H = 2/(1/2542 + 1/2466)
plot(sort(human_fst_observed), erf(sqrt(n.H*sort(human_fst_observed))), type = "l", xlim = c(0, 3/n.H), lwd = 3, col = "gray", xlab = "F.ST", ylab = "CDF")
lines(sort(human_fst_permuted), 1:length(human_fst_permuted)/length(human_fst_permuted))
lines(sort(human_fst_observed), 1:length(human_fst_permuted)/length(human_fst_observed), col = "RED")

n.H = 2/(1/2542 + 1/2466)
plot(sort(human_fst_observed_nomafFilter), erf(sqrt(n.H*sort(human_fst_observed_nomafFilter))), type = "l", xlim = c(0, 3/n.H), lwd = 3, col = "gray", xlab = "F.ST", ylab = "CDF")
lines(sort(human_fst_permuted_nomafFilter), 1:length(human_fst_permuted_nomafFilter)/length(human_fst_permuted_nomafFilter))
lines(sort(human_fst_observed_nomafFilter), 1:length(human_fst_permuted_nomafFilter)/length(human_fst_observed_nomafFilter), col = "RED")

#quantile plots (SuppMat Figs.)
plot(100*fly_percquantile_observed, pch = 20, col = "RED",ylim = c(0, 10))
points(100*fly_percquantile_permuted, pch = 20)
lines(rep(1, length(fly_percquantile_permuted)), col = "GRAY", lty = 2)

plot(100*fish_percquantile_observed, pch = 20, col = "RED", ylim = c(0, 10))
points(100*fish_percquantile_permuted, pch = 20)
lines(rep(1, length(fly_percquantile_permuted)), col = "GRAY", lty = 2)
plot(100*human_percquantile_observed, pch = 20, col = "RED", ylim = c(0,10))
points(100*human_percquantile_permuted, pch = 20)
lines(rep(1, length(fly_percquantile_permuted)), col = "GRAY", lty = 2)

plot(100*human_percquantile_observed_nomafFilter, pch = 20, col = "RED", ylim = c(0,10))
points(100*human_percquantile_permuted_nomafFilter, pch = 20)
lines(rep(1, length(fly_percquantile_permuted)), col = "GRAY", lty = 2)


#ratio plots (Fig. 2)
plot(fly_percquantile_observed/fly_percquantile_permuted, ylim = c(0,2), pch = 19)
lines(0:100, rep(1, 101))

plot(fish_percquantile_observed/fish_percquantile_permuted, ylim = c(0,2), pch = 19)
lines(0:100, rep(1, 101))

plot(human_percquantile_observed/human_percquantile_permuted, ylim = c(0,2), pch = 19)
lines(0:100, rep(1, 101))

