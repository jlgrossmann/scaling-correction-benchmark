sim_results = readRDS('out_permutation_all_cheese-Fattiness_2022-11-23_20.21.54.Rds')


if(sim_results$parameters$perm){
  # calculate p-values from permutation
  # with less than 10 exceedances, the p-value becomes less accurate
  pvals = sapply(as.data.frame(sim_results$f_table), function(x) {
    n_exceedances = sum(x[-1]>x[1], na.rm=T)
    data.frame(p_ecdf=signif((1+n_exceedances) / sum(!is.na(x))), exceedances=paste0(as.character(n_exceedances), {if(n_exceedances<10) ' (!)'}))}, simplify=F)
  pvals = do.call(rbind, pvals)
  pvals
}


# Figure 1
histobreaks = 0:200/20
pdf(file=paste0('Fig1.pdf'), height=6.6, width=4.1)
par(lend=1, mgp=c(2,.8,0), mar=c(4,3.8,3.3,1)*0.75, mfrow=c(2,1))
#F 
f_density = sapply(as.data.frame(sim_results$f_table), function(x) hist(x[x<10], breaks=histobreaks, plot=F)[c('density', 'mids')], simplify=F)
matplot(f_density$raw$mids, sapply(f_density, function(x) x$density)[,1:6], type='l', lty=1, xlab=expression(F[nu]), ylab='density', 
        main=expression(F[nu]~"null distributions (unrestricted)"), xlim=c(0.15,4.85), ylim=c(0.02,1.1), col=0)
polygon(f_density$raw$mids, df(f_density$raw$mids, df1=sim_results$parameters$n_samp-1, df2=(sim_results$parameters$n_samp-1)*(sim_results$parameters$n_ass-1)), col='#d4d4d4', lty=0)
matplot(f_density$raw$mids, sapply(f_density, function(x) x$density)[,1:6], type='l', lty=1, add=T, lwd=1.25)  
legend('topright', legend=c('raw', 'std', 'tenb', 'add', 'mult', expression(MAM[CA])), lty=1, col=1:6, cex=0.75, bg="white")


# p
histobreaks2 = 0:100/100
p_density = sapply(as.data.frame(sim_results$p_table), function(x) hist(x, breaks=histobreaks2, plot=F)$density)
matplot(histobreaks2[-1]-histobreaks2[2]/2, p_density[,1:6], type='l', lty=1, xlab=expression(p[nu]), ylab='density', 
        main=expression(p[nu]~"distributions (unrestricted)"), xlim=c(0.02,0.98), ylim=c(0,3), col=0)
polygon(c(0,0,1,1), c(0,1,1,0), col='#d4d4d4', lty=0)
matplot(histobreaks2[-1]-histobreaks2[2]/2, p_density[,1:6], type='l', lty=1, ylim=c(0,3), add=T)
legend('topright', legend=c('raw', 'std', 'tenb', 'add', 'mult', expression(MAM[CA])), lty=1, col=1:6, cex=0.75, bg="white")

dev.off()
