# scaling correction methods based on assessor model and generalised prodcrustes analysis


# generate data to try scaling adjustment methods
# I=4, J=5, K=2
test_data = list(readings=rep(unlist(sapply(1:10, function(x) 1:5*x, simplify=F)), each=2)+rnorm(100,1), 
                 assessors=as.factor(rep(1:10, each=10)), 
                 samples=as.factor(rep(rep(1:5, 10), each=2)), 
                 replicates=as.factor(rep(1:2, 50)))
# ten_berge(test_data$readings, test_data$assessors, test_data$samples, test_data$replicates)


ten_berge = function(readings, assessors, samples, replicates, na.remove=T, na.compensate=F, prnt=F){
    #' calculate ten Berge scaling factors
    #' input: 4 vectors of equal length
    #' intended for balanced datasets
    #' readings = attribute values
    #' alculations are done based on averaged replicates
    
    df_ = data.frame(readings=readings, assessors=assessors, samples=samples, replicates=replicates)
    df_ = df_[order(assessors, samples, replicates),]
    
    n_assessors = length(unique(assessors))
    n_samples = length(unique(samples))
    n_replicates = length(unique(replicates))
    
    if(nrow(df_) < n_assessors*n_samples*n_replicates) warning('ten Berge: dataset is not balanced')

    # average replicates
    df_$readings_mean = ave(df_$readings, paste(df_$assessors, df_$samples), FUN=mean)
    df_ = df_[!duplicated(paste(df_$assessors, df_$samples)),]
    
    if(nrow(df_) < n_assessors*n_samples) stop('ten Berge: dataset misses sample-assessor combinations')
    
    # rows = assessors, cols = samples
    Y = reshape2::acast(df_, assessors~samples, value.var='readings_mean')
    # subtract assessor means
    Y = sweep(Y, 1, rowMeans(Y), "-")
    rown = rownames(Y)
    
    
    RR = sum(Y^2, na.rm=T)  # sum of squares
    
    Y = sqrt(n_assessors/RR)*Y  # scale Y so that total sqrt(SSQ) = sqrt(n_assessors)
    
    # pairwise inner products between assessors
    YY = matrix(nrow=n_assessors, ncol=n_assessors)  
    
    for(ass in seq_along(rown))
        for(ass_ in seq_along(rown))
            if(ass<=ass_)
                YY[ass==unique(assessors), ass_==unique(assessors)] = 
                    YY[ass_==unique(assessors), ass==unique(assessors)] = 
                        sum(diag(t(Y[ass==rown,]) %*% Y[ass_==rown,]), na.rm=T)
            
    
    # for scaling by assessor
    Yd = diag(diag(YY), n_assessors, n_assessors)
    
    # scaling of rows and columns
    Q = expm::sqrtm(solve(Yd)) %*% YY %*% expm::sqrtm(solve(Yd))
    
    EQ = eigen(Q)
    PQ = EQ$vectors
    
    # factor G for each assessor
    # for correction, the readings of each assessor must be multiplied with the respective G
    G = sapply(1:n_assessors, function(x) (n_assessors/Yd[x,x])**(1/2)*PQ[x,1])
    names(G) = unique(assessors)
    if(mean(G)<0) G=-G
    
    return(G)
}


assessor_mod = function(readings, assessors, samples, replicates){
    #' estimate the assessor model according to Brockhoff&Skovgaard (1994) 
    #'     doi: 10.1016/0950-3293(94)90037-X
    #' and adjust readings according to Romano et al (2008)
    #' data is assumed to be ordered by (assessor, sample, replicate)      
    
    n_samples = length(unique(samples))
    n_assessors = length(unique(assessors))
    n_replication = length(unique(replicates))
    readings_mult = readings
    readings_add = readings
    
    # initialise product effects nu_p 
    # calculate product means
    Y.p. = aggregate(readings, by=list(samples), FUN=mean)[2]
    
    # mean intensity of assessors (alpha)
    alpha = aggregate(readings, by=list(assessors), FUN=mean)[2]
    
 
    ####### Step 1:Initialize Va (called nu_p in paper) ########
    
    nu_p = as.matrix((Y.p.-mean(readings)) / sqrt(sum((Y.p.-mean(Y.p.[,1]))^2)) )
    nu_p_est = matrix(nrow=n_samples, ncol=1)
    
    #############Step 2: Paramater estimations##################
    D_V = 2
    while(D_V>1e-7){
        
        ## estimate beta_a
        Beta = matrix(nrow=n_assessors, ncol =1)
        Yap = matrix(nrow=n_assessors, ncol=n_samples)
        
        for(a in 1:n_assessors){ 
            Sum_Beta = 0
            for(p in 1:n_samples){
                Yap[a,p]= mean(readings[assessors==a & samples==p])
                Sum_Beta = Sum_Beta + Yap[a,p]*nu_p[p,]
            }
            Beta[a,] = Sum_Beta
        }
        
        ## estimate sigma2_a
        sigma2_a = matrix(nrow=n_assessors, ncol =1)
        for(a in 1:n_assessors){
            sum_PR=0
            for(p in 1:n_samples){
                sum_PR = sum_PR + 
                    sum((readings[assessors==a & samples==p]-alpha[a,]-Beta[a,]*nu_p[p,])**2) 
            }
            sigma2_a[a,] = sum_PR/(n_samples*n_replication)
        }
        
        ## eq 11 (step 3)
        for(p in 1:n_samples){
            Sum_Beta_Sigma_Y = 0
            Sum_Beta_Sigma = 0
            
            for(a in 1:n_assessors){
                Sum_Beta_Sigma  = Sum_Beta_Sigma + Beta[a,]**2/sigma2_a[a,]
                Sum_Beta_Sigma_Y = Sum_Beta_Sigma_Y + (Beta[a,]**2/sigma2_a[a,])*((Yap[a,p] - alpha[a,])/(Beta[a,]))
            }
            nu_p_est[p,] = Sum_Beta_Sigma_Y/Sum_Beta_Sigma
        }
        
        ## step 4: new estimate for nu_p; calculate change
        SV = (nu_p_est - mean(nu_p_est))/sqrt(sum( (nu_p_est - mean(nu_p_est))^2 ))
        D_V = (t(SV-nu_p)%*%(SV-nu_p))
        nu_p = SV
        if(!is.finite(D_V)) return(list(status=F))
        # catch error
    }
    
    # multiplicative scaling (Romano et al., 2008, Eq. 10)
    for(i in 1:n_assessors)
        readings_mult[assessors==i] = (readings[assessors==i] - alpha[i,])/(Beta[i,])
    
    # additive scaling (remove only scaling part of interaction) (Romano et al., 2008, Eq. 11)
    for(i in 1:n_assessors){
        for(j in 1:n_samples){
            readings_add[assessors==i & samples==j] = 
                readings[assessors==i & samples==j] - (Beta[i,]-mean(Beta))*nu_p_est[j]
        }
    }
    return(list(mult=readings_mult, add=readings_add, Beta=Beta, status=T))
}        



fit_MAM = function(sim_dat){
    ################# MAM ###################
    
    # options for MAM: 
    #  1)  adjustedMAM or normal MAM: how to deal with negative scaling?
    #  2)  conditionalMAM: use MAM or AOV, depending on p(scaling)
    # resulting options: 
    # - MAM.AOV_default: adjusted MAM/AOV with p(scaling) = 0.2 [this is the default case] - called MAM_CA in paper
    # - MAM_adj: adjusted MAM with p(scaling) = 1, thus always use MAM, and adjust SS(scaling) - called MAM_A following paper nomenclature
    # - MAM.AOV_nonadj: non-adjusted MAM/AOV with p(scaling) = 0.2 - called MAM_C in paper
    # - MAM_nonadj: non-adjusted MAM/AOV with p(scaling) = 1, thus always MAM - called MAM following paper nomenclature
    #  alpha_conditionalMAM = 0 is identical to using conventional ANOVA aov()

    fit_MAM.AOV_default = MAManalysis(sim_dat[,c(2:4,1)], adjustedMAM=T, alpha_conditionalMAM=0.2)  # required column order: assessor, product, replicate, intensity
    fit_MAM_adj = MAManalysis(sim_dat[,c(2:4,1)], adjustedMAM=T, alpha_conditionalMAM=1)  # required column order: assessor, product, replicate, intensity
    fit_MAM.AOV_nonadj = MAManalysis(sim_dat[,c(2:4,1)], adjustedMAM=F, alpha_conditionalMAM=0.2)  # required column order: assessor, product, replicate, intensity
    fit_MAM_nonadj = MAManalysis(sim_dat[,c(2:4,1)], adjustedMAM=F, alpha_conditionalMAM=1)  # required column order: assessor, product, replicate, intensity

    
    mam_f = sapply(list(fit_MAM.AOV_default, fit_MAM_adj, fit_MAM.AOV_nonadj, fit_MAM_nonadj), function(x) x$aov[,,1][2,4])
    
    mam_p = sapply(list(fit_MAM.AOV_default, fit_MAM_adj, fit_MAM.AOV_nonadj, fit_MAM_nonadj), function(x) x$aov[,,1][2,5])
    
    # record details about the MAM: p(scaling) in MAM_CA, p(scaling) in MAM_C, p(disagreement) in MAM_CA, 
    #   p(disagreement) in MAM_C, number of negative scaling factors in MAM_CA
    mam_extra = c(fit_MAM.AOV_default$aov[,,1][3,5], fit_MAM.AOV_nonadj$aov[,,1][3,5], fit_MAM.AOV_default$aov[,,1][4,5], 
                 fit_MAM.AOV_nonadj$aov[,,1][4,5], fit_MAM.AOV_default$adj.MAM.negatives)  
    
    return(list(mam_f=mam_f, mam_p=mam_p, mam_extra=mam_extra))
}
    

scaling_fit_aov = function(sim_dat) {  
    # raw
    sim_dat$raw = sim_dat$intensity
    
    
    # standardisation
    sim_dat$std = (sim_dat$intensity - ave(sim_dat$intensity, sim_dat$assessor, FUN=function(x) mean(x))) / ave(sim_dat$intensity, sim_dat$assessor, FUN=function(x) sd(x))
    
    
    # Ten Berge 
    # calculate scaling factors for assessors
    tenb_facs = ten_berge(sim_dat$intensity, sim_dat$assessor, sim_dat$sample, sim_dat$repl)
    
    # transform according to tenberge factors
    sim_dat$tenb = sim_dat$intensity - ave(sim_dat$intensity, sim_dat$assessor, FUN=mean)  # subtract means
    sim_dat$tenb = (tenb_facs[sim_dat$assessor]*sim_dat$tenb) + mean(sim_dat$intensity)  # scale and restore original assessor means
    
    
    # Assessor model
    asm = assessor_mod(sim_dat$intensity, sim_dat$assessor, sim_dat$sample, sim_dat$repl)
    if(asm$status==F) {print('error in assessor model'); next}  # error in estimating assessor model: try new permutation
    sim_dat$add = asm$add
    sim_dat$mult = asm$mult

    
    method_names = c('raw', 'std', 'tenb', 'add', 'mult')
    
    # ANOVA 
    aov_fits = sapply(method_names, function(x) summary(aov(sim_dat[,x] ~ sim_dat$sample + Error(sim_dat$assessor/sim_dat$sample))), simplify=F)
    aov_f = sapply(aov_fits, function(x) x[[2]][[1]][1,'F value'])
    aov_p = sapply(aov_fits, function(x) x[[2]][[1]][1,'Pr(>F)'])
    return(list(aov_f=aov_f, aov_p=aov_p, status=asm$status))
}


