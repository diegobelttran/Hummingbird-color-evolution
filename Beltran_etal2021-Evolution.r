###### Script to replicate analyses in 
Speciation rates are positively correlated with the rate of plumage color evolution in hummingbirds
Beltrán D.F., Shultz A.J., Parra J.L., Evolution (2021)


library(ape)
library(RRphylo)

#### color descriptor data
vsM <- read.table('.../vsM.txt',h=T,sep='\t')
rownames(vsM) <- vsM$Species
vsF <- read.table('~/Documents/Datos_colibries/BaseDeDatos/vsF.txt',h=T,sep='\t')
rownames(vsF) <- vsF$Species

# keep common species for both sexes
vsM2 <- vsM[-which(is.na(match(vsM$Species,vsF$Species))),]
vsF2 <- vsF[-which(is.na(match(vsF$Species,vsM$Species))),]

#### phylogenetic tree with species in common for both sexes
phyU <- read.tree('.../treeU.tre')

# keep only species in the tree
vsM2 <- vsM2[-which(is.na(match(vsM2$Species,phyU$tip.label))),]
# keep only species in the tree
vsF2 <- vsF2[-which(is.na(match(vsF2$Species,phyU$tip.label))),]

### male evolution rates for each trait using RRphylo
ratesM <- data.frame(matrix(NA,ncol=dim(vsM2)[2],nrow=dim(vsM2)[1]))
names(ratesM) <- names(vsM2)
ratesM$Species <- vsM2$Species

for (i in 1:(dim(vsM2)[2]-1)){
trait <- setNames(scale(vsM2[,i+1])[,1],vsM2[,'Species'])
rr <- RRphylo(phyU,trait)

ratesM[,i+1] <- rr$rates[match(names(trait),rownames(rr$rates))]
}


### female evolution rates for each variable using RRphylo
ratesF <- data.frame(matrix(NA,ncol=dim(vsF2)[2],nrow=dim(vsF2)[1]))
names(ratesF) <- names(vsF2)
ratesF$Species <- vsF2$Species

for (i in 1:(dim(vsF2)[2]-1)){
trait <- setNames(scale(vsF2[,i+1])[,1],vsF2[,'Species'])
rr <- RRphylo(phyU,trait)

ratesF[,i+1] <- rr$rates[match(names(trait),rownames(rr$rates))]
}




#### ES-sim with ClaDS and reduced-tree trait evolution rates MALES
## ES-sim script adapted from: 
# Harvey M.G. and Rabosky D.L. 2018. Continuous traits and speciation rates: Alternatives to state-dependent diversification models. Methods in Ecology and Evolution 9: 984– 993. 

library(RPANDA)

# ClaDS run
load('.../RunTrial2_300000')

# extract the Maximum A Posteriori from the marginal posterior distribution
maps <- getMAPS_ClaDS(sp_sampler_14, burn = 1/4, thin = 1)
   
sp_clad <- setNames(maps[(which(phyU$edge[,2] <= length(phyU$tip.label)))+4],phyU$tip.label)

#### correlation between sp rates and trait evolution rates, compared against null BM model
is <- log(sp_clad) # log transform
is <- is[phyU$tip.label] # make sure 'is' is in the same order as phyU$tip.label

# table where es-sim results for each variable will be stored MALES
cladM <- data.frame(matrix(NA,ncol=9))
names(cladM) <- c('Trait','rho','upper','lower','p_value','FDR_p_value','meanSimulated','medianSimulated','varSimulated')

# table where null distribution or rho will be stored
rho <- data.frame(matrix(NA,ncol=1+10000))
names(rho)[1] <- 'Trait'

# new loop for each variable
library(mvtnorm)

for (u in 1:(dim(ratesM)[2]-1)){

# 'trait' is the log-transformed absolute value of the observed rates
trait <- setNames(log(abs(ratesM[,u+1])),ratesM$Species)
# 'traitUnt' is the observed un-transformed rate to obtain BM estimates
traitUnt <- setNames(ratesM[,u+1],ratesM$Species)

# make sure vectors are ordered
trait <- trait[phyU$tip.label]
traitUnt <- traitUnt[phyU$tip.label]
    
	# Pearsons correlation between claDS and observed trait

	res <- cor.test(is, trait, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates

	vv <- vcv.phylo(as.phylo(phyU))
	onev <- matrix(rep(1, length(traitUnt)), nrow=length(traitUnt), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% traitUnt))
	rate <- as.vector((t(traitUnt-root) %*% solve(vv) %*% (traitUnt-root))/length(traitUnt))

	# Brownian simulations 

	sims <- t(rmvnorm(10000, sigma=rate*vv))

	rownames(sims) <- rownames(vv)

	# Pearsons correlations of log-transformed absolute value of simulated datasets

	sim.r <- sapply(1:10000, function(x) cor.test(is[as.vector(rownames(sims))], log(abs(sims[,x])), method="pearson")$estimate)

	# Calculate the two-tailed p value

	corr <- res$estimate

	upper <- (length(sim.r[sim.r >= corr])+1)/(10000+1)

	lower <- (length(sim.r[sim.r <= corr])+1)/(10000+1)

	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed


	result <- as.vector(c(corr, pval))

	names(result) <- c("rho", "P Value")

# add results ofeach variable into a table
cladM[u,'Trait'] <- names(ratesM)[u+1]
cladM[u,'rho'] <- result[1]
cladM[u,'upper'] <- res$conf.int[2]
cladM[u,'lower'] <- res$conf.int[1]
cladM[u,'p_value'] <- result[2]
cladM[u,'meanSimulated'] <- mean(sim.r)
cladM[u,'medianSimulated'] <- median(sim.r)
cladM[u,'varSimulated'] <- var(sim.r)

# save null distribution
rho[u,'Trait'] <- names(ratesM)[u+1]
rho[u,2:dim(rho)[2]] <- sim.r
}

# False Discovery Rate
cladM$FDR_p_value <- p.adjust(cladM$p_value,method='fdr')


###### analyses separating positive and negative RRphylo rates MALES

# table where es-sim results for each phy will be stored
cladSignM <- data.frame(matrix(NA,ncol=7))
names(cladSignM) <- c('Trait','rhoP','p_valueP','FDRP','rhoN','p_valueN','FDRN')

# new loop for each variable

for (u in 1:(dim(ratesM)[2]-1)){

# separate positive and negative (absolute value for these) evolution rates
traitP <- setNames(log(ratesM[which(ratesM[,u+1] > 0),u+1]),ratesM[which(ratesM[,u+1] > 0),'Species'])
traitP_Unt <- setNames(ratesM[which(ratesM[,u+1] > 0),u+1],ratesM[which(ratesM[,u+1] > 0),'Species'])
traitN <- setNames(log(abs(ratesM[which(ratesM[,u+1] < 0),u+1])),ratesM[which(ratesM[,u+1] < 0),'Species'])
traitN_Unt <- setNames(ratesM[which(ratesM[,u+1] < 0),u+1],ratesM[which(ratesM[,u+1] < 0),'Species'])

# prune tree and es according to sign
phyP <- drop.tip(phyU,phyU$tip.label[which(is.na(match(phyU$tip.label,ratesM[which(ratesM[,u+1] > 0),'Species'])))])
isP<- is[phyP$tip.label]

phyN <- drop.tip(phyU,phyU$tip.label[which(is.na(match(phyU$tip.label,ratesM[which(ratesM[,u+1] < 0),'Species'])))])
isN<- is[phyN$tip.label]

# re-order traitP and tratN
traitP <- traitP[phyP$tip.label]
traitP_Unt <- traitP_Unt[phyP$tip.label]
traitN <- traitN[phyN$tip.label]
traitN_Unt <- traitN_Unt[phyN$tip.label]



	# Pearson's correlation between splits statistic and trait POSITIVE

	resP <- cor.test(isP, traitP, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates

	vv <- vcv.phylo(as.phylo(phyP))
	onev <- matrix(rep(1, length(traitP_Unt)), nrow=length(traitP_Unt), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% traitP_Unt))
	rate <- as.vector((t(traitP_Unt-root) %*% solve(vv) %*% (traitP_Unt-root))/length(traitP_Unt))

	# Brownian simulations 

	sims <- t(rmvnorm(10000, sigma=rate*vv))

	rownames(sims) <- rownames(vv)

	# Pearson's correlations of simulated datasets

	sim.r <- sapply(1:10000, function(x) cor.test(isP[as.vector(rownames(sims))], log(abs(sims[,x])), method="pearson")$estimate)

	# Calculate the two-tailed p value

	corr <- resP$estimate

	upper <- (length(sim.r[sim.r >= corr])+1)/(10000+1)

	lower <- (length(sim.r[sim.r <= corr])+1)/(10000+1)

	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed


	resultP <- as.vector(c(corr, pval))

	names(resultP) <- c("rho", "P Value")




	# Pearson's correlation between splits statistic and trait NEGATIVE

	resN <- cor.test(isN, traitN, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates

	vv <- vcv.phylo(as.phylo(phyN))
	onev <- matrix(rep(1, length(traitN_Unt)), nrow=length(traitN_Unt), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% traitN_Unt))
	rate <- as.vector((t(traitN_Unt-root) %*% solve(vv) %*% (traitN_Unt-root))/length(traitN_Unt))

	# Brownian simulations 

	sims <- t(rmvnorm(10000, sigma=rate*vv))

	rownames(sims) <- rownames(vv)

	# Pearson's correlations of simulated datasets

	sim.r <- sapply(1:10000, function(x) cor.test(isN[as.vector(rownames(sims))], log(abs(sims[,x])), method="pearson")$estimate)

	# Calculate the two-tailed p value

	corr <- resN$estimate

	upper <- (length(sim.r[sim.r >= corr])+1)/(10000+1)

	lower <- (length(sim.r[sim.r <= corr])+1)/(10000+1)

	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed


	resultN <- as.vector(c(corr, pval))

	names(resultN) <- c("rho", "P Value")



# add results ofeach variable into a table
cladSignM[u,'Trait'] <- names(ratesM)[u+1]
cladSignM[u,'rhoP'] <- resultP[1]
cladSignM[u,'p_valueP'] <- resultP[2]

cladSignM[u,'rhoN'] <- resultN[1]
cladSignM[u,'p_valueN'] <- resultN[2]
}

cladSignM[,'FDRP'] <- p.adjust(cladSignM$p_valueP,method='fdr')
cladSignM[,'FDRN'] <- p.adjust(cladSignM$p_valueN,method='fdr')




###### table where es-sim results for each variable will be stored FEMALES
cladF <- data.frame(matrix(NA,ncol=9))
names(cladF) <- c('Trait','rho','upper','lower','p_value','FDR_p_value','meanSimulated','medianSimulated','varSimulated')

# table where null distribution or rho will be stored
rho <- data.frame(matrix(NA,ncol=1+10000))
names(rho)[1] <- 'Trait'

# new loop for each variable
library(mvtnorm)

for (u in 1:(dim(ratesF)[2]-1)){

trait <- setNames(log(abs(ratesF[,u+1])),ratesF$Species)
traitUnt <- setNames(ratesF[,u+1],ratesF$Species)

trait <- trait[phyU$tip.label]
traitUnt <- traitUnt[phyU$tip.label]
    

	# Pearsons correlation between claDS and observed trait

	res <- cor.test(is, trait, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates

	vv <- vcv.phylo(as.phylo(phyU))
	onev <- matrix(rep(1, length(traitUnt)), nrow=length(traitUnt), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% traitUnt))
	rate <- as.vector((t(traitUnt-root) %*% solve(vv) %*% (traitUnt-root))/length(traitUnt))

	# Brownian simulations 

	sims <- t(rmvnorm(10000, sigma=rate*vv))

	rownames(sims) <- rownames(vv)

	# Pearsons correlations of log-transformed absolute value of simulated datasets

	sim.r <- sapply(1:10000, function(x) cor.test(is[as.vector(rownames(sims))], log(abs(sims[,x])), method="pearson")$estimate)

	# Calculate the two-tailed p value

	corr <- res$estimate

	upper <- (length(sim.r[sim.r >= corr])+1)/(10000+1)

	lower <- (length(sim.r[sim.r <= corr])+1)/(10000+1)

	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed


	result <- as.vector(c(corr, pval))

	names(result) <- c("rho", "P Value")

# add results ofeach variable into a table
cladF[u,'Trait'] <- names(ratesF)[u+1]
cladF[u,'rho'] <- result[1]
cladF[u,'upper'] <- res$conf.int[2]
cladF[u,'lower'] <- res$conf.int[1]
cladF[u,'p_value'] <- result[2]
cladF[u,'meanSimulated'] <- mean(sim.r)
cladF[u,'medianSimulated'] <- median(sim.r)
cladF[u,'varSimulated'] <- var(sim.r)

# save null distribution
rho[u,'Trait'] <- names(ratesF)[u+1]
rho[u,2:dim(rho)[2]] <- sim.r
}

# False Discovery Rate
cladF$FDR_p_value <- p.adjust(cladF$p_value,method='fdr')



#### analyses separating positive and negative RRphylo values FEMALES

# table where es-sim results for each phy will be stored
cladSignF <- data.frame(matrix(NA,ncol=7))
names(cladSignF) <- c('Trait','rhoP','p_valueP','FDRP','rhoN','p_valueN','FDRN')

# new loop for each variable

for (u in 1:(dim(ratesF)[2]-1)){

# separate positive and negative (absolute value for these) evolution rates
traitP <- setNames(log(ratesF[which(ratesF[,u+1] > 0),u+1]),ratesF[which(ratesF[,u+1] > 0),'Species'])
traitP_Unt <- setNames(ratesF[which(ratesF[,u+1] > 0),u+1],ratesF[which(ratesF[,u+1] > 0),'Species'])
traitN <- setNames(log(abs(ratesF[which(ratesF[,u+1] < 0),u+1])),ratesF[which(ratesF[,u+1] < 0),'Species'])
traitN_Unt <- setNames(ratesF[which(ratesF[,u+1] < 0),u+1],ratesF[which(ratesF[,u+1] < 0),'Species'])

# prune tree and es according to sign
phyP <- drop.tip(phyU,phyU$tip.label[which(is.na(match(phyU$tip.label,ratesF[which(ratesF[,u+1] > 0),'Species'])))])
isP<- is[phyP$tip.label]

phyN <- drop.tip(phyU,phyU$tip.label[which(is.na(match(phyU$tip.label,ratesF[which(ratesF[,u+1] < 0),'Species'])))])
isN<- is[phyN$tip.label]

# re-order traitP and tratN
traitP <- traitP[phyP$tip.label]
traitP_Unt <- traitP_Unt[phyP$tip.label]
traitN <- traitN[phyN$tip.label]
traitN_Unt <- traitN_Unt[phyN$tip.label]



	# Pearson's correlation between splits statistic and trait POSITIVE

	resP <- cor.test(isP, traitP, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates

	vv <- vcv.phylo(as.phylo(phyP))
	onev <- matrix(rep(1, length(traitP_Unt)), nrow=length(traitP_Unt), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% traitP_Unt))
	rate <- as.vector((t(traitP_Unt-root) %*% solve(vv) %*% (traitP_Unt-root))/length(traitP_Unt))

	# Brownian simulations 

	sims <- t(rmvnorm(10000, sigma=rate*vv))

	rownames(sims) <- rownames(vv)

	# Pearson's correlations of simulated datasets

	sim.r <- sapply(1:10000, function(x) cor.test(isP[as.vector(rownames(sims))], log(abs(sims[,x])), method="pearson")$estimate)

	# Calculate the two-tailed p value

	corr <- resP$estimate

	upper <- (length(sim.r[sim.r >= corr])+1)/(10000+1)

	lower <- (length(sim.r[sim.r <= corr])+1)/(10000+1)

	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed


	resultP <- as.vector(c(corr, pval))

	names(resultP) <- c("rho", "P Value")




	# Pearson's correlation between splits statistic and trait NEGATIVE

	resN <- cor.test(isN, traitN, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates

	vv <- vcv.phylo(as.phylo(phyN))
	onev <- matrix(rep(1, length(traitN_Unt)), nrow=length(traitN_Unt), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% traitN_Unt))
	rate <- as.vector((t(traitN_Unt-root) %*% solve(vv) %*% (traitN_Unt-root))/length(traitN_Unt))

	# Brownian simulations 

	sims <- t(rmvnorm(10000, sigma=rate*vv))

	rownames(sims) <- rownames(vv)

	# Pearson's correlations of simulated datasets

	sim.r <- sapply(1:10000, function(x) cor.test(isN[as.vector(rownames(sims))], log(abs(sims[,x])), method="pearson")$estimate)

	# Calculate the two-tailed p value

	corr <- resN$estimate

	upper <- (length(sim.r[sim.r >= corr])+1)/(10000+1)

	lower <- (length(sim.r[sim.r <= corr])+1)/(10000+1)

	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed


	resultN <- as.vector(c(corr, pval))

	names(resultN) <- c("rho", "P Value")



# add results ofeach variable into a table
cladSignF[u,'Trait'] <- names(ratesF)[u+1]
cladSignF[u,'rhoP'] <- resultP[1]
cladSignF[u,'p_valueP'] <- resultP[2]

cladSignF[u,'rhoN'] <- resultN[1]
cladSignF[u,'p_valueN'] <- resultN[2]
}

cladSignF[,'FDRP'] <- p.adjust(cladSignF$p_valueP,method='fdr')
cladSignF[,'FDRN'] <- p.adjust(cladSignF$p_valueN,method='fdr')

