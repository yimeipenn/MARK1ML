# Application of methods described in "Maximum Likelihood Based 
# Analysis of Equally Spaced Longitudinal Count Data with Specified 
# Marginal Means, First-order Antedependence, and Linear Conditional 
# Expectations

####################################################################
## Options
####################################################################

## Optional code for more decimal places
options(digits=10)

library(alabama)

####################################################################
## Data Processsing Functions
####################################################################

## Getting Info from the Data
## This function will remain all the subjects
## This function will not help order the subjects
## This function was written by Matt Guerra

cluster.size = function(id)
	{
 	clid = unique(id)
  	m = length(unique(id))
  	n = rep(0,m)
  	autotime = rep(0,0)
  	for(i in 1:m)
		{
    		n[i] = length(which(id == clid[i]))
    		autotime = c(autotime, 1:n[i])
  		}
  	id = rep(1:m, n)
  	return(list(m = m, n = n, id = id, autotime = autotime))
	}

## Data Process: This function will delete subjects with less or 
## equal to #=del.n observations
## This function was written by Matt Guerra with additions by Shaun Bender

data.proc = function(data,formula,time=NULL,id,del.n, binom = NULL)
	{
  	dat = data.frame(data)
  	col.name = names(dat)
  
  	cluster = cluster.size(id)
  	m = cluster$m
  	n = cluster$n
  	id = cluster$id
  	if(length(time)==0)
		{
    		time = cluster$autotime
  		}
  	autotime = cluster$autotime
  	index = order(id,time)  
  	if(ncol(dat) == 1)
		{
    		dat = dat[index,]
  		} else
		{
    		dat = dat[index,]
  		}
  	dat = data.frame(dat)
  	names(dat) = col.name
	if(Dist == "Binomial"){binomN = binom[index]} else 
		{binomN = NULL}
  
  	del = which(n <= del.n)
  	if(length(del) > 0)
		{
    		n = n[-del]
   	 	m = length(n)
    		mtch = match(id, del)
    		del.id = which(mtch != "NA")
    		dat = dat[-del.id,]
    		dat = data.frame(dat)
    		names(dat) = col.name
    		row.names(dat) = 1:nrow(dat)
    		time = time[-del.id]
    		autotime = autotime[-del.id]
    		id = rep(1:m, n)
		if(Dist == "Binomial"){binomN = binomN[-del.id]}
  		}
  
  	formula = as.formula(formula)
  	fml = as.formula(paste("~", formula[3], "+", formula[2], sep="")) 
  	dat = model.matrix(fml, data=dat)
  
  	return(list(data = dat, time = time, autotime = autotime, id = id, m = m, n = n, binomN = binomN))
	}

####################################################################
## Log likelihood
####################################################################

## Log Likelihood function.  The function returns the value of the log
## likelihood for the inputted parameter values.
## This function was written by Victoria Gamerman, with additions by
## Shaun Bender

drv.logl = function(start.values)
	{
	if(Dist == "Negative-Binomial"){Anc = 1}
	if(Dist != "Negative-Binomial"){Anc = 0}

	if(CorrStr == "AR(1)" | CorrStr == "Markov")
		{
 		alpha = start.values[1]
		beta = start.values[2:(length(start.values)-Anc)]
		if(Dist == "Negative-Binomial"){r = start.values[length(start.values)]}
		}
	if(CorrStr == "AD(1)")
		{
		alpha = start.values[1:max(n)-1]
		beta = start.values[max(n):(length(start.values)-Anc)]
		if(Dist == "Negative-Binomial"){r = start.values[length(start.values)]}
		}

	LogLik = 0

	for (i in 1:m)
		{
        	data_i = matrix(NA, nrow=n[i], ncol=dim(dataset$data)[2])
        	data_i[1:n[i],1:dim(dataset$data)[2]] = dataset$data[which(id==i),]
        	data.end = ncol(data_i)
        	x_i = matrix(NA, nrow=n[i], ncol=k+1)
        	x_i[1:n[i],1:(k+1)] = data_i[,-data.end]
        	y_i = data_i[,data.end]
        	n_i = nrow(data_i)
		time_i = dataset$time[which(id==i)]

		for (j in 1:n_i)
			{
			if (j == 1)
				{
				lam_ij = LinkInv(i, j, beta, x_i)

				if(Dist == "Negative-Binomial"){LogLik = LogLik + UnitLikelihood(i, j, y_i, lam_ij, r)}
				else{LogLik = LogLik + UnitLikelihood(i, j, y_i, lam_ij)}
				}
			if (j > 1)
				{
				lam = c(LinkInv(i, j, beta, x_i), LinkInv(i, j-1, beta, x_i))
				lamdot_i2 = MuSt_ij(i, j, y_i, time_i, alpha, lam, r)

				if(Dist == "Negative-Binomial"){LogLik = LogLik + UnitLikelihood(i, j, y_i, lamdot_i2, r)}
				else{LogLik = LogLik + UnitLikelihood(i, j, y_i, lamdot_i2)}
				}
			}
		}

	return(-LogLik)
	}

####################################################################
## Gradient
####################################################################

## Gradient function: It should take arguments matching those of 
## f and return a vector containing the gradient. 
## This function was written by Victoria Gamerman, with additions by
## Shaun Bender

drv.grad = function(start.values)
	{
	D_Beta = matrix(0, nrow = k+1, ncol = 1)
	D_R = matrix(0, nrow = 1, ncol = 1)

	if(Dist == "Negative-Binomial"){Anc = 1}
	if(Dist != "Negative-Binomial"){Anc = 0}

	if(CorrStr == "AR(1)" | CorrStr == "Markov")
		{
		alpha = start.values[1]
		beta = start.values[2:(length(start.values)-Anc)]
		if(Dist == "Negative-Binomial"){r = start.values[length(start.values)]}
		if(Dist != "Negative-Binomial"){r = NULL}
		D_Alpha = matrix(0, nrow = 1, ncol = 1)
		}
	if(CorrStr == "AD(1)")
		{
		alpha = start.values[1:max(n)-1]
		beta = start.values[max(n):(length(start.values)-Anc)]
		if(Dist == "Negative-Binomial"){r = start.values[length(start.values)]}
		if(Dist != "Negative-Binomial"){r = NULL}
		D_Alpha = matrix(0, nrow = max(n)-1, ncol = 1)
		}

	for (i in 1:m)
		{
		data_i = matrix(NA, nrow = n[i], ncol = dim(dataset$data)[2])
		data_i[1:n[i],1:dim(dataset$data)[2]] = dataset$data[which(id==i),]
		data.end = ncol(data_i)
		x_i = matrix(NA, nrow = n[i], ncol = k+1)
		x_i[1:n[i],1:(k+1)] = data_i[,-data.end]
		y_i = data_i[,data.end]
		n_i = nrow(data_i)
		time_i = dataset$time[which(id==i)]
		if(n_i>=1)
			{
			for(j in 1:n_i)
				{
				if(j == 1)
					{
					lam_ij = LinkInv(i, j, beta, x_i)
					if(Dist == "Negative-Binomial")
						{
						D_Beta = D_Beta + functBeta(i, j, x_i, y_i, time_i, alpha, lam_ij, r)
						D_R = D_R + functR(i, j, r, y_i, time_i, lam_ij, alpha)
						}
					if(Dist != "Negative-Binomial")
						{
						D_Beta = D_Beta + functBeta(i, j, x_i, y_i, time_i, alpha, lam_ij)
						}
					}
				if(j > 1)
					{
					lam = c(LinkInv(i, j, beta, x_i), LinkInv(i, j-1, beta, x_i))
					D_Alpha = D_Alpha + functAlpha(i, j, y_i, time_i, alpha, lam, r)
					if(Dist == "Negative-Binomial")
						{
						D_Beta = D_Beta + functBeta(i, j, x_i, y_i, time_i, alpha, lam, r)
						D_R = D_R + functR(i, j, r, y_i, time_i, lam, alpha)
						}
					if(Dist != "Negative-Binomial")
						{
						D_Beta = D_Beta + functBeta(i, j, x_i, y_i, time_i, alpha, lam)
						}
					}
				}
			}
		}

	Output = t(t(c(-D_Alpha, -D_Beta)))
	if(Dist == "Negative-Binomial"){Output = t(t(c(-D_Alpha, -D_Beta, -D_R)))}
	
	return(-Output)
	}

####################################################################
## Accessory functions used by Log-Likelihood and Gradient
####################################################################

## This function calculates the inverse of the link function
## This function was written by Shaun Bender

LinkInv = function(i, j, beta, x_i)
	{
	if(Dist == "Poisson")
		{
		lam_ij = exp(t(beta)%*%x_i[j,]) 
		}
	if(Dist == "Negative-Binomial")
		{
		lam_ij = exp(t(beta)%*%x_i[j,])
		if(lam_ij[1] == 0){lam_ij[1] = .001}
		}
	if(Dist == "Binomial")
		{
		N = binomN[which(id==i)][j]
		lam_ij = N * exp(t(beta)%*%x_i[j,]) / (1 + exp(t(beta)%*%x_i[j,]))
		}
	
	return(lam_ij[1])
	}

## This function calculates the unit log likelihood for a given assumed distribution
## This function was written by Shaun Bender

UnitLikelihood = function(i, j, y_i, lam_ij, r)
	{
	if(Dist == "Poisson")
		{
		Theta = log(lam_ij)
		B_Theta = lam_ij
		if(y_i[j] == 0){Const = 0}
		if(y_i[j] > 0){Const = sum(log(seq(from = 1, to = y_i[j], by = 1)))}
		}
	if(Dist == "Negative-Binomial")
		{
		Theta = log(r * lam_ij / (1 + r * lam_ij))
		B_Theta = 1 / r * log(1 + r * lam_ij)
		Const = -lgamma(y_i[j] + 1/r) + lgamma(y_i[j] + 1) + lgamma(1/r)
		}
	if(Dist == "Binomial")
		{
		N = binomN[which(id==i)][j]
		Theta = log(lam_ij / (N - lam_ij))
		B_Theta = -N * log((N - lam_ij) / N)
		Const = -log(choose(N, y_i[j]))
		}

	Unit = y_i[j]*Theta - B_Theta - Const
	
	return(Unit)
	}

## This function calculates the derivative of the log likelihood with respect to r, under the
## assumption of a Negative-Binomial distribution

functR = function(i, j, r, y_i, time_i, lam, alpha)
	{
	lam_ij = lam[1]
	
	if(j == 1){DTheta = 1 / r / (1 + r * lam_ij)}
	if(j > 1)
		{
		lam_ij_1 = lam[2]
	
		Corr = FindCorr(j, time_i, alpha)
		C_ij = Corr[1]
		if(j>2){C_ij_1 = Corr[2]}

		E_V = FindE_V(i, j, lam, Corr, r)
		E_V_ij_1 = E_V[1]
		E_V_ij = E_V[2]
		
		if(j == 2){DE_V_ij_1 = lam_ij_1^2 / 2 / sqrt(lam_ij_1 + r * lam_ij_1^2)}
		if(j > 2){DE_V_ij_1 = 1/2/sqrt(E_V_ij_1) * (lam_ij_1^2+lam_ij_1*C_ij_1^2/(1 - C_ij_1^2)) / 
						(1-r*C_ij_1^2/(1-C_ij_1^2))^2}
		DE_V_ij = 1/2/sqrt(E_V_ij) * (lam_ij^2+lam_ij*C_ij^2/(1 - C_ij^2)) / (1-r*C_ij^2/(1-C_ij^2))^2

		Num = sqrt(E_V_ij_1) * DE_V_ij - E_V_ij * DE_V_ij_1
	
		if(j == 2){DMuSt = C_ij / sqrt(1-C_ij^2) * (y_i[j-1]-lam_ij_1) * Num / E_V_ij_1}
		if(j > 2){DMuSt = C_ij * sqrt(1-C_ij_1^2)/sqrt(1-C_ij^2)*(y_i[j-1]-lam_ij_1) * Num / E_V_ij_1}

		MuSt = MuSt_ij(i, j, y_i, time_i, alpha, lam, r)

		DTheta = (r * DMuSt + MuSt) / r / MuSt / (1 + r * MuSt)
		}

	DC = 0
	if(y_i[j] != 0){for(i in 0:(y_i[j]-1)){DC = DC + (1/r^2) * 1 / (1/r + i)}}

	Output = (y_i[j]-lam_ij) * DTheta - DC

	return(Output)
	}

## This function calculates the term to be added for the derivative of the log-likelihood
## This function was written by Shaun Bender

functAlpha = function(i, j, y_i, time_i, alpha, lam, r)
	{
	MuSt = MuSt_ij(i, j, y_i, time_i, alpha, lam, r)
	Added = (y_i[j]-MuSt)*DerivG(i, MuSt, r)*DMuStarAlpha(i, j, y_i, time_i, alpha, lam, r)

	return(Added)
	}

## This function calculates the value of Mu Star.  This function was written by
## Shaun Bender

MuSt_ij = function(i, j, y_i, time_i, alpha, lam, r)
	{
	Corr = FindCorr(j, time_i, alpha)
	Corr_ij = Corr[1]
	Corr_ij_1 = Corr[2]
	
	E_V = FindE_V(i, j, lam, Corr,r)
	E_V_ij = E_V[1]
	E_V_ij_1 = E_V[2]

	lam_ij = lam[1]
	lam_ij_1 = lam[2]
	
	if(j == 2)
		{
		MuSt = lam_ij + Corr_ij/sqrt(1-Corr_ij^2) * sqrt(E_V_ij / E_V_ij_1) * (y_i[j-1] - lam_ij_1)
		}
	if(j > 2)
		{
		MuSt = lam_ij + Corr_ij * sqrt(E_V_ij / E_V_ij_1) * sqrt((1-Corr_ij_1^2)/(1-Corr_ij^2)) * (y_i[j-1] - lam_ij_1)
		}

	constr = sqrt(lam_ij / (lam_ij_1 + lam_ij))
	if(is.finite(constr) == FALSE){constr = 0.2}
	if(is.finite(MuSt) == FALSE){MuSt = 0.5*constr}
	if(MuSt < 0){MuSt = 0.5*constr}

	return(MuSt)
	}

## This function calculates the derivative of Mu Star with respect to Alpha.
## This function was written by Shaun Bender

DMuStarAlpha = function(i, j, y_i, time_i, alpha, lam, r)
	{
	Corr = FindCorr(j, time_i, alpha)
	Corr_ij = Corr[1]
	Corr_ij_1 = Corr[2]
	
	E_V = FindE_V(i, j, lam, Corr,r)
	E_V_ij = E_V[1]
	E_V_ij_1 = E_V[2]

	lam_ij = lam[1]
	lam_ij_1 = lam[2]

	if(CorrStr == "AR(1)")
		{
		if(j == 2){Output = sqrt(E_V_ij/E_V_ij_1) * (y_i[1] - lam_ij_1) / (1-alpha^2)^(3/2)}
		if(j > 2){Output = sqrt(E_V_ij/E_V_ij_1) * (y_i[j-1] - lam_ij_1)}
		}
	if(CorrStr == "Markov")
		{
		if(j == 2)
			{
			Output = (y_i[1]-lam_ij_1) / sqrt(1-alpha^(2*time_i[2]-2*time_i[1])) * sqrt(E_V_ij/E_V_ij_1) *
					(time_i[2]-time_i[1])*alpha^(time_i[2]-time_i[1]-1) * (1 + alpha^(2*time_i[2]-2*time_i[1])/
					(1-alpha^(2*time_i[2]-2*time_i[1])))
			}
		if(j > 2)
			{
			Output = (y_i[j-1]-lam_ij_1)*sqrt(1-alpha^(2*time_i[j-1]-2*time_i[j-2]))/sqrt(1-alpha^(2*time_i[j]-2*time_i[j-1]))*
					sqrt(E_V_ij/E_V_ij_1)*alpha^(time_i[j]-time_i[j-1]-1)*((time_i[j]-time_i[j-1])/
					(1-alpha^(2*time_i[j]-2*time_i[j-1]))-(time_i[j-1]-time_i[j-2])*alpha^(2*time_i[j-1]-2*time_i[j-2])/
					(1-alpha^(2*time_i[j-1]-2*time_i[j-2])))
			}
		}
	if(CorrStr == "AD(1)")
		{
		if(j == 2)
			{
			One = c(1, rep(0,(length(alpha)-1)))
			Output = sqrt(E_V_ij/E_V_ij_1) * (y_i[1] - lam_ij_1) * One / (1-alpha[1]^2)^(3/2)
			}
		if(j > 2)
			{
			J_1 = c(rep(0,j-2),1,rep(0,length(alpha)-j+1))
			J_2 = c(rep(0,j-3),1,rep(0,length(alpha)-j+2))
			Output = sqrt(E_V_ij/E_V_ij_1) * (y_i[1] - lam_ij_1) * (J_1 * sqrt(1-alpha[j-2]^2)/(1-alpha[j-1]^2)^(3/2) - 
					J_2 * alpha[j-1] * alpha[j-2] / sqrt(1-alpha[j-1]^2) / sqrt(1-alpha[j-2]^2))
			}
		}

	return(Output)
	}

## This function calculates the term to be added for the derivative of the log-likelihood
## This function was written by Shaun Bender

functBeta = function(i, j, x_i, y_i, time_i, alpha, lam, r)
	{
	if(j == 1){Added = (y_i[j]-lam) * DerivG(i, lam, r) * Dlam_ij(j, lam, x_i)}
	if(j > 1)
		{
		MuSt = MuSt_ij(i, j, y_i, time_i, alpha, lam, r)
		Added = (y_i[j]-MuSt)*DerivG(i, MuSt, r)*DMuStarBeta(i, j, x_i, y_i, time_i, alpha, lam, r)
		}

	return(Added)
	}

## This function calculates the derivative of Mu Star with respect to Beta.
## This function was written by Shaun Bender

DMuStarBeta = function(i, j, x_i, y_i, time_i, alpha, lam,r)
	{
	Corr = FindCorr(j, time_i, alpha)
	Corr_ij = Corr[1]
	Corr_ij_1 = Corr[2]

	E_V = FindE_V(i, j, lam, Corr, r)
	E_V_ij = E_V[1]
	E_V_ij_1 = E_V[2]

	lam_ij = lam[1]
	lam_ij_1 = lam[2]

	if(j == 2)
		{
		Output = Dlam_ij(j, lam_ij, x_i) + Corr_ij/sqrt(1-Corr_ij^2)*sqrt(E_V_ij/E_V_ij_1)*((y_i[1]-lam_ij_1)/2*
				(DE_V(i,j,lam_ij,x_i,time_i,alpha,r)/E_V_ij - DE_V(i,j-1,lam_ij_1,x_i,time_i,alpha,r)/E_V_ij_1) - 
				Dlam_ij(j-1, lam_ij_1, x_i))
		}
	if(j > 2)
		{
		Output = Dlam_ij(j, lam_ij, x_i) + Corr_ij*sqrt(1-Corr_ij_1^2)/sqrt(1-Corr_ij^2)*sqrt(E_V_ij/E_V_ij_1)*
				((y_i[j-1]-lam_ij_1)/2*(DE_V(i,j,lam_ij,x_i,time_i,alpha,r)/E_V_ij - 
				DE_V(i,j-1,lam_ij_1,x_i,time_i,alpha,r)/E_V_ij_1)-Dlam_ij(j-1, lam_ij_1, x_i))
		}

	return(Output)
	}

## This function calculates the correlation structure at time j (and j-1 for j>2)
## This function was written by Shaun Bender

FindCorr = function(j, time_i, alpha)
	{
	if(CorrStr == "AR(1)")
		{
		Corr_ij = alpha
		if(j > 2){Corr_ij_1 = alpha}
		}
	if(CorrStr == "Markov")
		{
		Corr_ij = alpha^(time_i[j]-time_i[j-1])
		if(j > 2){Corr_ij_1 = alpha^(time_i[j-1]-time_i[j-2])}
		}
	if(CorrStr == "AD(1)")
		{
		Corr_ij = alpha[j-1]
		if(j > 2){Corr_ij_1 = alpha[j-2]}
		}

	if(j == 2){Output = Corr_ij}
	if(j > 2){Output = c(Corr_ij, Corr_ij_1)}

	return(Output)
	}

## This function calculates the E(Var(Y_ij|Y_ij-1)) (and if j>2, E(Var(Y_ij_1|Y_ij_2)))
## This function was written by Shaun Bender

FindE_V = function(i, j, lam, Corr, r)
	{

	lam_ij = lam[1]
	lam_ij_1 = lam[2]
	
	Corr_ij = Corr[1]
	if(j > 2){Corr_ij_1 = Corr[2]}

	if(Dist == "Poisson")
		{
		E_V_ij = lam_ij
		E_V_ij_1 = lam_ij_1
		}
	if(Dist == "Negative-Binomial")
		{
		C_ij = 1 - r * Corr_ij^2 / (1 - Corr_ij^2)
		if(j > 2){C_ij_1 = 1 - r * Corr_ij_1^2 / (1 - Corr_ij_1^2)}

		E_V_ij = (lam_ij + r * lam_ij^2) / C_ij
		if(j == 2){E_V_ij_1 = lam_ij_1 + r * lam_ij_1^2}
		if(j > 2){E_V_ij_1 = (lam_ij_1 + r * lam_ij_1^2) / C_ij_1}
		}
	if(Dist == "Binomial")
		{
		N = binomN[which(id==i)][j]
		if(j > 1){N_1 = binomN[which(id==i)][j-1]}
		
		C_ij = 1 + Corr_ij^2 / (1 - Corr_ij^2) / N
		if(j > 2){C_ij_1 = 1 + Corr_ij_1^2 / (1 - Corr_ij_1^2) / N_1}
		
		E_V_ij = lam_ij * (N - lam_ij) / N / C_ij
		if(j == 2){E_V_ij_1 = lam_ij_1 * (N_1 - lam_ij_1) / N_1}
		if(j > 2){E_V_ij_1 = lam_ij_1 * (N_1 - lam_ij_1) / N_1 / C_ij_1}
		}
	Output = c(E_V_ij, E_V_ij_1)
	
	return(Output)
	}

## This function calculates the derivative of E(Var(Y_ij|Y_ij-1)) with respect to beta
## This function was written by Shaun Bender

DE_V = function(i, j, lam_ij, x_i, time_i, alpha, r)
	{
	if(Dist == "Poisson"){Output = Dlam_ij(j, lam_ij, x_i)}
	if(Dist == "Negative-Binomial")
		{
		if(j == 1){Output = (2 * r * lam_ij +1) * Dlam_ij(j, lam_ij, x_i)}
		if(j > 1)
			{
			Corr_ij = FindCorr(j, time_i, alpha)[1]
			C_ij = 1 - r * Corr_ij^2 / (1 - Corr_ij^2)
			Output = (2 * r * lam_ij + 1) * Dlam_ij(j, lam_ij, x_i) / C_ij
			}
		}
	if(Dist == "Binomial")
		{
		N = binomN[which(id==i)][j]
		if(j == 1){Output = (K - 2 * lam_ij) * Dlam_ij(j, lam_ij, x_i) / N}
		if(j > 1)
			{
			Corr_ij = FindCorr(j, time_i, alpha)[1]
			C_ij = 1 + Corr_ij^2 / (1 - Corr_ij^2) / N
			Output = (N - 2 * lam_ij) * Dlam_ij(j, lam_ij, x_i) / N / C_ij
			}
		}
	return(Output)
	}

## This function calculates the derivative of Lamda_i1 with respect to beta.
## This function was written by Shaun Bender

Dlam_ij = function(j, lam_ij, x_i)
	{
	if(Dist == "Poisson")
		{
		Dlam = x_i[j,] * lam_ij
		}
	if(Dist == "Negative-Binomial")
		{
		Dlam = x_i[j,] * lam_ij
		}
	if(Dist == "Binomial")
		{
		Dlam = x_i[j,] * lam_ij / (1 + exp(t(beta)%*%x_i[j,]))
		}

	return(Dlam)
	}

## This function calculates the derivative of the Link Function evaluated at "Eval".
## This function was written by Shaun Bender

DerivG = function(i, Eval, r)
	{
	if(Dist == "Poisson")
		{
		Output = 1 / Eval
		}
	if(Dist == "Negative-Binomial")
		{
		Output = 1 / Eval / (1 + r * Eval)
		}
	if(Dist == "Binomial")
		{
		N = binomN[which(id==i)][j]
		Output = N / Eval / (N - Eval)
		}
	
	return(Output)
	}

####################################################################
## Compilation of results
####################################################################

## This function organizes and output, computes several statistics of
## interest.
## This function was written by Shaun Bender, based on code by Victoria 
## Gamerman

CompileResults = function(model, formula, N_Alp, N_Var, N_Subjects)
	{
	mle.alpha = model$par[1:N_Alp]
	mle.beta = model$par[(N_Alp+1):(N_Alp+N_Var)]
	if(Dist == "Negative-Binomial"){mle.r = model$par[N_Alp+N_Var+1]}
	mle.full = -model$value
	mle.cov = solve(model$hessian)

	AIC = 2*(N_Var+1)-2*(mle.full) 
	BIC = log(N_Subjects)*(length(mle.beta)+1)-2*(mle.full)
	
	Stderr = matrix(NA, nrow = N_Var, ncol = 1)
	Wald = matrix(NA, nrow = N_Var, ncol = 1)
	pval = matrix(NA, nrow = N_Var, ncol = 1)
	for(p in (N_Alp+1):(N_Alp+N_Var))
		{
		Stderr[p-N_Alp,] = sqrt(mle.cov[(p),(p)])
		Wald[p-N_Alp,] = (mle.beta[p-N_Alp] / sqrt(mle.cov[(p),(p)]))^2
		pval[p-N_Alp,] = 1-pchisq(Wald[p-N_Alp,1], df = 1, lower.tail = TRUE, log.p = FALSE)
		}
	results = cbind(mle.beta, Stderr, Wald, pval)
	Alpha_Cov = NULL
	for(p in 1:N_Alp)
		{
		Alpha_Cov = c(Alpha_Cov, sqrt(mle.cov[p,p]))
		}
	alpha_results = cbind(mle.alpha,Alpha_Cov)
	fit_stats = rbind(mle.full, AIC, BIC)

	if(Dist == "Negative-Binomial")
		{
		r_Cov = NULL
		r_Cov = sqrt(mle.cov[N_Alp+N_Var+1,N_Alp+N_Var+1])
		r_results = cbind(mle.r,r_Cov)
		colnames(r_results) = c("Estimate", "Std.err")
		rownames(r_results) = "r"
		}

	#format output
	rownames(fit_stats) = c("Log-Likelihood:", "AIC:", "BIC:")
	colnames(fit_stats) = c("")
	colnames(results) = c("Estimate", "Std.err", "Wald", "Pr(>|W|)")
	
	rownames(results) = c("(Intercept)", all.vars(formula[[3]]))

	colnames(alpha_results) = c("Estimate", "Std.err")
	rownames(alpha_results) = c(rep("alpha",N_Alp)) 

	if(Dist != "Negative-Binomial"){return(list(fit_stats, results, alpha_results))}
	if(Dist == "Negative-Binomial"){return(list(fit_stats, results, alpha_results, r_results))}
	}

####################################################################
## Constraints function
####################################################################

Constraints = function(Values)
	{
	Last = length(Values)
	Output = rep(NA, 1)
	if(CorrStr == "AR(1)")
		{
		Output[1] = 1 - Values[1]
		Output[2] = Values[1] + 1
		if(Dist == "Negative-Binomial"){Output[3] = - Values[1]^2 + 1/(Values[Last]+1)}
		}
	if(CorrStr == "Markov")
		{
		Output[1] = 1-Values[1]
		Output[2] = Values[1] + 1
		if(Dist == "Negative-Binomial"){Output[3] = -Values[1]^2 + 1/(Values[Last]+1)}
		}
	if(CorrStr == "AD(1)")
		{
		for(i in 0:(N_Alp-1))
			{
			Output[2*i+1] = 1 - Values[i+1]
			Output[2*i+2] = Values[i+1] + 1
			}
		if(Dist == "Negative-Binomial")
			{
			for(i in 1:N_Alp)
				{
				Output[2*N_Alp+i] = - Values[i]^2 + 1/(Values[Last]+1)
				}
			}
		}
	return(Output)
	}

ConstraintsJacobian = function(Values)
	{
	Last = length(Values)
	if(CorrStr == "AR(1)")
		{
		if(Dist != "Negative-Binomial"){Output = matrix(0, 2, Last)}
		if(Dist == "Negative-Binomial")
			{
			Output = matrix(0, 3, length(Values))
			Output[3,1] = -2*Values[1]
			Output[3,Last] = 1 / (Values[Last] + 1)^2
			}
		Output[1,1] = -1
		Output[2,1] = 1
		}
	if(CorrStr == "Markov")
		{
		if(Dist != "Negative-Binomial"){Output = matrix(0, 2, Last)}
		if(Dist == "Negative-Binomial")
			{
			Output = matrix(0, 3, length(Values))
			Output[3,1] = -2*Values[1]
			Output[3,Last] = 1 / (Values[Last] + 1)^2
			}
		Output[1,1] = -1
		Output[2,1] = 1
		}
	if(CorrStr == "AD(1)")
		{
		if(Dist != "Negative-Binomial"){Output = matrix(0, 2*N_Alp, Last)}
		if(Dist == "Negative-Binomial")
			{
			Output = matrix(0, 3*N_Alp, Last)
			for(i in 1:N_Alp)
				{
				Output[2*N_Alp+i,i] = -2*Values[i]
				Output[2*N_Alp+i,Last] = 1 / (Values[Last] + 1)^2
				}
			}
		for(i in 0:(N_Alp-1))
			{
			Output[2*i+1,i+1] = -1
			Output[2*i+2,i+1] = 1
			}
		}
	return(Output)
	}

####################################################################
## Analysis function to call when running code
####################################################################

EndResults = function(formula, CorrInput, Dist, DatasetInput, IDInput, TimeInput, start.values, DistOpts)
	{
	id = IDInput
	t = TimeInput
	d = dim(DatasetInput)
	k <<- length(all.vars(formula))-1
	dt.fm = data.frame(DatasetInput)

	Dist <<- Dist

	if(Dist == "Binomial")
		{
		dataset <<- data.proc(data = dt.fm, formula = formula, time = t, id = id, del.n = 0, binom = DistOpts)
		}
	if(Dist != "Binomial")
		{
		dataset <<- data.proc(data = dt.fm, formula = formula, time = t, id = id, del.n = 0, binom = NULL)
		}
  	m <<- dataset$m
  	n <<- dataset$n
  	id <<- dataset$id
  	.GlobalEnv$time <- dataset$time
  	autotime <<- dataset$autotime
	binomN <<- dataset$binomN

	CorrStr <<- CorrInput
	
	#if(missing(DistOpts)){K = NULL}
	#if(!missing(DistOpts))
	#	{
	#	#K = DistOpts[order(IDInput,TimeInput)]
	#	K = DistOpts
	#	}
	#K <<- K

	N_Var = length(all.vars(formula[[3]]))+1

	if(CorrStr == "AR(1)" | CorrStr == "Markov"){N_Alp <<- 1}
	if(CorrStr == "AD(1)"){N_Alp <<- length(unique(time))-1}
	
	lb = rep(-Inf, length(start.values))
	ub = rep(Inf, length(start.values))
	lb[1:N_Alp] = rep(-1, N_Alp)
	ub[1:N_Alp] = rep(1, N_Alp)
	full.ml = auglag(par = start.values, fn = drv.logl, hin = Constraints, hin.jac = ConstraintsJacobian,
			control.outer = list("itmax" = 1000, "trace" = FALSE, ilack.max = 10, 
			eps = 10^-8))
	
	Output = CompileResults(full.ml, formula, N_Alp, N_Var, m)
	
	rm(list = c('k','dataset','m','n','time','autotime','CorrStr','K'), pos = ".GlobalEnv")

	return(Output)
	}
