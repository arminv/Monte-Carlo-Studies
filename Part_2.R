# Armin Varshokar

# ======================================================================================================
# 1) Compare the following methods of optimization: 1) Newton-Raphson ; 2) Basic MC (Uniform) ; 3) MC (other than unifrom) ; 4) Simulated Annealing:

# We choose this function : f(x) =  abs(x*exp(sin(x)^2)  and aim to find its minimum location on the interval : [-15 , 15] 

f1		<- function(x){ abs(x*exp(sin(x)^2)) }

curve(f1, xlim= c(-15,15), main= " Min of f(x) =  abs(x*exp(sin(x)^2) on [-15,15] ")		# plotting the chosen function	
abline(v=0, col="red")											# this is the true value of the global minimum on interval [-15,15] which is found analytically and is equal to zero

# it seems like a good candidate for demonstrating the advantages/disadvantages of each of the optimization methods because it has many local minima and one global minimum on the interval [-15, 15].

# We can simply find the minimum on a given interval using built-in functions: Newton-Raphson, nlm() and optimize() :

optimize(function(x) abs(x*exp(sin(x)^2)), interval=c(-15,15), maximum=FALSE)	# returns : argmin = min = 1.776357e-15 which is reasonable

#==========================================================
# ------> METHODS: 

# 1) Newton-Raphson: due to complexity of this candidate function we use the built-in nlm() function to implement the method:

nlm(function(x) abs(x*exp(sin(x)^2)), p = -1)		# returns : argmin = 1.541264e-07  & min = -1.541264e-07

# This seems fine as it is very close to zero.

nlm(function(x) abs(x*exp(sin(x)^2)), p = -5)		# returns : argmin = 9.398114 & min = 	-9.371327 ; ** NOTE: it gets stuck in a local minimum due to a far starting point relative to the global minimum on [-15,15]

# Most numerical optimization methods (including Newton-Raphson) do not always do a good job of finding the extrema. This is especially true if the function has many local extrema and the method gets stuck in them. This is the case here and therefore we need better methods for this function.

#===========================
# 2) The Basic MC: we implement this as a function so that we can use it for simulation studies later on:

Basic_1	<- function(m){
	u			<- runif(m, -15, 15)									# Unifrom variables on [-15, 15] generated 
	temp		<- sapply(u, f1)										# feeding the uniform RVs into the function f1(x)
	loc		<- which.min(temp)								# which value in temp is the smallest?
	return(u[loc])														# which value in u does it correspond to? This is our answer hence return it
} 

set.seed(12345)														# making it reproducible
m		<- 10^5
Basic_1(m)																# returns : 0.0003044028 which seems to be fairly acceptable. We will check this later by doing a simulation study

#===========================
# 3) MC Other Than Uniform: 


# A) We propose density : g(x) = x^2 

curve(f1, xlim=c(-15, 15), main= " A:  g(x) = x^2 ")
curve(x^2, xlim=c(-15, 15), col="red", add=T)
abline(v=0, h=f1(0), col="green")

# Based on the graph this density seems like a reasonable choice as it has some desirable similarities with our function. In other words we need a distribution that explores more in areas where the function f(x) is large. Therefore we proceed with this and will check its performance with a simulation study later.

# Similar to Basic MC we implement this as a function. Note also that we find the CDF to be G(x) = x^3/3 and the inverse CDF to be G_inv(x) = (3*x)^(1/3). We use Inverse Transform Method to sample from this distribution:

MC_1 		<- function(m){
	u			<- runif(m, 0, 1)										# uniform RVs
	inv		<- (3*u)^(1/3)											# sampling from g(x) by Inverse Transform method
	temp		<- sapply(inv, f1)										# feeding the RVs from density g(x) into the function f1(x)
	loc		<- which.min(temp)								# which value in temp is the smallest?
	return(u[loc])														# which value in u does it correspond to? This is our answer hence return it
}

set.seed(12345)														# making it reproducible
m		<- 10^5															# number of samples
MC_1(m)																# returns : 4.265457e-06 which seems reasonable and more accurate than the basic MC. We will have to check this by using a simulation study.

#===========
# B) We propose another density : g(x) = x^2 / 6 

curve(f1, xlim=c(-15, 15), main= " B:  g(x) = x^2 / 6 ")
curve(x^2/6, xlim=c(-15, 15), col="red", add=T)
abline(v=0, h=f1(0), col="green")

# Based on the graph this density seems like a more reasonable choice than A. We will check whether it performs better than A later in the simulation study.

# Note that we find the CDF to be G(x) = x^3 / 18 and the inverse CDF to be G_inv(x) = (18*x)^(1/3). We use Inverse Transform Method to sample from this distribution:

MC_2 		<- function(m){
	u			<- runif(m, 0, 1)										# uniform RVs
	inv		<- (18*u)^(1/3)										# sampling from g(x) by Inverse Transform method
	temp		<- sapply(inv, f1)										# feeding the RVs from density g(x) into the function f1(x)
	loc		<- which.min(temp)								# which value in temp is the smallest?
	return(u[loc])														# which value in u does it correspond to? This is our answer hence return it
}

set.seed(123456)														# making it reproducible
m		<- 10^5
MC_2(m)																# returns :  4.192349e-06 which seems reasonable 

#===========
# C) We propose another density : g(x) = x^4 / 200 

curve(f1, xlim=c(-15, 15), main= " C:  g(x) = x^4 / 200 ")
curve(x^4/200, xlim=c(-15, 15), col="red", add=T)
abline(v=0, h=f1(0), col="green")

# Based on the graph this density seems like a more reasonable choice than A. We will check whether it performs better than A later in the simulation study.

# Note that we find the CDF to be G(x) = x^5 / 1000 and the inverse CDF to be G_inv(x) = (1000*x)^(1/5). We use Inverse Transform Method to sample from this distribution:

MC_3 		<- function(m){
	u			<- runif(m, 0, 1)										# uniform RVs
	inv		<-	 (1000*u)^(1/5)									# sampling from g(x) by Inverse Transform method
	temp		<- sapply(inv, f1)										# feeding the RVs from density g(x) into the function f1(x)
	loc		<- which.min(temp)								# which value in temp is the smallest?
	return(u[loc])														# which value in u does it correspond to? This is our answer hence return it
}

set.seed(1234567)													# making it reproducible
m		<- 10^5
MC_3(m)																# returns : 1.689233e-05 which seems reasonable 

#===========================
# 4) Simulated Annealing : we implement this method as a function so that we can use it for simulation studies later on:

# First, in order to demonstrate the method, we use a temperature schedule and a scale that is likely to work and we also use the Uniform distribution as the proposal density (similar to basic MC):

set.seed(12345)
tol			<-	1e-20
min.iter	<-	50
max.iter	<-	500
temp			<-	1/ log(1+1:max.iter)^4					 # temperature schedule
scale			<-	sqrt(temp)		

x				<-	runif(1)											# this is the starting point, we choose it randomly for the purpose of this demonstration
fval			<-	f1(x)
fcur			<-	fval

iter			<-	1
diff			<-	1
Counter 	<- 0														# keeping track of the number of times we reject each time

for(i in 1:max.iter){
	
	prop		<-	x[i] + scale[i]*(runif(1, -1, 1))			
	rho		<-	exp((fcur-f1(prop))/temp[i])   		# fcur-f1(prop) because we want the minimum
	rho		<-	min(rho,1)
	u			<-	runif(1)

	# if REJECT, then reset proposal	
	if((u>rho)||(prop>1)||(prop<0))		prop <-	x[i] ; Counter <- Counter + 1
	x			<-	c(x, prop)
	fcur		<-	f1(prop)
	fval		<-	c(fval, fcur)

	# check STOPPING RULE
	diff		<-	max(fval)-max(fval[1:(i/2)])
	if((iter>min.iter)&&(length(unique(x[(i/2):i]))>1)&&(diff < tol)) break
	
	iter		<-	iter+1
}
	
list(x= x, fval= fval, Counter=Counter)				# respectively: argmin, min, number of rejects before ending

# plotting the results of the optimization:
grid	<-	seq(-1.5,1.5, length.out=1001)
plot(grid, f1(grid), type="l", col=gray(0.7), main=" Simulated Annealing : f(x) =abs(x*exp(sin(x)^2)) " )
points(x[1], fval[1], pch=19, col="red")
lines(x, fval, col="red", type="b")
points(x[length(x)], fval[length(x)], pch=19, col="blue")

#===========
# Now that we have shown the method above works, we proceed to program it as a function and find other variations of it for the simulation study:

# First case:

SA_1	<- function(min.iter=50, max.iter=500, tol=1e-20){
	
	tol			<-	tol
	min.iter	<-	min.iter
	max.iter	<-	max.iter

	temp			<-	1/ log(1+1:max.iter)^4				# Temperature Schedule
	scale				<-	10*sqrt(temp)		

	x			<-	1						# we take e.g. x = 1 for all the annealing methods in order to control for this factor in our simulation study					
	fval		<-	f1(x)
	fcur		<-	fval

	iter			<-	1
	diff			<-	1
	Counter 	<- 0					# keeping track of the number of times we reject each time

	for(i in 1:max.iter){
	
		prop		<-	x[i] + scale[i]*(runif(1, -1, 1))			
		rho		<-	exp((fcur-f1(prop))/temp[i])  			# fcur-f1(prop) because we want the minimum
		rho		<-	min(rho,1)
		u			<-	runif(1)

		# if REJECT, then reset proposal	
		if((u>rho)||(prop>1)||(prop<0))		prop <-	x[i] ; Counter <- Counter + 1
		x			<-	c(x, prop)
		fcur		<-	f1(prop)
		fval		<-	c(fval, fcur)

		# check STOPPING RULE
		diff		<-	max(fval)-max(fval[1:(i/2)])
		if((iter>min.iter)&&(length(unique(x[(i/2):i]))>1)&&(diff < tol)) break
	
		iter		<-	iter+1
}
	return(list(x= x, fval= fval, Counter=Counter))			# return the following: argmin, min, number of rejects before ending
}

SA_1()							# making sure the function works fine

#===========
# Second case: a different temperature schedule:

SA_2	<- function(min.iter=50, max.iter=500, tol=1e-20){
	
	tol			<-	tol
	min.iter	<-	min.iter
	max.iter	<-	max.iter

	temp			<-	10/ log(1+1:max.iter)^3			# Temperature Schedule
	scale				<-	sqrt(temp)		

	x			<-	1						# we take e.g. x = 1 for all the annealing methods in order to control for this factor in our simulation study				
	fval		<-	f1(x)
	fcur		<-	fval

	iter			<-	1
	diff			<-	1
	Counter 	<- 0					# keeping track of the number of times we reject each time

	for(i in 1:max.iter){
	
		prop		<-	x[i] + scale[i]*(runif(1, -1, 1))			
		rho		<-	exp((fcur-f1(prop))/temp[i])  		 # fcur-f1(prop) because we want the minimum
		rho		<-	min(rho,1)
		u			<-	runif(1)

		# if REJECT, then reset proposal	
		if((u>rho)||(prop>1)||(prop<0))		prop <-	x[i] ; Counter <- Counter + 1
		x			<-	c(x, prop)
		fcur		<-	f1(prop)
		fval		<-	c(fval, fcur)

		# check STOPPING RULE
		diff		<-	max(fval)-max(fval[1:(i/2)])
		if((iter>min.iter)&&(length(unique(x[(i/2):i]))>1)&&(diff < tol)) break
	
		iter		<-	iter+1
}
	return(list(x= x, fval= fval, Counter=Counter))
}

SA_2()							# making sure the function works fine

#===========
# Third case: yet another different temperature schedule:

SA_3	<- function(min.iter=50, max.iter=500, tol=1e-20){
	
	tol			<-	tol
	min.iter	<-	min.iter
	max.iter	<-	max.iter

	temp			<-	300/ (1+1:max.iter)^6			# Temperature Schedule
	scale				<-	sqrt(temp)		

	x			<-	1						# we take e.g. x = 1 for all the annealing methods in order to control for this factor in our simulation study				
	fval		<-	f1(x)
	fcur		<-	fval

	iter			<-	1
	diff			<-	1
	Counter 	<- 0					# keeping track of the number of times we reject each time

	for(i in 1:max.iter){
	
		prop		<-	x[i] + scale[i]*(runif(1, -1, 1))			
		rho		<-	exp((fcur-f1(prop))/temp[i])  		 # fcur-f1(prop) because we want the minimum
		rho		<-	min(rho,1)
		u			<-	runif(1)

		# if REJECT, then reset proposal	
		if((u>rho)||(prop>1)||(prop<0))		prop <-	x[i] ; Counter <- Counter + 1
		x			<-	c(x, prop)
		fcur		<-	f1(prop)
		fval		<-	c(fval, fcur)

		# check STOPPING RULE
		diff		<-	max(fval)-max(fval[1:(i/2)])
		if((iter>min.iter)&&(length(unique(x[(i/2):i]))>1)&&(diff < tol)) break
	
		iter		<-	iter+1
}
	return(list(x= x, fval= fval, Counter=Counter))
}

SA_3()						# making sure the function works fine


#==========================================================
# ------> COMPARISON:

# 1) Fixed Sample size: we can now compare all the methods (Note: The simulation may take a while to finish):

set.seed(123456)
m	<- 1e5												# sample size
B <- 10^3											# number of simulations / replications

res_Newt_1 <- rep(0, B)						# containers for results of: Newton-Raphson with different starting points:
res_Newt_2 <- rep(0, B)
res_Newt_3 <- rep(0, B)

res_Basic_1 <- rep(0, B)						# container for result of: Basic MC 

res_MC_1 	<- rep(0, B)						# containers for results of: Monte Carlo with different proposal densities:
res_MC_2 	<- rep(0, B)
res_MC_3	<- rep(0, B)

res_SA_1 <- rep(0, B)							# containers for results of: Simulated Annealing with different temperature schedules:
res_SA_2 <- rep(0, B)
res_SA_3 <- rep(0, B)

for (i in 1:B){
# these are deterministic therefore they do not change, but we leave them here for the purpose of demonstration:
	res_Newt_1[i] <- 	nlm(f1, p = -1)	$minimum		
	res_Newt_2[i] <- 	nlm(f1, p = 2)$minimum
	res_Newt_3[i] <- 	nlm(f1, p = 4)$minimum

# Uniform proposal density:
	res_Basic_1[i] <- 	Basic_1(m)					 

# Different: proposal densities
	res_MC_1[i] 	<- 	MC_1(m)					
	res_MC_2[i] 	<- 	MC_2(m)
	res_MC_3[i]	<- 	MC_3(m)

# Same: starting points (x=1), tolerances, min/max no of iterations. Different: temperature schedules
	res_SA_1[i] 	<- 	SA_1(min.iter=50, max.iter=500, tol=1e-20)	$x[length(x)]	
	res_SA_2[i] 	<- 	SA_2(min.iter=50, max.iter=500, tol=1e-20)$x[length(x)]
	res_SA_3[i] 	<- 	SA_3(min.iter=50, max.iter=500, tol=1e-20)$x[length(x)]
	
}

# producing the comparative boxplots in one figure:

boxplot(res_Newt_1, res_Newt_2, res_Newt_3, res_Basic_1, res_MC_1, res_MC_2, res_MC_3, res_SA_1, res_SA_2, res_SA_3, main="Comparison of All Optimization Methods", names=c("Newt_1", "Newt_2","Newt_3","Basic_1","MC_1","MC_2","MC_3","SA_1","SA_2","SA_3"), col="grey")
abline(h= 0, col="red")			# true value of global argmin on [-15,15]

# due to large differences in values, we should take simulated annealing results and plot them together (i.e. to 'zoom in'):

boxplot(res_SA_1, res_SA_2, res_SA_3, main="Comparison of Simulated Annealing Methods", names=c("SA_1","SA_2","SA_3"), col="grey", ylim=c(0,1))
abline(h= 0, col="red")

# similarly for the rest but we leave out the results for Newt_2 as it gets stuck in local minima and produces outlier results:

boxplot(res_Newt_1, res_Newt_3, res_Basic_1, res_MC_1, res_MC_2, res_MC_3, main="Comparison of MC & Newton-Raphson Methods", names=c("Newt_1", "Newt_3","Basic_1","MC_1","MC_2","MC_3"), col="grey", range=0)
abline(h= 0, col="red")

# we leave out Newt_1 and Newt_2 and Basic_1 results out this time:

boxplot(res_MC_1, res_MC_2, res_MC_3, main="Comparison of MC Methods", names=c("MC_1","MC_2","MC_3"), col="grey", range=0)
abline(h= 0, col="red")

# Conclusion: please see the report.

#===========
# 2) Increasing Sample size: we now increase the sample size with each iteration for MC optimization methods in order to see how they evolve:

set.seed(123456)
m1	<- 1e3												
m2	<- 10^4										

res2_Basic_1 <- rep(0, m2-m1)					# container for result of: Basic MC 

res2_MC_1 	<- rep(0, m2-m1)						# containers for results of: Monte Carlo with different proposal densities:
res2_MC_2 	<- rep(0, m2-m1)
res2_MC_3	<- rep(0, m2-m1)

for (i in (m1+1):m2){
	
# Uniform proposal density:
	res2_Basic_1[i-m1] <- 	Basic_1(i)					 

# Different: proposal densities
	res2_MC_1[i-m1] 	<- 	MC_1(i)					
	res2_MC_2[i-m1] 	<- 	MC_2(i)
	res2_MC_3[i-m1]	<- 	MC_3(i)
	
}

par(mfrow=c(4,1))
plot((m1+1):m2, res2_Basic_1, type= "l", col="red", ylim=c(-0.005, 0.005), main="Comparison of MC Methods, Increasing Sample Size", ylab="Basic_1", xlab="")
abline(h=0)
plot((m1+1):m2, res2_MC_1, type= "l", col="blue", ylim=c(0, 0.002), ylab="MC_1",xlab="")
abline(h=0, col="red")
plot((m1+1):m2, res2_MC_2, type= "l", col="green", ylim=c(0, 0.002), ylab="MC_2", xlab="")
abline(h=0, col="red")
plot((m1+1):m2, res2_MC_3, type= "l", col="orange", ylim=c(0, 0.002), ylab="MC_3")
abline(h=0, col="red")

# Conclusion: the performance of the basic MC does not seem to be as sensitive to sample size as other MC methods.

# ======================================================================================================
# ======================================================================================================