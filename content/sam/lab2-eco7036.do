clear all
set more off
*cd "C:\Users\Samuel\Google Drive\ecole\université\travail\eco7036\lab 2"

capture log close
*log using "log\lab2.log",replace

*use ""

set obs 1000

g beta=.
g rejet=.


mata
j=1000 //Nombre de réplication

/*Initiatialisation de nos vecteurs de résultats*/
beta = J(j,1,0)
rejet = J(j,1,0)

i=1
while (i<=j){
	// Paramètres du modèle
	n = 10
	beta0 = 0.8
	beta1 = 1
	sigma2_b = 4
	sigma2_e = 2

	meanx = J(n,1,4)
	alpha = J(n,1,1)
	
	//Construction du modèle
	x1=sqrt(sigma2_b)*invnormal(uniform(n,1))
	x=meanx+x1
	er=invnormal(uniform(n,1))
	y=beta0*alpha+beta1*x+sqrt(sigma2_e)*er

	//Estimation du modèle
	X = (alpha,x)
	k = cols(X)

	B_OLS = invsym(X'*X)*(X'*y)
	be = B_OLS[2]
	beta[i]=be
		
		
	er1 = y - X*B_OLS
	sig2 =(er1'*er1)/(n-k)

	/*Variance de alpha et de beta*/
	varb = sig2*luinv(X'*X)
	/*Écart-type de alpha et de beta*/
	std = diagonal(sqrt(varb))
	be_std = std[2]

	t = (be)/be_std
	p = 2*(1-t(n-k,abs(t)))
	
	if (p < 0.05) {
		rejet[i]=1
	}
	i++
	}


st_store(.,"beta",beta)
st_store(.,"rejet",rejet)

end

sum beta rejet

hist(beta)

*log close
