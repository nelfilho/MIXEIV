stan_codeCN <- "
data {
int<lower=0> N;
real xobs[N];
real y[N]; //y é um conjunto de n vetores de tamanho k
}
parameters {
real beta;
real alpha;
real<lower=0> omega2_1;
real<lower=0> omega2_0;
vector[2] mu;
real<lower=0> gamma2[2];
real<lower=0> f;
real<lower=0, upper=1> peso;
real<lower=0, upper=1> tau;
real<lower=0, upper=1> rho;
real x[N];
real<lower=0> uni;
}
model {
beta~normal(0,100);
alpha~normal(0,100);
omega2_0 ~ inv_gamma(0.01,0.01);
omega2_1 ~ inv_gamma(0.01,0.01);
mu~normal(0,100);
gamma2 ~ inv_gamma(2,f);
f ~ gamma(0.1,0.1);  
peso ~ beta(1, 1);
tau ~ beta(1, 1);
rho ~ beta(1, 1);

for(i in 1:N){
uni ~ uniform(0,1);
if(uni < rho){
target += log_mix(peso, normal_lpdf(x[i] | mu[1], sqrt(gamma2[1]/tau)),normal_lpdf(x[i] | mu[2], sqrt(gamma2[2]/tau)));
}else{
target += log_mix(peso, normal_lpdf(x[i] | mu[1], sqrt(gamma2[1])),normal_lpdf(x[i] | mu[2], sqrt(gamma2[2])));
}
xobs[i] ~ normal(x[i],sqrt(omega2_0));
y[i] ~ normal(alpha+beta*x[i], sqrt(omega2_1));
}
}
"
