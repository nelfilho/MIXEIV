stan_codeSN <- "
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
vector[2] Delta;
real<lower=0, upper=1> peso;
real x[N];
real<lower=0> t[N];
}
model {
beta~normal(0,100);
alpha~normal(0,100);
omega2_0 ~ inv_gamma(0.01,0.01);
omega2_1 ~ inv_gamma(0.01,0.01);
mu~normal(0,100);
gamma2 ~ inv_gamma(0.01,0.01);
Delta ~ normal(0,100);
peso ~ beta(1, 1);

for(i in 1:N){
t[i] ~ normal(0,1)T[0,];
target += log_mix(peso, normal_lpdf(x[i] | mu[1] + Delta[1]*t[i], sqrt(gamma2[1])),normal_lpdf(x[i] | mu[2] + Delta[2]*t[i], sqrt(gamma2[2])));
xobs[i] ~ normal(x[i],sqrt(omega2_0));
y[i] ~ normal(alpha+beta*x[i], sqrt(omega2_1));
}
}
"