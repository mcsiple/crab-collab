data {
int N;
vector[N] x;
vector[N] y;
}
parameters {
vector[2] log_a;
ordered[2] log_b;
real<lower=0> sigma;
}
transformed parameters {
vector<lower=0>[2] a;
vector<lower=0>[2] b;
a = exp(log_a);
b = exp(log_b);
}
model {
vector[N] ypred;
ypred = a[1]*exp(-b[1]*x) + a[2]*exp(-b[2]*x);
y ~ lognormal(log(ypred), sigma);
log_a ~ normal(0, 1); // I’ve hardcoded the prior mean and sd of log_a and log_b.
log_b ~ normal(0, 1); // If model is well-scaled, these priors are weakly informative.
}
