function P = Prix(sigma1,sigma2,rho,alphaS10,betaS20,T)
sigma = sqrt(sigma1*sigma1+sigma2*sigma2-2*sigma1*sigma2*rho);
d1 = (log(betaS20/alphaS10)+sigma*sigma*T/2)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
P = alphaS10*normcdf(-d2)-betaS20*normcdf(-d1);
end