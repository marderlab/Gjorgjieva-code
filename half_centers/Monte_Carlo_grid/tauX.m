function tauX = tauX(V,A,B,D,E)

tauX = A - B./(1+exp((V+D)./E));

end