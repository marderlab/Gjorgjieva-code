function tauhNa = tau_hNa(V)

tauhNa = (0.67./(1+exp((V+62.9)./-10.0))).*(1.5 + 1./(1+exp((V+34.9)./3.6)));

end