function tauhCaS = tau_hCaS(V)

tauhCaS = 60 + (150./((exp((V+55)./9))+(exp((V+65)./-16))));

end