function mKCainf = m_KCainf(V,Ca)

mKCainf = (Ca./(Ca+3)).*(1./(1+exp((V+28.3)./-12.6)));

end