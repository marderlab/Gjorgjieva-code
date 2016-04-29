mu = 1; %mean
s_c(1) = 0; %IC
s_c1(1) = 0;
s_c2(1) = 0;
dt = 0.01;
tau_c = 1; %0 is white noise
sigma = 1;
N = 10000000;
for ii=2:N
    rand_g1 = randn;
    ds = (mu-s_c(ii-1))*(1-exp(-dt/tau_c)) + sigma*sqrt(1-exp(-2*dt/tau_c))*rand_g1; 
    s_c(ii) = s_c(ii-1) + ds;
    ds1 = (0-s_c1(ii-1))*(1-exp(-dt/tau_c)) + 1*sqrt(1-exp(-2*dt/tau_c))*rand_g1; 
    s_c1(ii) = s_c1(ii-1) + ds1;
    ds2 = (0-s_c2(ii-1))*(1-exp(-dt/tau_c)) + sigma*sqrt(1-exp(-2*dt/tau_c))*rand_g1; 
    s_c2(ii) = s_c2(ii-1) + ds2;
end


s_c = s_c(0.1*N:N);
s_c1 = s_c1(0.1*N:N);
s_c2 = s_c2(0.1*N:N);
%figure;plot(0.1*N*dt:dt:N*dt,s_c1*sigma+mu)
%hold on;plot(0.1*N*dt:dt:N*dt,s_c,'o')

std(s_c)
%std(s_c1*sigma+mu)
%std(s_c2)

%figure;plot(dt:dt:N*dt,s_c2./s_c)