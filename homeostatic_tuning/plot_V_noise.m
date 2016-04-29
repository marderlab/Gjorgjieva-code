% mean = 0;
% sig = 1;
% tau = 500;
% coef = 0;
% result = 0;
% for ii=1:10
%     V1=dlmread(sprintf('V1_mean_%g_sig_%g_tau_%g_coef_%g_num_%d.dat',mean,sig,tau,coef,ii));
%     V2=dlmread(sprintf('V2_mean_%g_sig_%g_tau_%g_coef_%g_num_%d.dat',mean,sig,tau,coef,ii));
%     figure;
%     plot(0:0.1:5000,V1,'b');
%     hold on;plot(0:0.1:5000,V2,'r--');
%     xlim([4500 5000]);
%     prompt = 'Enter 0 or 1 if complete patern?';
%     x = input(prompt);
%     result = result+x;
%     close
% end
% result

%%
mean = 5;
sig = 5;
tau = 100;
coef = 0.5;
result = 0;
Ntot = 50;
for ii=1:Ntot
    V1=dlmread(sprintf('V1_mean_%g_sig(reset)_%g_tau_%g_coef_%g_num_%d.dat',mean,sig,tau,coef,ii));
    V2=dlmread(sprintf('V2_mean_%g_sig(reset)_%g_tau_%g_coef_%g_num_%d.dat',mean,sig,tau,coef,ii));
    figure;
    plot(0:0.1:5000,V1,'b');
    hold on;plot(0:0.1:5000,V2,'r--');
    xlim([4500 5000]);
    title(sprintf('Run %g',ii));
    prompt = 'Enter 0 or 1 if complete patern?';
    x = input(prompt);
    result = result+x;
    close
end
result
%%
mean = 0;
sig = 5;
tau = 500;
coef = 0.5;
result = 0;
Ntot = 50;
for ii=1:50
    V1=dlmread(sprintf('V1_mean_%g_sig(reset)_%g_tau_%g_coef_%g_num_%d.dat',mean,sig,tau,coef,ii));
    V2=dlmread(sprintf('V2_mean_%g_sig(reset)_%g_tau_%g_coef_%g_num_%d.dat',mean,sig,tau,coef,ii));
    figure;
    plot(0:0.1:5000,V1,'b');
    hold on;plot(0:0.1:5000,V2,'r--');
    xlim([4400 5000]);
    prompt = 'Enter 0 or 1 if complete patern?';
    x = input(prompt);
    result = result+x;
    close
end
result
