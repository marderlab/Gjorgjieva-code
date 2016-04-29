result = 0;
for ii=1:50
    V1=dlmread(sprintf('V1_mean_0_sig_0_num_%d.dat',ii));
    V2=dlmread(sprintf('V2_mean_0_sig_0_num_%d.dat',ii));
    figure;
    plot(0:0.1:5000,V1,'b');
    hold on;plot(0:0.1:5000,V2,'r--');
    xlim([4500 5000]);
    prompt = 'Enter 0 or 1 if complete patern?';
    x = input(prompt);
    result = result+x;
    close
end
result
%%
result = 0;
for ii=1:50
    V1=dlmread(sprintf('V1_mean_5_sig_0_num_%d.dat',ii));
    V2=dlmread(sprintf('V2_mean_5_sig_0_num_%d.dat',ii));
    figure;
    plot(0:0.1:5000,V1,'b');
    hold on;plot(0:0.1:5000,V2,'r--');
    xlim([4500 5000]);
    prompt = 'Enter 0 or 1 if complete patern?';
    x = input(prompt);
    result = result+x;
    close
end
result