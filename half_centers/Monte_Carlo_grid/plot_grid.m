function plot_grid(A,g_min_max)
%plots the entries in the matrix A as a grid of points
%here A is a matrix of conductances

figure;
N = size(A,2);
for ii=1:N
    for jj=ii+1:N
        subplot(N-1,N-1,(ii-1)*N+jj - ii);plot(A(:,ii),A(:,jj),'k.');
        xlim([g_min_max(ii,1) g_min_max(ii,2)])
        ylim([g_min_max(jj,1) g_min_max(jj,2)])
    end
end
subplot(6,6,1);xlabel('g_{Na}');
subplot(6,6,1);ylabel('g_{CaT}');
subplot(6,6,8);xlabel('g_{CaT}');
subplot(6,6,8);ylabel('g_{CaS}');
subplot(6,6,15);xlabel('g_{CaS}');
subplot(6,6,15);ylabel('g_{A}');
subplot(6,6,22);xlabel('g_{A}');
subplot(6,6,22);ylabel('g_{KCa}');
subplot(6,6,29);xlabel('g_{KCa}');
subplot(6,6,29);ylabel('g_{Kd}');
subplot(6,6,36);xlabel('g_{Kd}');
subplot(6,6,36);ylabel('g_{H}');
