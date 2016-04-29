function c = get_color(m,N)

if (m <= N/3)
    c = [m/(N/3) 0 0];
elseif (m <= 2*N/3)
    c = [1 (m-N/3)/(N/3) 0];
else
    c = [1 1 (m-2*N/3)/(N/3)];
end