%create new colormap

num_points = 30;

%red to white to blue
%[1 0 0] -> [1 1 1] -> [0 0 1]
map = zeros(num_points+1,3);
for jj = 1:num_points+1
    ii = jj-1;
    x1 = [1*(ii<=num_points/2)+(1+(num_points/2-ii)/num_points*2)*(ii>num_points/2)];
    x2 = [ii/num_points*2*(ii<=num_points/2)+(1+(num_points/2-ii)/num_points*2)*(ii>num_points/2)];
    x3 = [ii/num_points*2*(ii<=num_points/2)+1*(ii>num_points/2)];
    map(jj,:) = [x1 x2 x3];
end