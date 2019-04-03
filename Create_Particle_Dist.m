function [] = Create_Particle_Dist(d_mean, sigma, dmin, dmax, width, flShow)

d_mean = log(d_mean);
sigma = (log(sigma));
yp = @(x) 0.5 + 0.5*erf((log(x)-d_mean)/sigma/sqrt(2));

dists = [];
sizes = [];
dcur = dmin;
i = 1;
while(dcur <= dmax)
    dists(i) = yp(dcur+width/2)-yp(dcur-width/2);
    sizes(i) = dcur;
    dcur = dcur + width;
    i = i+1;
end

tot_dist = sum(dists);
dists = dists + (1 - tot_dist)/length(dists);

out = zeros(length(sizes)*2,1);
% sizes = sizes*1e-9;
for i = 1:length(sizes)
    out((i-1)*2+1) = sizes(i)*1e-9;
    out(i*2) = dists(i);
end
% fprintf('Cumulative Sum = %g\n', sum(dists));
save2file(out, 'load\\pdist.txt');
if(flShow)
    return;
end
figure
plot(sizes, dists);