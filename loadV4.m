function [v4] = loadV4(name)

v = load(name);
lx = v(:,1:4:end);
ly = v(:,2:4:end);
lz = v(:,3:4:end);
lw = v(:,4:4:end);
[x,y] = size(lx);
v4 = zeros(x,y,4);
v4(:,:,1) = lx;
v4(:,:,2) = ly;
v4(:,:,3) = lz;
v4(:,:,4) = lw;




