Acsrinds = Acsr ~= 0;
Ainds = A~= 0;
Adiff = Ainds - Acsrinds;
[x,y] = size(Adiff);
rindicies = zeros(sum(sum(Adiff)),2);
cindicies = zeros(sum(sum(Adiff)),2);

ii = find(Adiff);

for i = 1:length(ii)
    iit = ii(i);
    iir = mod(iit,x);
    iic = (iit-iir)/x;
    rindr = mod(iir,64);
    rindc = (iir-rindr/64);
    cindr = mod(iic,64);
    cindc = (iic-cindr/64);
    rindicies(i,1) = rindr;
    rindicies(i,2) = rindc;
    cindicies(i,1) = cindr;
    cindicies(i,2) = cindc;
end
    
    
    
    
    
    
    
    
    