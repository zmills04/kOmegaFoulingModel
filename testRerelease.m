load releaseDebug.txt;

% rd = releaseDebug(1:end/2);
rd = releaseDebug;
Yvals = rd(1:8:end);
Urand = rd(2:8:end);
Ufluid = rd(3:8:end);
Uind = rd(4:8:end);
% Yvals = 1.5+[0.0, 0.25, 0.5, 0.75, 1.5, 2, 199, 199.25, 199.5, 199.75, 200]';

Uvals = [Yvals, Yvals];
Uinds = [Yvals, Yvals];
for i = 1:length(Yvals)
    YindVal = floor(Yvals(i))+1;
    Ulow = lbux(90,YindVal);
    Uhigh = lbux(90,YindVal+1);
    if(Yvals(i) < 2)
        ddY = (Yvals(i) - 1.5) / 0.5;
    elseif(Yvals(i) > 201)
        ddY = (Yvals(i) - 201)/0.5;
    else
        ddY = Yvals(i) - floor(Yvals(i));
    end
    Uvals(i,:) = [Uhigh * ddY + Ulow*(1-ddY),Ufluid(i)];
    
    Uinds(i,:) = [90 + (YindVal-1)*2048, Uind(i)];    
end
Uinds = [Uinds,Uinds(:,2) - Uinds(:,1)];
Uvals = [Uvals,Uvals(:,2) - Uvals(:,1)];

Errs = [];

for i = 1:length(Yvals)
    if((abs(Uvals(i,3))/Uvals(i,1) > 0.001))
        Errs = [Errs;[Yvals(i),Uvals(i,1),Uvals(i,2)]];
    end
end



