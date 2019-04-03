function []=MakeChannel_Test(H, Len, Re, devId)
if(~exist('load','dir'))
    mkdir('load');
end

if(nargin == 3)
  devId = 0;
endif

%  ptypes
dbl_t = 1;
int_t = 2;
str_t = 3;
bool_t = 4;
dblvec_t = 5;
intvec_t = 6;
systemParamNum = 1;
fluidParamNum = 2;
thermalParamNum = 3;
lsParamNum = 5;
trParamNum = 6;
flParamNum = 7;

%docType and final docCells array will be created right before calling yaml_emitter
%docType = [];
%docCells = {};
%docInd = 1;


sysCells = {};
lbCells = {};
fdCells = {};
lsCells = {};
trCells = {};
flCells = {};






% close all
Yeq = H*2;
dX = 1;
Umean = 0.05;
visc = Umean*Yeq/2/Re

flShow = 1;

tau = 1/(visc*3+0.5);

nX = Len;
nY = Yeq+4;

sysCells = appendCell(sysCells, int_t, "nX", nX);
sysCells = appendCell(sysCells, int_t, "nY", nY);
sysCells = appendCell(sysCells, dbl_t, "Pipe Radius", H);

lbCells = appendCell(lbCells, dbl_t, "Tau", tau);
lbCells = appendCell(lbCells, dbl_t, "Re", Re);

lsCells = appendCell(lsCells, dbl_t, "LS Spacing", dX);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


err = 1e-3;
Cx = [];
Cy = [];
BL = [];

flag = 0;
xval = -0.75;
yval = -Yeq/2;
while(flag == 0)
    Cx = [Cx;xval];
    Cy = [Cy;yval];

    if(xval > nX-0.5)
        flag = 1;
    else
        xval = xval+dX;
    end
end

for i = 1:length(Cx)-1
    BL = [BL;[i,i+1]];
end
startind = length(Cx)+1;
yval = Yeq/2;
flag = 0;
while(flag == 0)
    Cx = [Cx;xval];
    Cy = [Cy;yval];
    
    if(xval <= -0.75)
        flag = 1;
    else
        xval = xval-dX;
    end
end

for i = startind:length(Cx)-1
    BL = [BL;[i,i+1]];
end

BN = 1:length(Cx);

Yini = nY/2 - dX + 0.5;

Cy = Cy + Yini;



%for i = 1:length(Cx)
%    xval = Cx(i)+0.5;
%    if(xval < Period)
%        dispp = Amp;
%    elseif(xval > Period*(num_per + 1))
%        dispp = Amp;
%    else
%        dispp = Amp*cos((xval)*2*pi/Period);
%    end
%    Cy(i) = Cy(i)+dispp;
%end
%
%%flag = testdistances_2D(Cx,Cy,BL,dX,1e-3);
%%     flag = 0;
%attempt = 0;
%while(flag == 1)
%      dispy = 1 - 2*rand(1,1);
%      Cy = Cy+dispy;
%      flag = testdistances_2D(Cx,Cy,BL,dX,err);
%      attempt = attempt+1;
%      if(mod(attempt,1000) == 0)
%          err = err/2;
%          fprintf('%d\n',attempt);
%      end
%end
%
%middle = 0.5*(max(Cy) + min(Cy));
%disp = middle - nY/2;
%disp = round(disp);
%Cy = Cy - (disp);





save2file([Cx Cy], 'load\\lsc0.txt');
save2file(BN-1, 'load\\lsbn.txt');
save2file([nX nY], 'load\\lbsize.txt');
save2file(BL-1, 'load\\lsbl.txt');
[Dp_size, Dp_dist] = Create_Particle_Dist_yaml(130, 1.4, 50, 290, 20, flShow);

lsCells = appendCell(lsCells, int_t, "nN", length(Cx));

trCells = appendCell(trCells, bool_t, "Use Par Solver", false);
trCells = appendCell(trCells, int_t, "Num Par Sizes", length(Dp_size));
trCells = appendCell(trCells, dblvec_t, "Par Diameters", Dp_size);
trCells = appendCell(trCells, dblvec_t, "Dp Dists", Dp_dist);

docType = [];
docCells = {};
docInd = 1;

if(length(sysCells) > 0)
  docType(docInd) = systemParamNum;
  docCells{docInd} = sysCells;
  docInd += 1;
else
  exit("no system parameters defined");
end

if(length(lsCells) > 0)
  docType(docInd) = lsParamNum;
  docCells{docInd} = lsCells;
  docInd += 1;
else
  exit("no ls parameters defined");
end

if(length(lbCells) > 0)
  docType(docInd) = fluidParamNum;
  docCells{docInd} = lbCells;
  docInd += 1;
else
  exit("no fluid parameters defined");
end


% remain parameters files are not required depending on the type of simulation
% if they are empty, a dummy key is added to avoid throwing error in yaml-cpp
if(length(fdCells) == 0)
  fdCells = appendCell(fdCells, int_t, "Dummy", 0);
end

if(length(trCells) == 0)
  trCells = appendCell(trCells, int_t, "Dummy", 0);
end

if(length(flCells) == 0)
  flCells = appendCell(flCells, int_t, "Dummy", 0);
end

docType(docInd) = thermalParamNum;
docCells{docInd} = fdCells;
docInd += 1;

docType(docInd) = trParamNum;
docCells{docInd} = trCells;
docInd += 1;


docType(docInd) = flParamNum;
docCells{docInd} = flCells;
docInd += 1;

yaml_emitter("load/RunParam.yaml", docCells, docType); 










if(flShow == 0)
    return
end


figure
hold on

plot([-0.5,-0.5,nX-0.5,nX-0.5,-0.5],[-0.5,nY-0.5,nY-0.5,-0.5,-0.5],'r');
plot(Cx,Cy,'ro')

axis([-2,nX+2,-1,nY])
endfunction



function doc = appendCell(doc, type_, name_, val_)
  doc{length(doc)+1} = {type_, name_, val_}; 

endfunction
