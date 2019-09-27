shortenArray = 1;

xstart = 300;
xstop = 400;
ystart = 1;
ystop = 100000;


load nType.txt;
[x,y] = size(nType);
figure
hold on
load lsc.txt;
load dXarr_0.txt;
load dXarr_1.txt;
load dXarr_2.txt;
load dXarr_3.txt;
dX = zeros(x,y,8);
dX(:,:,1) = dXarr_0;
dX(:,:,2) = dXarr_1;
dX(:,:,3) = dXarr_2;
dX(:,:,4) = dXarr_3;
dX(:,:,5) = (dX(:,:,1).*dX(:,:,3));
dX(:,:,6) = (dX(:,:,2).*dX(:,:,4));
dX(:,:,7) = (dX(:,:,1).*dX(:,:,4));
dX(:,:,8) = (dX(:,:,2).*dX(:,:,3));



Solid = zeros(x,y);
Solid0 = Solid;
Fluid = Solid;
Fluid0 = Fluid;
solidBound = Fluid;
solidBound0 = Fluid;
neswBound0 = Fluid;
bound0 = Fluid;
eBound = Fluid;
wBound = Fluid;
nBound = Fluid;
sBound = Fluid;
neBound = Fluid;
nwBound = Fluid;
seBound = Fluid;
swBound = Fluid;
for i = 1:x
    for j = 1:y
        val = nType(i,j);
        
        if(mod(val,2) == 1)
            Solid(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            Fluid(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            Solid0(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            Fluid0(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            solidBound(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            solidBound0(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            neswBound0(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            bound0(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            eBound(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            wBound(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            nBound(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            sBound(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            neBound(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            swBound(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            seBound(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
        if(mod(val,2) == 1)
            nwBound(i,j) = 1;
            val = (val-1);
        end
        val = val/10;
        
    end
end

boundaries = zeros(x,y,8);
boundaries(:,:,1) = eBound;
boundaries(:,:,2) = wBound;
boundaries(:,:,3) = nBound;
boundaries(:,:,4) = sBound;
boundaries(:,:,5) = neBound;
boundaries(:,:,6) = swBound;
boundaries(:,:,7) = seBound;
boundaries(:,:,8) = nwBound;
cx = [1,-1,0,0,1,-1,1,-1];
cy = [0,0,1,-1,1,-1,-1,1];



if(shortenArray)
    if(ystart >= y || xstart >= x)
        error("cutoff start values too large\n");
    end
    ystop = min(ystop,y);
    xstop = min(xstop,x);
    
    lsc((lsc(:,1) < (xstart-3)) | (lsc(:,1) > (xstop+3)),:) = [];
    lsc((lsc(:,2) < (ystart-3)) | (lsc(:,2) > (ystop+3)),:) = [];
    plot(lsc(:,1),lsc(:,2),'r')
else
plot(lsc(end/2+1:end,1),lsc(end/2+1:end,2),'r');
plot(lsc(end/2+1:end,1),lsc(end/2+1:end,2),'r');
end

solidNodes = zeros(sum(sum(Solid)),2);
fluidNodes = zeros(sum(sum(Fluid)),2);
sboundNodes = zeros(sum(sum(solidBound)),2);
solidInd = 1;
fluidInd = 1;
sboundInd = 1;
for i = xstart:xstop
    for j = ystart:ystop
        if(solidBound(i,j) == 1)
            sboundNodes(sboundInd,:) = [i-1,j-1];
            sboundInd = sboundInd+1;
        end
        if(Solid(i,j) == 1)
            solidNodes(solidInd,:) = [i-1,j-1];
            solidInd = solidInd+1;
        end
        if(Fluid(i,j) == 1)
            fluidNodes(fluidInd,:) = [i-1,j-1];
            fluidInd = fluidInd+1;
        end
    end
end

plot(solidNodes(1:solidInd-1,1),solidNodes(1:solidInd-1,2),'ro')
plot(fluidNodes(1:fluidInd-1,1),fluidNodes(1:fluidInd-1,2),'r*')
plot(sboundNodes(1:sboundInd-1,1),sboundNodes(1:sboundInd-1,2),'k.')


for i = xstart:xstop
    for j = ystart:ystop
        for k = 1:8
            if(boundaries(i,j,k) == 1)
%                 if(k < 5)
                    plot([i-1,i-1+cx(k)*dX(i,j,k)],[j-1,j-1+cy(k)*dX(i,j,k)],'b')
%                 else
%                     plot([i-1,i-1+cx(k)],[j-1,j-1+cy(k)],'b')
%                 end
            end
        end
    end
end

for i = 2:x-1
    for j = 2:y-1
        if(Solid(i,j) == 1)
            continue;
        end
        for k = 1:8
            if(Solid(i+cx(k),j+cy(k)) == 1)
                if(boundaries(i,j,k) ~= 1)
                    sprintf("error at (%d, %d), dir %d. Should have bit flag set, but isnt\n",i,j,k);
                end
            end
            if(boundaries(i,j,k) == 1)
                if(Solid(i+cx(k),j+cy(k)) ~= 1)
                    sprintf("error at (%d, %d), dir %d. Have bitflag set when shouldnt be\n",i,j,k);
                end
            end
        end
    end
end