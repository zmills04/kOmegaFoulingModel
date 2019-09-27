load lba_0.txt;
load lba_1.txt;
load lba_2.txt;
load lba_3.txt; 
load lba_4.txt;
load lba_5.txt;
load lba_6.txt;
load lba_7.txt;
load lba_8.txt;
[x,y] = size(lba_0);
lba = zeros(x,y,9);
lba(:,:,1) = lba_0;
lba(:,:,2) = lba_1;
lba(:,:,3) = lba_2;
lba(:,:,4) = lba_3;
lba(:,:,5) = lba_4;
lba(:,:,6) = lba_5;
lba(:,:,7) = lba_6;
lba(:,:,8) = lba_7;
lba(:,:,9) = lba_8;

load lbf_0.txt;
load lbf_1.txt;
load lbf_2.txt;
load lbf_3.txt; 
load lbf_4.txt;
load lbf_5.txt;
load lbf_6.txt;
load lbf_7.txt;
load lbf_8.txt;
[x,y] = size(lbf_0);
lbf = zeros(x,y,9);
lbf(:,:,1) = lbf_0;
lbf(:,:,2) = lbf_1;
lbf(:,:,3) = lbf_2;
lbf(:,:,4) = lbf_3;
lbf(:,:,5) = lbf_4;
lbf(:,:,6) = lbf_5;
lbf(:,:,7) = lbf_6;
lbf(:,:,8) = lbf_7;
lbf(:,:,9) = lbf_8;

test = lbf == lba;
sum(sum(sum(~test)))



