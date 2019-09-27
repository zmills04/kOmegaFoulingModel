load debugOutData.txt
pos = debugOutData(1:7:end,:);
uadv = debugOutData(2:7:end,:);
uth = debugOutData(3:7:end,:);
U00 = debugOutData(4:7:end,:);
U10 = debugOutData(5:7:end,:);
U01 = debugOutData(6:7:end,:);
U11 = debugOutData(7:7:end,:);

plot(pos(:,1),pos(:,2),'r.')
fprintf("Max Adv Vel = (%g,%g), Max Thermo Vel = (%g, %g)\n",max(abs(uadv(:,1))), max(abs(uadv(:,2))),max(abs(uth(:,1))), max(abs(uth(:,2))));