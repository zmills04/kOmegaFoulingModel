post = load('Par_Pos.txt');loct = load('Par_Loc.txt');timert = load('Par_Timer.txt');depTimert = load('Par_DepTimer.txt');numRept = load('Par_NumRep.txt');typet = load('Par_Type.txt');posx = [posx,post(:,1)];posy = [posy,post(:,2)];loc = [loc,loct];timer = [timer,timert];depTimer = [depTimer,depTimert];numRep = [numRep,numRept];type = [type,typet];


posx = [];posy = [];loc = [];timer = [];depTimer = [];numRep = [];type = [];



Partemp = struct('pos',pos,'loc',loc,'timer',timer,'depTimer',depTimer,...
    'numRep',numRep,'type',type);

loc = floor(Par2.pos(:,1)) + floor(Par2.pos(:,2))*1856 - 88;

[~,lt] = sort(ParT.loc);

ParS = ParT;
ParS.pos = ParT.pos(lt,:);
ParS.loc = ParT.loc(lt);
ParS.timer = ParT.timer(lt);
ParS.depTimer = ParT.depTimer(lt);
ParS.numRep = ParT.numRep(lt);
ParS.type = ParT.type(lt);








ParComp = struct('posx',[ParS.pos(:,1),Par2.pos(:,1)],...
    'posy',[ParS.pos(:,2),Par2.pos(:,2)],...
    'loc',[ParS.loc,Par2.loc],...
    'timer',[ParS.timer,Par2.timer],...
    'depTimer',[ParS.depTimer,Par2.depTimer],...
    'numRep',[ParS.numRep,Par2.numRep],...
    'type',[ParS.type,Par2.type]);

