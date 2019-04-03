function [] = saveParamFile(nX, nY, PipeRadius, DeviceID, nN_ls)

fname = "load/RunParam.yaml";

%  ptypes
dbl_t = 1;
int_t = 2;
str_t = 3;
bool_t = 4;

% doc types
systemParamNum = 1;
fluidParamNum = 2;
thermalParamNum = 3;
kOmegaParamNum = 4;
lsParamNum = 5;
trParamNum = 6;
flParamNum = 7;

% docCells = { document1 (n1 element cell array), document2 (n2 element cell array), ...}
% document = {param1 (3 element cell array), param2 (3 element cell array),...}
% param = {ptype (int), name (string), val (ptype)}




docType = [];
docCells = {};
docInd = 1;

%system params
doc = {};
ind = 1;
docType(docInd) = systemParamNum;


doc{ind} = {int_t, "Device ID", DeviceID}; 
ind += 1;

doc{ind} = {int_t, "nX", nX}; 
ind += 1;

doc{ind} = {int_t, "nY", nY}; 
ind += 1;

doc{ind} = {dbl_t, "Pipe radius", PipeRadius}; 
ind += 1;

docCells{docInd} = doc;
docInd += 1;

%LS params
doc = {};
ind = 1;
docType(docInd) = lsParamNum;

doc{ind} = {int_t, "nN", nN_ls}; 
ind += 1;

docCells{docInd} = doc;
docInd += 1;

yaml_emitter(fname, docCells, docType); 





