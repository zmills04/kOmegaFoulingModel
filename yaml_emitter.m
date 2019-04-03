function [] = yaml_emitter(fname, docCells, docTypes) 
  
%  ptypes
  dbl_t = 1;
  int_t = 2;
  str_t = 3;
  bool_t = 4;
  dblvec_t = 5;
  intvec_t = 6;
  % doc types
  systemParamNum = 1;
  fluidParamNum = 2;
  thermalParamNum = 3;
  kOmegaParamNum = 4;
  lsParamNum = 5;
  trParamNum = 6;
  flParamNum = 7;
  
  docNames = {'System', 'Fluid', 'Thermal',...
   'kOmega', 'LS', 'Particle', 'Fouling layer'};
  
  fid = fopen(fname, 'w');
  
  fprintf(fid, "# Parameters for Fouling simulation\n");
  
  for i = 1:length(docTypes)
    docVals = docCells{i};
        
    fprintf(fid, "--- # %s Parameters\n", docNames{docTypes(i)});
    fprintf(fid, "ParamsName: %d\n", docTypes(i))
    for j = 1:length(docVals)
      param = docVals{j};
      ptype = param{1};
      pname = param{2};
      pval = param{3};
      fprintf(fid, "%s: ", pname);
      switch ptype
        case dbl_t
          fprintf(fid, "!!float %.8e\n", pval);
        case int_t
          fprintf(fid, "!!int %s\n", int2str(pval));
        case bool_t
          if(pval)
            fprintf(fid, "True\n");
          else
            fprintf(fid, "False\n");
          endif
        case str_t
          fprintf(fid, "!!str %s\n", pval);
        case dblvec_t
          fprintf(fid, "[%.8e", pval(1));
          for i = 2:length(pval)
            fprintf(fid, ", %.8e", pval(i));
          end
          fprintf(fid, "]\n");
        case intvec_t
          fprintf(fid, "[%s", int2str(pval(1)));
          for i = 2:length(pval)
            fprintf(fid, ", %s", int2str(pval(i)));
          end
          fprintf(fid, "]\n");
        otherwise
          exit("Incorrect data type identifier supplied\n");  
      endswitch
    endfor
    fprintf(fid, "\n");
  endfor
  
  fclose(fid);
endfunction      
 
         
      
      
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  