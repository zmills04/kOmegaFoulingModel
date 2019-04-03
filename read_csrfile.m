function [Inds] = read_csrfile(prefix)
  
  JA = load(sprintf("%s_JA.txt",prefix)) + 1;
  IA = load(sprintf("%s_IA.txt",prefix));
  rows = length(IA) - 1;
  Inds = zeros(rows);
  curtot = 0;
  for i = 2:length(IA)
    numel = IA(i) - IA(i-1);
    if(numel == 0)
      continue;
    end
    for j = IA(i-1)+1:IA(i)
      Inds(i-1, JA(j)) = 1;
    end
 end   
    
  
  
  
  
  
  