function [out] = readbin(fname, y)
  if(nargin == 1)
    y = 1;
  end
  
  fname = ['load/',fname,'.bin'];

  bfile = fopen(fname);
  bdata = fread(bfile, 'double');
  x = length(bdata)/y;
  if((round(x) - x) < 1e-10)
    out = reshape(bdata,x,y); 
  else
    out = bdata;
  end
  fclose(bfile);