function out = getVector(name, len)
  out = readbin(name);
  out = out(1:len);
  
  
  