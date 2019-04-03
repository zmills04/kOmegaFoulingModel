function convert_to_rcv(Fmat, name)
  
  [Row,Col,Val] = find(Fmat);
  
  Out = [Row,Col,Val];
  dlmwrite(sprintf("%s.txt",name),Out);
  
  
  
  
  
  
  
  
  
  
  
  
  