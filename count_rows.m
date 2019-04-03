a = load("TFD_RA.txt");


setvals = zeros(1,10);
for i = 1:max(a)
  tot = sum(a == i);
  if(tot > 0)
    setvals(tot) += 1;
  end
end

