function [t1] = getDbgInfo(num, name)
name2 = sprintf("koDbg%d_%s.txt",num,name); 
t1 = load(name2);
t1max = max(max(t1));
t1min = min(min(t1));
t1absMax = max(max(abs(t1)));
fprintf("For %s: max = %g, min = %g, maxabs = %g\n", name, t1max, t1min, t1absMax);
