[Rc, Cc, Vc] = find(Acsr);

[Rt, Ct, Vt] = find(A);
Rdiff = Rc - Rt;
Cdiff = Cc - Ct;
Vdiff = Vc-Vt;
sum(Rdiff)
sum(Cdiff)
max(Vdiff)

