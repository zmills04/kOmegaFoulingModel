function [Uxx] = plot_uvels(Ux, Re, Visc, Uxx_prev)
[nx, ny] = size(Ux);
H = ny/2;
close all
Uxx = mean(Ux, 1);
hold on


Yvals = (0.5:ny-0.5) - H;
F = Re*Re*Visc*Visc/H/H/H;
Umax = H*H*F/2/Visc;



plot(Yvals, Uxx,'k')
UUU = Umax * (1. - Yvals.*Yvals/H/H);
plot(Yvals, UUU,'r')

if(nargin == 4)
  plot(Yvals,Uxx_prev,'r');
  Ux_bot_prev = Uxx_prev(1:end/2);
  Ux_top_prev = Uxx_prev(end:-1:end/2+1);
end

Ux_bot = Uxx(1:end/2);
Ux_top = Uxx(end:-1:end/2+1);
Ux_bot_theo = UUU(1:end/2);
Ux_top_theo = UUU(end:-1:end/2+1);


kappa = 0.41;
Cp = 5.5;
U0 = 0.05;

g = Re*Re*Visc*Visc/H/H/H;
u_tau = Re*Visc/H;

U_plus_bot = Ux_bot(1:end)/u_tau;
U_plus_top = Ux_top(1:end)/u_tau;
U_plus_bot_theo = Ux_bot_theo(1:end)/u_tau;
U_plus_top_theo = Ux_top_theo(1:end)/u_tau;

y = 0.5:H-0.5;
y2 = logspace(-1,2,50);
y3 = linspace(0.5,H-0.5,100);

y_plus = y*u_tau/Visc;
y2_plus = y2*u_tau/Visc;
y3_plus = y3*u_tau/Visc;

U_plus = log(y2_plus)/kappa + Cp;
% U_plus_lin = log(y3_plus) / kappa + Cp;

a = find(y2_plus < 11.44532166);
% a_lin = find(y3_plus < 11.44532166);

U_plus(a) = y2_plus(a);
% U_plus_lin(a_lin) = y3_plus(a_lin);
% U_plus_lin *= u_tau;
% U_plus_lin = [U_plus_lin,U_plus_lin(end:-1:1)];
% plot([y3-H,y3], U_plus_lin,'c')
hold off

% printf("Ulbm = %g, Utheo = %g", mean(Uxx), mean(U_plus_lin));
% if(nargin == 4)
% 	printf(", Ulbm_old = %g\n", mean(Uxx_prev));
% else
% 	printf("\n");
% end



figure
hold on
semilogx(y_plus, U_plus_top,'k')
semilogx(y_plus, U_plus_bot,'r')
semilogx(y_plus, U_plus_top_theo,'k--')
semilogx(y_plus, U_plus_bot_theo,'r--')
semilogx(y2_plus, U_plus, 'b')
if(nargin == 4)
  semilogx(y_plus, Ux_top_prev(1:end)/u_tau,'k--')
  semilogx(y_plus, Ux_bot_prev(1:end)/u_tau,'r--')
end
hold off
  
  
  
  
  
  
  