function x = get_x(g,t,t_min,t_max)
% here use the values in g to solve the IVP and return the solution at
% times t
% the ouput should be size(x) = [length(t),2]

Q = g(1);
P = g(2);
m = g(3);
a = g(4);

% sol2 = ode15s(@(t,p) ...
%               (Q+P-p)*m*p/(Q+P-p+m/a), ...
%               [t_min,t_max],P);
% p = deval(sol2,t);
% q = Q+P-p;
% x = [q',p'];

QP = Q+P;
M = m/a/QP;

sol2 = ode15s(@(t,x) m*x(1)*x(2)/(x(2)+M)*[-1;+1],...
              [t_min,max(t)],...
              [Q,P]/QP,...
              odeset('NonNegative',1:2));
x = QP * deval(sol2,t)';

end
