function [fob] = objective_SHWS(x,u_opt,distp,N)
% [fob] = objective_SHWS(xref,wx,wu)% 
% Input dimentions are given as below:
% xref denoted the references and are explained as xref  = [Tcout-setpoint, Ttout_setoint]
% wx are the weights on the states and are given as xref = [Tcout-weight, Ttout_weight]
% wu are weights on inputs and given as wu = [Tcout-weight, Ttout_weight]

fob = 0; 
ts =1; 
x_N =[];
 % reference setpoints for x(1), x(3)
 xref = repmat([65 51.5],N,1);
 wx = [1e0,1e0];
 wu = [1e1;1e1]*0;
 
 % calculation of objective function
    for i = 2: N
         d= distp(:,i);
         [t,x_pred]  = ode45_dynamics(@dynamics_SHWS,u_opt(i-1:i),d,[0 ts],x);
          x_N = [x_N x_pred(end,:)'];
%           fob = fob + wx(1)*norm((xref(1)-x_pred(end,1)), 2)  + wx(2)*norm((xref(2)-x_pred(end,3)), 2)+   wu(1)*norm(u_opt(i-1), 2)+   wu(2)* norm(u_opt(i), 2);
    end

     fob = wx(1)*norm((xref(:,1) - x_N(1,:)), 2) + wx(2)*norm((xref(:,2) - x_N(3,:)), 2);
end

