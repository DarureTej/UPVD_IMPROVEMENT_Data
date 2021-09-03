function [C,Ceq] = constraints_SHWS(x,u_opt,distp,N)
% This function defines the constraints on the system:
    % [xmin, xmax] denote lower and upper boundaries for the manipualed
    % variables. 


x_min = [88;62.5;49.5];
x_max =  [91;64;51.5];

x0 =x; 
x_N= x; 
tspan = [0 N];
  
for i = 2: N
     [t,x_pred]  = ode45_dynamics(@dynamics_SHWS,u_opt(i-1:i),distp(:,i),tspan,x0);
     x_N = [ x_pred(1,:)'; x_N];
end

C=[repmat(x_min,N,1) - x_N;x_N-repmat(x_max,N,1) ];
Ceq = 0;
end

