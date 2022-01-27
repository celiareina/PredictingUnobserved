%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of the derivative of total potential 
% V(x1, ..., xN, lambda) = V_s(x1-0) + V_s(x2-x1) + ... + V_s(xN-x(N-1)) + V_s(lambda-xN)
% with quartic interatomic potential
% V_s(u) = 1/2 * k_2 * u^2 + 1/4 * k_4 * u^4
%
% Derivative wrt x1, ..., xN and lambda
% x=[x1, ... xN, lambda];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function force=DQuarticPotential(k2,k4,N,x)

    DV=@(u) k2*u+k4*u.^3;
    force = zeros(N+1,1); 
    force(1:N)= DV([x(1); x(2:N)-x(1:N-1)])- DV(x(2:N+1)-x(1:N));
    force(N+1)= DV(x(N+1)-x(N));    
    
end