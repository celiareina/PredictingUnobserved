%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of the gradient of quartic potential 
% at x= [x1(t1),	x1(t2),     ..., x1(tNt);
%        x2(t1),    x2(t2),     ..., x2(tNt);
%        ...
%        xN(t1),    xN(t2),     ..., xN(tNt)]
%        lambda(t1),lambda(t2), ..., lambda(tNt)] size of (N+1, Nt)
% derivative wrt x1, ..., xN
% The result has size (N, Nt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dVx=DV_x_Quartic(k2, k4, N, x)
    
    if ismatrix(x) == 0
        disp('x is not a matrix!');
        stop
    end
    [Nx, Nt]=size(x);
    if Nx ~= N+1
        disp('Size of x is wrong!');
        stop
    end
    
    dVbs=@(u) k2 * u + k4 * u.^3;
    
    dVx = zeros(N,Nt); 
    dVx(1,:) = dVbs( x(1,:) - 0 )- dVbs( x(2,:) - x(1,:) );
    if N > 1
        dVx(2:N,:)= dVbs( x(2:N,:) - x(1:N-1,:) ) - dVbs( x(3:N+1,:) - x(2:N,:) );
    end
    
end