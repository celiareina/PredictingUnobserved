%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of the gradient of quartic potential 
% at x= [x1(t1),	x1(t2),     ..., x1(tNt);
%        x2(t1),    x2(t2),     ..., x2(tNt);
%        ...
%        xN(t1),    xN(t2),     ..., xN(tNt)]
%        lambda(t1),lambda(t2), ..., lambda(tNt)] size of (N+1, Nt)
% derivative wrt x1, ..., xN
% The result has size (N, N, Nt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d2Vx=D2V_x_Quartic(k2, k4, N, x)

    if ismatrix(x) == 0
        disp('x is not a matrix!');
        stop
    end
    [Nx, Nt]=size(x);
    if Nx ~= N+1
        disp('Size of x is wrong!');
        stop
    end
    
    d2Vbs=@(u) k2 + 3 * k4 * u.^2;
    
    d2Vx = zeros(N,N,Nt); 
    
    d2Vx(1,1,:) = d2Vbs( x(1,:) - 0 ) + d2Vbs( x(2,:) - x(1,:) );
    if N > 1
        for i = 2:N
            d2Vx(i,i,:) = d2Vbs( x(i,:) - x(i-1,:) ) + d2Vbs( x(i+1,:) - x(i,:) );
        end
        for i = 2:N-1
            d2Vx(i,i+1,:) = - d2Vbs( x(i+1,:) - x(i,:) );
            d2Vx(i+1,i,:) = -d2Vbs( x(i,:) - x(i-1,:) );
        end
    end
    
end