%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return to the reference trajectory that is linear in space with given 
% protocol (lambda) at different times for quadratic potential expansion
% xr with size (N, Nt) where Nt = length(lambda)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xr = xr_t_linear(N, lambda)
    
    Nt = length(lambda);
    
    xr = zeros(N,Nt); 
    for i = 1:N
        xr(i,:) = i/(N+1) * lambda;
    end
    
end