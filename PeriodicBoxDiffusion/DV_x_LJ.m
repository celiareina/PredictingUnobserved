%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of the gradient of total potential 
% at x= [x1(t1),	x1(t2),     ..., x1(tNt);
%        x2(t1),    x2(t2),     ..., x2(tNt);
%        ...
%        xNDof(t1),    xNDof(t2),     ..., xNDof(tNt)] size of (NDof, Nt)
% derivative wrt x1, ..., xNDof, where NDof = N * NDim
% The result has size (NDof, Nt)

% L-J interatomic potential 
% V(r) = 4 * Vm * ( (sigma / |r|)^12 - (sigma / |r|)^6 )
% with perodic boundary condition in box with size L
% Notice that the code only counts the nearest distance between two 
% particles in the periodic box.
% dV(r)/dri = - 24 * Vm * ( 2* (sigma / |r|)^12 - (sigma / |r|)^6 ) * ri / |r|^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dVx=DV_x_LJ(Vm, sigma, N, index_A, NDim, x, L)
    
    if ismatrix(x) == 0
        disp('x is not a matrix!');
        stop
    end
    [NDoF, Nt]=size(x);
    if NDoF ~= N * NDim
        disp('Size of x is wrong!');
        stop
    end
    
    sigma_AA = sigma(1);
    sigma_AB = sigma(2);
    sigma_BB = sigma(3);    
    
    Vm_AA = Vm(1);
    Vm_AB = Vm(2);
    Vm_BB = Vm(3);
    
    if max([sigma_AA, sigma_AB, sigma_BB]) > L/2
        disp('Periodic boundary condition needs to be modified!');
        stop
    end
    
    dVs = @(dr, norm_dr, Vm, sigma) - 24 * Vm * ( 2* (sigma ./ norm_dr).^12 - (sigma ./ norm_dr).^6 ) .* dr ./ norm_dr.^2;
     
    dVs_AA=@(dr, norm_dr) dVs(dr, norm_dr, Vm_AA, sigma_AA);
    dVs_AB=@(dr, norm_dr) dVs(dr, norm_dr, Vm_AB, sigma_AB);
    dVs_BB=@(dr, norm_dr) dVs(dr, norm_dr, Vm_BB, sigma_BB);
    
    dVx = zeros(N * NDim,Nt);
    for i = 1:N-1
        i_index = (i-1)*NDim+1:i*NDim;
        for j = i+1:N
            j_index = (j-1)*NDim+1:j*NDim;
            ri = x(i_index, :);
            rj = x(j_index, :);
            drij = rj - ri;
            
            % apply periodic boundary condition
%             for k = 1:NDim
%                 if drij(k) > L/2
%                     drij(k,iT) = drij(k,iT) - L;
%                 elseif drij(k) < -L/2
%                     drij(k,iT) = drij(k,iT) + L;
%                 end
%             end
            drij = mod(drij + L/2, L) - L/2;
            
            norm_drij = sqrt(sum(drij.^2, 1));
            
            % force by i-th particle on j-th particle
            if ismember(i, index_A)
                if ismember(j, index_A)
                    dV_ij = dVs_AA(drij, norm_drij + 1e-50);     % force between A-A
                else
                    dV_ij = dVs_AB(drij, norm_drij + 1e-50);     % force between A-B
                end
            else
                if ismember(j, index_A)
                    dV_ij = dVs_AB(drij, norm_drij + 1e-50);     % force between A-B
                else
                    dV_ij = dVs_BB(drij, norm_drij + 1e-50);     % force between B-B
                end
            end

            dVx(j_index, :) = dVx(j_index, :) + dV_ij;
            dVx(i_index, :) = dVx(i_index, :) - dV_ij;
        end
    end
    
end