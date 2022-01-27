%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of the gradient of quartic potential 
% at x= [x1(t1),	x1(t2),     ..., x1(tNt);
%        x2(t1),    x2(t2),     ..., x2(tNt);
%        ...
%        xN(t1),    xN(t2),     ..., xN(tNt)]
%        lambda(t1),lambda(t2), ..., lambda(tNt)] size of (N+1, Nt)
% derivative wrt x1, ..., xN
% The result has size (N, N, Nt)

% L-J interatomic potential 
% V(r) = 4 * Vm * ( (sigma / |r|)^12 - (sigma / |r|)^6 )
% with perodic boundary condition in box with size L
% dV(r)/dri = - 24 * Vm * ( 2 * (sigma / |r|)^12 - (sigma / |r|)^6 ) * ri / |r|^2
% d2V(r)/dridrj = 24 * Vm * ( 26 * (sigma / |r|)^12 - 7 * (sigma / |r|)^6 ) * ri rj / |r|^4
%               - 24 * Vm * ( 2 * (sigma / |r|)^12 - (sigma / |r|)^6 ) / |r|^2 * ( delta_ij - ri rj / |r|^2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d2Vx=D2V_x_LJ(Vm, sigma, N, index_A, NDim, x, L)

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
    
    d2Vs=@(drdr, norm_dr, Vm, sigma) 24 * Vm * ( 26 * (sigma ./ norm_dr).^12 - 7 * (sigma ./ norm_dr).^6 ) ...
        .* drdr ./ norm_dr.^4  -  24 * Vm * ( 2 * (sigma ./ norm_dr).^12 - (sigma ./ norm_dr).^6 ) ...
        ./ norm_dr.^2 .* ( eye(NDim) - drdr ./ norm_dr.^2 );

    d2Vs_AA=@(drdr, norm_dr) d2Vs(drdr, norm_dr, Vm_AA, sigma_AA);
    d2Vs_AB=@(drdr, norm_dr) d2Vs(drdr, norm_dr, Vm_AB, sigma_AB);
    d2Vs_BB=@(drdr, norm_dr) d2Vs(drdr, norm_dr, Vm_BB, sigma_BB);
    
    
    d2Vx = zeros(NDoF,NDoF,Nt); 
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
            
            drdr = zeros(NDim,NDim,Nt); 
            for p=1:NDim
                for q=1:NDim
                    drdr(p,q,:) = drij(p,:) .* drij(q,:);
                end
            end
            
            % force by i-th particle on j-th particle
            d2V_ij = zeros(NDim, NDim, Nt);
            if ismember(i, index_A)
                if ismember(j, index_A)
                    for iT=1:Nt
                        d2V_ij(:,:,iT) = d2Vs_AA(drdr(:,:,iT), norm_drij(iT) + 1e-50);     % force between A-A
                    end
                else
                    for iT=1:Nt
                        d2V_ij(:,:,iT) = d2Vs_AB(drdr(:,:,iT), norm_drij(iT) + 1e-50);     % force between A-B
                    end
                end
            else
                if ismember(j, index_A)
                    for iT=1:Nt
                        d2V_ij(:,:,iT) = d2Vs_AB(drdr(:,:,iT), norm_drij(iT) + 1e-50);     % force between A-B
                    end
                else
                    for iT=1:Nt
                        d2V_ij(:,:,iT) = d2Vs_BB(drdr(:,:,iT), norm_drij(iT) + 1e-50);     % force between B-B
                    end
                end
                
                d2Vx(i_index, i_index, :) = d2Vx(i_index, i_index, :) + d2V_ij;
                d2Vx(i_index, j_index, :) = d2Vx(i_index, j_index, :) - d2V_ij;
                d2Vx(j_index, i_index, :) = d2Vx(j_index, i_index, :) - d2V_ij;
                d2Vx(j_index, j_index, :) = d2Vx(j_index, j_index, :) + d2V_ij;
            end
        end
    end
    
end