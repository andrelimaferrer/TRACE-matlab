function base = trace_somp_prepare(training,Nr,theta_grid)
%TRACE_SOMP_PREPARE Build the joint EM/spatial dictionary for the SOMP baseline.
%
%   base = TRACE_SOMP_PREPARE(training,Nr,theta_grid)
%
% Assembles the dictionary that the greedy baseline of Section V searches.
% Classical hybrid-array channel estimation poses the problem as sparse
% recovery over an angular dictionary and solves it by simultaneous OMP
% [alkhateeb2014channel]. Transplanted to the tri-hybrid model, the
% unknown is jointly sparse in the EM state index m and the angular grid
% index g, so the dictionary is the Kronecker product
%
%     D = F_EMbar kron (F_RF A_grid)     in C^((T*NRF) x (M*G)),
%
% whose column (m-1)*G + g is kron(fbar_{EM,m}, F_RF a(theta_g)).
%
% THIS IS THE POINT OF THE COMPARISON
%   The baseline must search all M*G atoms jointly. TRACE fits the tensor
%   first, after which eqs. (23) and (24) DECOUPLE into one M-atom and
%   one G-atom problem per path. For the configuration of Fig. (b)
%   (M = 24, G = 1001) that is 24024 atoms versus 24 + 1001. The
%   dictionary built here is also what makes the baseline memory-hungry:
%   D has (T*NRF) x (M*G) entries.
%
% The ordering matches TRACE_UNFOLD(.,3).': rows are indexed with the RF
% chain fast and the EM probing state slow, consistent with eq. (7).
%
% INPUTS
%   training    struct with .FEMbar (T x M) and .FRF (NRF x Nr)
%   Nr          number of antennas
%   theta_grid  angular grid Theta, G = numel(theta_grid); the SAME grid
%               is given to TRACE_PREPARE_RECOVERY so that neither method
%               enjoys a finer grid than the other
%
% OUTPUT
%   base.D          (T*NRF) x (M*G) joint dictionary
%   base.Dnorm      the same with unit-norm columns (matched filter)
%   base.theta_grid, base.G, base.M   for decoding the selected indices
%
% NOTE ON kron
%   kron here is the MATLAB built-in. Some tensor toolboxes ship a
%   shadowing kron.m; this package depends on none of them.
%
% PAPER CROSS-REFERENCE
%   Section V, greedy dictionary baseline; contrast with eqs. (23)-(24)
%
% See also TRACE_SOMP_BASELINE, TRACE_PREPARE_RECOVERY.

Agrid=ula_steering(Nr,theta_grid);
PhiA=training.FRF*Agrid;
D=kron(training.FEMbar,PhiA); % columns kron(FEM(:,m),PhiA(:,g))
base.theta_grid=theta_grid(:).'; base.G=numel(theta_grid); base.M=size(training.FEMbar,2);
base.D=D; base.Dnorm=D./max(vecnorm(D),eps);
end
