function cache = trace_prepare_recovery(training,Nr,theta_grid)
%TRACE_PREPARE_RECOVERY Precompute the dictionaries used by Algorithm 1.
%
%   cache = TRACE_PREPARE_RECOVERY(training,Nr,theta_grid)
%
% Builds, once per experiment, the two normalized dictionaries searched by
% the recovery operators of Steps 3 and 4 of Algorithm 1:
%
%   P_EM (eq. (23)) matches each column of T_H_hat against the M columns
%        of F_EMbar                                  -> cache.FEMnorm
%   P_A  (eq. (24)) matches each column of T_A_hat against the G
%        RF-projected steering vectors F_RF a(theta) -> cache.PhiRFnorm
%
% Both dictionaries depend only on the training matrices and the angular
% grid, never on the data, so they are computed ONCE and reused across
% every Monte-Carlo realization and SNR point. Rebuilding them inside the
% loop would dominate the runtime without changing any result.
%
% WHY THE COLUMNS ARE NORMALIZED
%   Eqs. (23) and (24) maximize the magnitude of an inner product. With
%   unnormalized atoms the search would be biased toward high-energy
%   columns; normalizing turns each maximization into a correlation and
%   makes the two operators genuinely matched filters. Note that
%   F_RF a(theta) does NOT have constant norm over theta (this is the
%   deep-null effect flagged by TRACE_DESIGN_TRAINING), so the
%   normalization is essential for P_A, not merely cosmetic.
%
% SEPARABILITY
%   This is the structural advantage of fitting the tensor first: the two
%   searches decouple into one M-atom and one G-atom problem per path,
%   whereas the greedy baseline must search the joint G*M-atom dictionary
%   built by TRACE_SOMP_PREPARE.
%
% INPUTS
%   training    struct with .FEMbar (T x M) and .FRF (NRF x Nr)
%   Nr          number of antennas
%   theta_grid  angular grid Theta, in radians (G = numel(theta_grid))
%
% OUTPUT
%   cache.FEMnorm    T x M, unit-norm columns of F_EMbar
%   cache.theta_grid 1 x G grid
%   cache.Agrid      Nr x G steering dictionary
%   cache.PhiRFnorm  NRF x G, unit-norm columns of F_RF A_grid
%
% PAPER CROSS-REFERENCE
%   eqs. (23), (24); Algorithm 1, Steps 3-4
%
% See also TRACE_PHYSICAL_RECOVERY, TRACE_SOMP_PREPARE, ULA_STEERING.

F=training.FEMbar;
cache.FEMnorm=F./max(vecnorm(F),eps);
cache.theta_grid=theta_grid(:).';
cache.Agrid=ula_steering(Nr,cache.theta_grid);
Phi=training.FRF*cache.Agrid;
cache.PhiRFnorm=Phi./max(vecnorm(Phi),eps);
end
