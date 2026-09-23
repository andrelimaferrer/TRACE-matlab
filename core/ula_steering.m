function A = ula_steering(Nr,theta)
%ULA_STEERING Half-wavelength ULA steering vectors a(theta).
%
%   A = ULA_STEERING(Nr,theta)
%
% Builds the steering matrix A(theta) = [a(theta_1),...,a(theta_L)] that
% appears in eq. (2) of the paper, with
%
%     [a(theta)]_n = exp(-j*pi*n*sin(theta)) / sqrt(Nr),   n = 0,...,Nr-1,
%
% i.e. a uniform linear array with half-wavelength element spacing. The
% 1/sqrt(Nr) factor makes every steering vector unit norm, so that the
% path gains alpha carry all the channel energy; this normalization is
% what makes the SNR definition of Section V and the NMSE comparable
% across array sizes.
%
% INPUTS
%   Nr     number of receive antenna elements
%   theta  angles of arrival in radians; scalar, row or column vector
%
% OUTPUT
%   A      Nr x numel(theta) steering matrix, unit-norm columns
%
% PAPER CROSS-REFERENCE
%   a(theta)      : eq. (2)
%   A(theta)      : eq. (2), and factor T_A = F_RF A(theta) in eq. (12)
%   grid version  : the dictionary {a(theta) : theta in Theta} searched by
%                   the beamspace operator P_A in eq. (24)
%
% See also TRACE_PREPARE_RECOVERY, TRACE_GENERATE_CLEAN_DATA.

n=(0:Nr-1).';
theta=theta(:).';
A=exp(-1j*pi*n*sin(theta))/sqrt(Nr);
end
