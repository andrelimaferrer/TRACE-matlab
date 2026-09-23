function data = trace_add_noise(clean,SNRdB)
%TRACE_ADD_NOISE Add complex AWGN at a prescribed SNR.
%
%   data = TRACE_ADD_NOISE(clean,SNRdB)
%
% Adds the noise tensor N of eq. (11) to a noiseless realization produced
% by TRACE_GENERATE_CLEAN_DATA.
%
% SNR DEFINITION (Section V)
%   "The SNR is the ratio of the average power of the noiseless pilot
%   tensor to the noise power per entry", i.e.
%
%       sigma^2 = E{|Y_clean|^2} / 10^(SNR/10),
%
%   with the average taken over ALL NRF*T*P entries of the realization.
%   The noise is circularly symmetric complex Gaussian, so each of the
%   real and imaginary parts gets variance sigma^2/2 -- hence the
%   sqrt(sigma2/2) factor below.
%
%   Because sigma^2 is recomputed from the realization at hand, the SNR is
%   exact per realization rather than only on average. This is why the
%   Monte-Carlo curves of Section V can be read as conditional on SNR.
%
% REUSING ONE CHANNEL ACROSS SNR POINTS
%   The reproduction scripts call this function repeatedly on the SAME
%   `clean` structure to sweep the SNR. That is deliberate: it makes the
%   SNR curves paired across realizations, removing channel-draw variance
%   from the comparison between operating points.
%
% INPUTS
%   clean   output of TRACE_GENERATE_CLEAN_DATA
%   SNRdB   scalar SNR in dB
%
% OUTPUT
%   data    copy of `clean` with the extra fields
%             data.Y       noisy tensor, eq. (11)
%             data.SNRdB   the SNR used
%             data.sigma2  the resulting noise variance per entry
%
% PAPER CROSS-REFERENCE
%   eq. (11) (the term calligraphic N); Section V, SNR definition
%
% See also TRACE_GENERATE_CLEAN_DATA.

data=clean;
Ps=mean(abs(clean.Yclean(:)).^2);
sigma2=Ps/10^(SNRdB/10);
noise=sqrt(sigma2/2)*(randn(size(clean.Yclean))+1j*randn(size(clean.Yclean)));
data.Y=clean.Yclean+noise; data.SNRdB=SNRdB; data.sigma2=sigma2;
end
