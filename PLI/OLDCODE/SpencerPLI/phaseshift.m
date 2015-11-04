function x_shifted = phaseshift(x,phi,delta_imag)
% PHASESHIFT Frequency-domain phase shift (time-domain delay) of a signal
%
%   X_SHIFTED = PHASESHIFT(X,PHI)
%   Shift the signal in X by PHI radians across all frequencies, returning
%   the result in X_SHIFTED.  If X is a matrix, PLI operates along columns.
%   PHI should have one entry per column of X.  If X a single column vector
%   and PHI has multiple entries, X_SHIFTED will contain one column per
%   element of PHI.  Likewise for the case where X is a matrix but PHI is
%   scalar.
%
%   PHASESHIFT(X,PHI,DELTA_IMAG)
%   Specify the maximum value allowable for the imaginary part of the
%   inverse-FFT of the phase-shifted signal (default 1e-10).

% check for corner case: shifting by 0, 2*pi, 4*pi, ...
if isscalar(phi)&&mod(phi,2*pi)==0
    x_shifted = x;
    return;
end

% phi is either scalar or one element per column
if isscalar(phi) && size(x,2)>1
    phi = repmat(phi,1,size(x,2));
elseif size(x,2)==1 && length(phi)>1
    x = repmat(x,1,length(phi));
end
assert(length(phi)==size(x,2),'Phi must be scalar, or have one entry per column of X');

% basic information 
nfft = size(x,1); % size of the FFT

% default value of imaginary threshold
if nargin<3||isempty(delta_imag),delta_imag=nfft*1e-2;end

% frequencies (Hz)
f = -1/2 + (0:(nfft-1))'/nfft;

% calculate the FFT of the reference signals
X_ref = fft(x,nfft);

% shift in frequency domain
sgnf = sign(f(:));
p = exp(1i*sgnf*phi);
X_shifted = X_ref.*ifftshift(p,1); % ifftshift has different default behavior for matrix input, must specify dimension

% convert back to the time domain
x_shifted = ifft(X_shifted);

% make sure the imaginary part is not larger than threshold
if max(imag(x_shifted(:)))>delta_imag
    warning('Maximum imaginary part of shifted signal (%g) larger than specified threshold (%g)',max(imag(x_shifted(:))),delta_imag);
end

% return only the real part
x_shifted = real(x_shifted);