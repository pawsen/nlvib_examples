function R = HB_residual_Duffing(X,mu,zeta,kappa,gamma,P,H,N)
% residual of the harmonically forced duffing oscillator
%
%  μq̈ + ζq̇ + κq + γq^3 = P/2*exp(iΩt) + P/2*exp(−iΩt) = P*cos(Ωt)
%
% INPUTS:
% X = [q̂ .'_{sc_H},  Ω].':  sine-cosine fourier coeff. and ext. force freq.
% SYSTEM parameters
% μ, ζ, κ, γ, P
% HB parameters
% H: harmonic truncation order
% N: number of AFT sampling points per period

% NLvib uses the sine–cosine representation externally (for x and R).
% Internally, the †-variant of the complex-exponential representation is
% used (as introduced in Chap. 2, c.f. Eq. 2.73). The conversion between 
% these are given by eq 2.75.

% Conversion of sine-cosine to complex-exponential representation
% Q_ce = [q(-H)...q(-1), q(0), q(1)...q(H)]. eq. 2.62
% We know from eq 2.3 that the negative spectrum is the complex conjugate
% of the positive: q(−k) = conj(q(k)). eq 2.7 gives the conversion for k=>1
% and q(0) is identical for both representations.
Q_ce = [flipud(X(2:2:end-1)+1i*X(3:2:end-1)) / 2; ...
    X(1); ...
    (X(2:2:end-1)-1i*X(3:2:end-1)) / 2 ];

% Excitation frequency
Om = X(end);
% P is the magnitude of the cosine forcing
Fex_ce = [ zeros(H-1,1) ; P/2 ; 0 ; P/2; zeros(H-1,1) ];

% Specify time samples along period
tau = (0:2*pi/N:2*pi-2*pi/N)';  % 2*pi/N*(0:N-1)'
% Inverse discrete Fourier transform matrix, eq 2.64. N samples, H harmonics
E_NH = exp(1i*tau*(-H:H));

% Apply inverse discrete Fourier transform. We are now in gen. coor in time
% domain. 
q = real(E_NH * Q_ce);
% Evaluate nonlinear force in the time domain
fnl = gamma * q.^3;

% Apply discrete Fourier transform, eq 2.65
Fnl_ce = E_NH'/N * fnl;

% Dynamic force equilibrium, eq. 3.7
R_ce = ( -((-H:H)'*Om).^2*mu + 1i*(-H:H)'*Om*zeta + kappa) ...
    .*Q_ce + Fnl_ce - Fex_ce;

% Conversion of complex-exponential to sine-cosine residual
% conversion rules eqs. 2.5-2.7. Since r(-k)=conj(r(k)) it is sufficient to
% use only R_ce(0) through R_ce(H). 
R = [ real(R_ce(H+1)); real(R_ce(H+2:end)); -imag(R_ce(H+2:end))];
