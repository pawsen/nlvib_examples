%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of a single-degree-of-freedom (SDOF)
% oscillator with cubic spring, i.e. the Duffing oscillator, under harmonic
% external forcing and light linear viscous damping. The dynamics are
% governed by the second-order ordinary diferential equation,
% 
%   mu \ddot q + zeta \dot q + kappa q + \gamma q^3 = P cos( Om * t ).

clearvars;
close all;
clc;
srcpath = '~/src/matlab/nlvib/SRC';
addpath(srcpath);
%% Parameters of the Duffing oscillator
mu = 1;
zeta = 0.05;
kappa = 1;
gamma = 0.1;
P = 0.18;
%% Compute frequency response using numerical implementation of the 
% harmonic balance method

% Analysis parameters
H = 7;      % harmonic order
N = 2^6;    % number of time samples per period
Om_s = .5;  % start frequency
Om_e = 1.6; % end frequency

% Initial guess (from underlying linear system)
Q = (-Om_s^2*mu+1i*Om_s*zeta+kappa)\P;
x0 = [0;real(Q);-imag(Q);zeros(2*(H-1),1)];

% Solve and continue w.r.t. Om
ds = .01;                       % Path continuation step size
Sopt = struct('jac','none');    % No analytical Jacobian provided here
X = solve_and_continue(x0,...
    @(X) HB_residual_Duffing(X,mu,zeta,kappa,gamma,P,H,N),...
    Om_s,Om_e,ds,Sopt);

% Determine excitation frequency and amplitude (magnitude of fundamental
% harmonic)
Om = X(end,:);
a = sqrt(X(2,:).^2 + X(3,:).^2);
%% Analytically calculate frequency response using single-term harmonic 
% balance.

% In this case, the excitation frequencies are determined depending on the
% amplitude. Hence, we specify the amplitude range and samples for which we
% determine the frequencies.
a_min = .1;
a_max = 3;
a_ana = linspace(a_min,a_max,50);

% For each amplitude, determine the two associated excitation frequencies
Om_ana = zeros(length(a_ana),2);
for i=1:length(a_ana)
    Om_ana(i,1) = sqrt(1-zeta^2/2+3*gamma*a_ana(i)^2/4 + ...
        sqrt(P^2/a_ana(i)^2+zeta^4/4-zeta^2-3*zeta^2*gamma*a_ana(i)^2/4));
    Om_ana(i,2) = sqrt(1-zeta^2/2+3*gamma*a_ana(i)^2/4 - ...
        sqrt(P^2/a_ana(i)^2+zeta^4/4-zeta^2-3*zeta^2*gamma*a_ana(i)^2/4));
end
% Only the real-valued solutions exist. Let us store the information which
% ones are valid.
valid_ana = imag(Om_ana(:,1))==0 & imag(Om_ana(:,2))==0;
%% Illustrate results
figure; hold on;
plot(Om,a,'k-');
plot(Om_ana(valid_ana,:),a_ana(valid_ana),'gx');
set(gca,'xlim',[Om_s Om_e]);
xlabel('excitation frequency \Omega'); ylabel('response amplitude a');
legend(['numerical HB, H=' num2str(H)],'analytical HB, H=1',...
    'LOCATION','NW');

%% Plot of periodic solution in time domain

% max amplitude and index
[am, idx] = max(a);
tau = (0:2*pi/N:2*pi-2*pi/N)';  % 2*pi/N*(0:N-1)'
E_NH = exp(1i*tau*(-H:H));

Q_ce = [flipud(X(2:2:end-1,idx)+1i*X(3:2:end-1,idx)) / 2; ...
    X(1,idx); ...
    (X(2:2:end-1,idx)-1i*X(3:2:end-1,idx)) / 2 ];
q = real(E_NH * Q_ce);

% differential operator
nabla = diag(1i*(-H:H));
dq = real(E_NH *X(end,idx) * nabla * Q_ce);
ddq = real(E_NH *X(end,idx)^2 * nabla^2 * Q_ce);

% Analytically for one-term HB expansion
% use eq. 1.10 for transform from polar coords to fourier coeff and eq 1.21
[~,idx_ana] = max(a_ana(valid_ana));
Om_ana_max = max(Om_ana(idx_ana,:));
theta = asin(zeta*Om_ana_max*a_ana(idx_ana)/P);
qc = a_ana(idx_ana)*cos(theta);
qs = a_ana(idx_ana)*sin(theta);
% find the phase for q_ana = 0
% -qc/qs=sin(x)/cos(x)=tan(x)
phi = atan(-qc/qs);

% eq 1.3-1.5
% T = 2*pi/Om_ana_max;
% tau_ana = Om_ana_max * (0:T/N:T-T/N);  % equal to tau
q_ana = qc*cos(tau+phi) + qs*sin(tau+phi);
dq_ana = Om_ana_max*(-qc*sin(tau+phi) + qs*cos(tau+phi));
ddq_ana = Om_ana_max^2*(-qc*cos(tau+phi) - qs*sin(tau+phi));

% Mismatch between HB=7 and HB_ana due to phase. 
figure;
hold on
plot(tau, q,'b')
plot(tau, dq,'b')
plot(tau, ddq,'b')
plot(tau, q_ana,'r--')
plot(tau, dq_ana,'r--')
plot(tau, ddq_ana,'r--')
