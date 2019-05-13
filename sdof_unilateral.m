clearvars;
close all;
% clc;
srcpath = '~/src/matlab/nlvib/SRC';
addpath(genpath(srcpath));

%% Define system
% Parameters of the underlying linear system
mu = 1;
kappa = 1;
zeta = .03;
P = .1;

% Unilateral spring
nonlinear_elements{1} = struct(...
    'type','unilateralSpring','stiffness',10,'gap',1,'force_direction',1) ;

% Define system as single mass oscillator
oscillator = SingleMassOscillator(mu,zeta,kappa,nonlinear_elements,P);

%% Compute frequency response using harmonic balance
analysis = 'FRF';

% Analysis parameters
H = 7;         % harmonic order
N = 2^10;       % number of time samples per period
Om_s = .8;     % start frequency
Om_e = 1.6;    % end frequency

% Initial guess (from underlying linear system)
Q = (-Om_s^2*mu+1i*Om_s*zeta+kappa)\P;
y0 = [0;real(Q);-imag(Q);zeros(2*(H-1),1)];

% Solve and continue w.r.t. Om
ds = .02;
% Sopt = struct('jac','none');    % No analytical Jacobian provided here
Sopt = struct('Dscale',[1e-0*ones(size(y0));Om_s]);
% Sopt = struct();

[X_HB,Solinfo_HB] = solve_and_continue(y0,...
    @(X) HB_residual(X,oscillator,H,N,analysis),...
    Om_s,Om_e,ds,Sopt);

% Interpret solver output
Om_HB = X_HB(end,:);
Q_HB = X_HB(1:end-1,:);

%% Compute frequency response using shooting method
% Number of time samples per period
Ntd = 2^10;

% Initial guess (solution of underlying linear system). ys = [y0, yd0]
ys = [real(Q);-Om_s*imag(Q)];

% Solve and continue w.r.t. Om
ds = .02;
Sopt = struct('Dscale',[1e0*ones(size(ys));Om_s],'dynamicDscale',1);
Np = 1; % we seek period-one solutions
[X_shoot,Solinfo] = ...
    solve_and_continue(ys,...
    @(X) shooting_residual(X,oscillator,Ntd,Np,analysis),...
    Om_s,Om_e,ds,Sopt);

% Interpret solver output
Om_shoot = X_shoot(end,:);
Ys = X_shoot(1:end-1,:);
%% Determine and compare amplitude-frequency curves for both harmonic 
% balance and shooting method

% Define amplitude as magnitude of the fundamental harmonic of displacement
a_HB = sqrt(Q_HB(2,:).^2 + Q_HB(3,:).^2);

% RMS = sqrt( 1/T int_0^T f^2(t) dt )
% the integral is evaluated using Parseval’s identity, eq 2.18. We need the
% fourier coeffs. in exp. format
Q_ce = [flipud(Q_HB(2:2:end,:)+1i*Q_HB(3:2:end,:)) / 2; ...
    Q_HB(1,:); ...
    (Q_HB(2:2:end,:)-1i*Q_HB(3:2:end,:)) / 2 ];
a_HB_rms = sqrt(sum(abs(Q_ce(1:end,:)).^2,1));

% Calculate amplitude also for the results of the shooting method, and
% determine the asymptotic stability according to Floquet theory
a_shoot = zeros(size(Om_shoot));
a_shoot_rms = zeros(size(Om_shoot));
stable = zeros(size(Om_shoot));
mucrit = zeros(size(Om_shoot));
tic;
for i=1:length(Om_shoot)    
    % Evaluate solution and monodromy matrix
    [~,~,~,Y,dye_dys] = shooting_residual(X_shoot(:,i),...
        oscillator,Ntd,Np,analysis);
    
    % Determine fundamental harmonic magnitude
    Qc = fft(Y(:,1))/Ntd;
    a_shoot(i) = 2*abs(Qc(2));
%     a_shoot(i) = (max(Y(:,1)) + abs(min(Y(:,1))) )/2;
    a_shoot_rms(i) = rms(Y(:,1));
    
    % Determine stability in accordance with Floquet theory: a periodic
    % solution is stable, if all eigenvalues of the monodromy matrix remain
    % within the unit circle in the complex plane
    mucrit(i) = eigs(dye_dys,1,'lm');   % leading Floquet multiplier
    stable(i) = abs(mucrit(i))<=1; % allow for some tolerance
end
disp(['A posteriori Floquet stability analysis required ' ...
    num2str(toc) ' s.']);

% Select interesting points for further investigation:
label{1} = 'Point in regime beyond period doubling bifurcation';
[~,ind(1)] = max(a_shoot);
label{2} = 'Point in regime of stable period-one solutions';
[~,ind(2)] =  min(abs(Om_shoot-.635));
% label{3} = 'Point in regime beyond Neimark-Sacker bifurcation';
% [~,ind(3)] = min(abs(Om_shoot-.6));

%% Plot results

figure; hold on;
plot(Om_HB,a_HB,'g-');
plot(Om_shoot,a_shoot,'k.','markersize',10);
% NOTE: Some actually unstable points on the overhanging branch might be
% identified as stable, e.g. because the time discretization is still
% comparatively coarse for a reliable stability analysis.
plot(Om_shoot(~stable),a_shoot(~stable),'rx');
plot(Om_shoot(ind),a_shoot(ind),'bo');
xlabel('excitation frequency'); ylabel('response amplitude');
legend('HB','Shooting','unstable','selected','LOCATION','NW');

%% RMS amplitude
figure; hold on;
plot(Om_HB,a_HB_rms,'g-');
plot(Om_shoot,a_shoot_rms,'k.','markersize',10);
plot(Om_shoot(~stable),a_shoot_rms(~stable),'rx');
plot(Om_shoot(ind),a_shoot_rms(ind),'bo');
xlabel('excitation frequency'); ylabel('RMS response amplitude');
legend('HB','Shooting','unstable','selected','LOCATION','NW');

%% Investigate the selected response points
idx_HB = ones(length(ind),1);
[~, idx_HB(1)] = max(a_HB);
for i=1:length(ind)
    Om = Om_shoot(ind(i));
    ys = Ys(:,ind(i));
    % Simulate a few periods using direct time step integration
    tic
    n = size(oscillator.M,1);
    [t,Y] = ode45(@(t,y) [y(n+1:2*n); oscillator.M\( P*cos(Om*t) -...
        oscillator.K*y(1:n) - oscillator.D*y(n+1:2*n) - ...
        oscillator.nonlinear_elements{1}.force_direction*...
        oscillator.nonlinear_elements{1}.stiffness*(...
        oscillator.nonlinear_elements{1}.force_direction'*y(1:n)-...
        oscillator.nonlinear_elements{1}.gap)*double(...
        oscillator.nonlinear_elements{1}.force_direction'*y(1:n)-...
        oscillator.nonlinear_elements{1}.gap > 0) )],[0 2*pi/Om*1e2],...
        ys,odeset('maxstep',2*pi/Om/2e2));
    disp(['Conventional forward numerical integration for a single (!) ' ...
        'frequency from initial conditions of Shooting required ' ...
    num2str(toc) ' s.']);
    % Illustrate results in phase projection
    figure; hold on; title(label{i});
    
    % Plot computed periodic solution (stable or unstable).
    % We take here the results from the shooting method, but could of
    % course also synthesize the results from the harmonic balance method
    [~,~,~,Y_shoot] = shooting_residual([ys;Om],...
        oscillator,Ntd,Np,analysis);
    plot(Y_shoot(:,1),Om*Y_shoot(:,n+1),'g-','linewidth',2);
    
    % Select last few periods of direct time step integration
    ILP = t>(t(end)-5*2*pi/Om);
    plot(Y(ILP,1),Y(ILP,n+1),'k--');
    
    % HB
    idx = idx_HB(i);
    tau = (0:2*pi/N:2*pi-2*pi/N)';  % 2*pi/N*(0:N-1)'
    E_NH = exp(1i*tau*(-H:H));
    q = real(E_NH * Q_ce(:,idx));
    % differential operator
    nabla = diag(1i*(-H:H));
    qdot = real(E_NH *Om_HB(idx) * nabla * Q_ce(:,idx));
    plot(q,qdot,'--r')
    
    xlabel('q_1'); ylabel('u_1');
    legend('Shooting','forward num. integration', 'HB');
end
