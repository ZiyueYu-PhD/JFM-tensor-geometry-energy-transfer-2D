%% demo_energy_flux.m
% Minimal, self-contained reproduction of the two core quantities of the paper
% from one archived velocity set:
%   (1) the filter-space total energy flux  Pi(L) = -tau_ij * S_ij
%       (space- and time-averaged) as a function of the filter scale L;
%   (2) the PDF of the alignment angle theta between the extensional
%       eigenvectors of the filtered rate-of-strain tensor and of the
%       sub-grid stress tensor, at L = 7 cm.
%
% Run from the repository root.  Only MATLAB base functions are needed.
% Note: the paper evaluates these quantities on an analysis subdomain
% (dashed box in figure 1) and on the full-length records; here the full
% 20 cm x 20 cm field and the 21 archived frames are used, so the numbers
% illustrate the method rather than reproduce the captions exactly.

clear; clc; close all
addpath('4_functions');

file  = '2_data/cellular/Cellular_Inverse_set27.mat';   % any file in 2_data/
S     = load(file);
nf    = numel(S.frames);
n     = size(S.px,1);
d0    = 20/n;          % grid spacing used by the filter (field width / points), cm
sigma = 0.5;           % Gaussian conditioning of the sharp spectral filter

%% (1) energy flux Pi(L) vs filter scale
Llist = 1:12;                              % filter scale L (cm)
Pi    = zeros(size(Llist));
for iL = 1:numel(Llist)
    acc = 0;
    for k = 1:nf
        Q   = fst_geo_efficiency_2D(S.px, S.py, S.pvx(:,:,k), S.pvy(:,:,k), Llist(iL), sigma);
        acc = acc + mean(Q, 'all', 'omitnan');
    end
    Pi(iL) = acc/nf;
end

figure('Name','Energy flux');
plot(Llist, Pi, 'o-', 'LineWidth', 1.2); hold on; yline(0, 'k:');
xlabel('L (cm)'); ylabel('\Pi(L)  (cm^2 s^{-3})');
title(strrep(file, '_', '\_'), 'FontSize', 9);
grid on

%% (2) PDF of the stress / strain alignment angle theta at L = 7 cm
L     = 7;                                 % cm
theta = [];
for k = 1:nf
    [~,~,~,~,~,ang_strain] = Cal_eigvector_eigvalue_strainandstress( ...
        S.px, S.py, S.pvx(:,:,k), S.pvy(:,:,k), 1, L, d0, 1);   % type 1: rate of strain
    [~,~,~,~,~,ang_stress] = Cal_eigvector_eigvalue_strainandstress( ...
        S.px, S.py, S.pvx(:,:,k), S.pvy(:,:,k), 0, L, d0, 1);   % type 0: sub-grid stress
    th            = abs(ang_stress - ang_strain);               % degrees
    th(th > 90)   = 180 - th(th > 90);                          % fold to [0, 90]
    theta         = [theta; th(:)/180];                         % in units of pi
end

figure('Name','PDF of theta');
histogram(theta, 0:0.025:0.5, 'Normalization', 'pdf');
xlabel('\theta / \pi'); ylabel('PDF');
title(sprintf('L = %g cm', L));

fprintf('Pi(L) [cm^2/s^3]:\n'); disp([Llist(:) Pi(:)]);
fprintf('mean theta/pi = %.3f\n', mean(theta));
