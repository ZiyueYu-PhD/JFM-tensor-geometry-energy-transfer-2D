%% demo_energy_flux.m
% Minimal, self-contained reproduction of the core quantities of the paper
% from the archived velocity sets of one experiment day (cellular by default):
%   (1) the filter-space total energy flux  Pi(L) = -tau_ij * S_ij
%       (space- and time-averaged) as a function of the filter scale L,
%       for the four configurations: pure perturbation (rod array only),
%       inverse, forward and neutral;
%   (2) the PDF of the alignment angle theta^(L) between the extensional
%       eigenvectors of the filtered rate-of-strain tensor and of the
%       sub-grid stress tensor, at L = 7 cm, for the same four configurations;
%   (3) for inverse / forward / neutral, one frame of the fields behind
%       these numbers: the local flux Pi^(L)(x,y), the largest eigenvalue
%       and extensional eigenvector of the filtered rate of strain, and the
%       largest eigenvalue and eigenvector of the sub-grid stress.
%
% Run from the repository root.  Only MATLAB base functions are needed.
% Note: the paper evaluates these quantities on an analysis subdomain
% (dashed box in figure 1) and on the full-length records; here the full
% 20 cm x 20 cm field and the 21 archived frames are used, so the numbers
% illustrate the method rather than reproduce the captions exactly.

clear; clc; close all
addpath('4_functions');

cfg   = 'cellular';                       % 'cellular' or 'shear'
sigma = 0.5;                              % Gaussian conditioning of the sharp spectral filter
Llist = 1:12;                             % filter scales for the flux scan (cm)
Lpdf  = 7;                                % filter scale for theta and the field plots (cm)
kshow = 1;                                % frame used for the field plots (camera frame 101)

% the four configurations archived for each experiment day
tags   = {'PurePerturbation', 'Inverse', 'Forward', 'Neutral'};
labels = {'pure perturbation (rod array only)', 'inverse (\theta \approx 0)', ...
          'forward (\theta \approx \pi/2)',     'neutral (\theta \approx \pi/4)'};
setno  = struct('cellular', [21 27 32 39], 'shear', [24 30 39 48]);
Cfg    = [upper(cfg(1)) cfg(2:end)];
nset   = numel(tags);

D = cell(1, nset);
for i = 1:nset
    file = sprintf('2_data/%s/%s_%s_set%d.mat', cfg, Cfg, tags{i}, setno.(cfg)(i));
    D{i} = load(file);
    D{i}.file = file;
    D{i}.d0   = 20/size(D{i}.px, 1);     % grid spacing used by the filter (field width / points), cm
end
nf = numel(D{1}.frames);

%% (1) energy flux Pi(L) vs filter scale, all configurations
Pi = zeros(nset, numel(Llist));
for i = 1:nset
    S = D{i};
    for iL = 1:numel(Llist)
        acc = 0;
        for k = 1:nf
            Q   = fst_geo_efficiency_2D(S.px, S.py, S.pvx(:,:,k), S.pvy(:,:,k), Llist(iL), sigma);
            acc = acc + mean(Q, 'all', 'omitnan');
        end
        Pi(i, iL) = acc/nf;
    end
end

figure('Name', 'Energy flux');
plot(Llist, Pi', 'o-', 'LineWidth', 1.2); hold on; yline(0, 'k:');
xlabel('L (cm)'); ylabel('\Pi(L)  (cm^2 s^{-3})');
legend(labels, 'Location', 'best');
title(sprintf('%s background, space- and time-averaged flux', cfg));
grid on

%% (2) PDF of the stress / strain alignment angle theta at L = Lpdf, all configurations
theta = cell(1, nset);
for i = 1:nset
    S  = D{i};
    th_all = [];
    for k = 1:nf
        [~,~,~,~,~,ang_strain] = Cal_eigvector_eigvalue_strainandstress( ...
            S.px, S.py, S.pvx(:,:,k), S.pvy(:,:,k), 1, Lpdf, S.d0, 1);   % type 1: rate of strain
        [~,~,~,~,~,ang_stress] = Cal_eigvector_eigvalue_strainandstress( ...
            S.px, S.py, S.pvx(:,:,k), S.pvy(:,:,k), 0, Lpdf, S.d0, 1);   % type 0: sub-grid stress
        th          = abs(ang_stress - ang_strain);                     % degrees
        th(th > 90) = 180 - th(th > 90);                                % fold to [0, 90]
        th_all      = [th_all; th(:)/180];                              % in units of pi
    end
    theta{i} = th_all;
end

figure('Name', 'PDF of theta');
hold on
for i = 1:nset
    histogram(theta{i}, 0:0.025:0.5, 'Normalization', 'pdf', ...
              'DisplayStyle', 'stairs', 'LineWidth', 1.2);
end
xline(0.25, 'k:');                       % theta = pi/4: no net flux for given sigma, gamma
xlabel('\theta^{(L)} / \pi'); ylabel('PDF');
legend(labels, 'Location', 'best');
title(sprintf('%s background, L = %g cm', cfg, Lpdf));

%% (3) fields at one frame: local flux, strain and stress eigenstructure
% columns: Pi^(L)(x,y) | rate of strain: lambda_max + extensional eigenvector
%          | sub-grid stress: lambda_max + eigenvector
% Eigenvectors are axial (sign-free), so they are drawn as short bars.
cmap_div = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ...
            ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];    % blue - white - red
for i = 2:nset                                        % inverse, forward, neutral
    S = D{i};  k = kshow;
    Q = fst_geo_efficiency_2D(S.px, S.py, S.pvx(:,:,k), S.pvy(:,:,k), Lpdf, sigma);
    [~,~,ex_s,ey_s,lam_s] = Cal_eigvector_eigvalue_strainandstress( ...
        S.px, S.py, S.pvx(:,:,k), S.pvy(:,:,k), 1, Lpdf, S.d0, 1);
    [~,~,ex_t,ey_t,lam_t] = Cal_eigvector_eigvalue_strainandstress( ...
        S.px, S.py, S.pvx(:,:,k), S.pvy(:,:,k), 0, Lpdf, S.d0, 1);

    figure('Name', sprintf('Fields: %s', tags{i}), 'Position', [100 100 1500 480]);
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s %s, camera frame %d, L = %g cm', cfg, labels{i}, S.frames(k), Lpdf), ...
          'Interpreter', 'tex');

    ax = nexttile;                                          % local energy flux
    imagesc(S.px(1,:), S.py(:,1), Q); axis xy equal tight
    cl = robust_max(Q); clim([-cl cl]); colormap(ax, cmap_div); colorbar
    title('\Pi^{(L)}(x,y)  (cm^2 s^{-3})'); xlabel('x (cm)'); ylabel('y (cm)')

    ax = nexttile;                                          % rate of strain
    imagesc(S.px(1,:), S.py(:,1), lam_s); axis xy equal tight
    clim([0 robust_max(lam_s)]); colormap(ax, parula); colorbar; hold on
    draw_bars(S.px, S.py, ex_s, ey_s, 6, 0.8);
    title('rate of strain: \lambda_{max} (s^{-1}), extensional eigenvector'); xlabel('x (cm)')

    ax = nexttile;                                          % sub-grid stress
    imagesc(S.px(1,:), S.py(:,1), lam_t); axis xy equal tight
    clim([0 robust_max(lam_t)]); colormap(ax, parula); colorbar; hold on
    draw_bars(S.px, S.py, ex_t, ey_t, 6, 0.8);
    title('sub-grid stress: \lambda_{max} (cm^2 s^{-2}), eigenvector'); xlabel('x (cm)')
end

%% console summary
fprintf('Pi(L) [cm^2/s^3], columns = L, %s:\n', strjoin(tags, ', '));
disp([Llist(:) Pi']);
for i = 1:nset
    fprintf('%-18s  Pi(L=%g) = %+.3e   mean theta/pi = %.3f\n', ...
            tags{i}, Lpdf, Pi(i, Llist == Lpdf), mean(theta{i}));
end

%% local helpers
function draw_bars(px, py, ex, ey, step, len)
% draw the axial unit vectors (ex,ey) as centred bars of length len (cm),
% one every `step` grid points
    xs = px(1:step:end, 1:step:end);  ys = py(1:step:end, 1:step:end);
    ex = ex(1:step:end, 1:step:end);  ey = ey(1:step:end, 1:step:end);
    plot([xs(:) - 0.5*len*ex(:), xs(:) + 0.5*len*ex(:)]', ...
         [ys(:) - 0.5*len*ey(:), ys(:) + 0.5*len*ey(:)]', 'k', 'LineWidth', 0.8);
end

function m = robust_max(A)
% 99th percentile of |A| (ignoring NaN), used as a colour-scale limit
    a = sort(abs(A(~isnan(A))));
    m = a(max(1, round(0.99*numel(a))));
    if m == 0, m = 1; end
end
