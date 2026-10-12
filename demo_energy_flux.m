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
%
% Analysis subdomain.  As in the paper, the filter is applied to the full
% 20 cm x 20 cm field, but all averages and statistics are taken over the
% analysis subdomain (dashed box in the figures): the quadrilateral that stays
% inside the zone swept by the reciprocating rod array.  Its corner
% coordinates (cm, image axes) are those of the original analysis script and
% are listed below for each set.  The field plots show the full field with
% the subdomain outlined.  The 21 archived frames replace the full-length
% records, so the numbers illustrate the method rather than reproduce the
% figure captions exactly.

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

% analysis subdomain of each set: 4 corners [x1 y1 x2 y2 x3 y3 x4 y4] in cm
% (same order as the rows of `tags`; values from the original analysis script)
corners = struct( ...
    'cellular', [ 6.45  7.82  12.85  7.82  13.00 11.63   6.61 11.32;    % set21
                  8.09  7.52  11.39  7.52  11.68 12.25   8.09 12.11;    % set27
                  7.07  7.86  12.48  7.86  12.33 11.10   7.07 11.10;    % set32
                  7.60 10.64  10.93  7.16  13.40  9.77   9.92 12.96], ... % set39
    'shear',    [ 8.72  7.76  14.76  7.76  14.76 11.70   8.72 11.76;    % set24
                  4.09 10.38   7.99  6.60  10.51  9.12   6.48 13.03;    % set30
                  4.09 10.38   7.99  6.60  10.51  9.12   6.48 13.03;    % set39
                  9.46  7.76  15.50  7.76  15.50 11.70   9.46 11.76]);  % set48

Cfg  = [upper(cfg(1)) cfg(2:end)];
nset = numel(tags);

D = cell(1, nset);
for i = 1:nset
    file = sprintf('2_data/%s/%s_%s_set%d.mat', cfg, Cfg, tags{i}, setno.(cfg)(i));
    D{i} = load(file);
    D{i}.file   = file;
    D{i}.d0     = 20/size(D{i}.px, 1);   % grid spacing used by the filter (field width / points), cm
    D{i}.corner = corners.(cfg)(i, :);
    D{i}.mask   = subdomain_mask(D{i}.px, D{i}.corner);   % logical n x n, true inside the subdomain
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
            acc = acc + mean(Q(S.mask), 'omitnan');            % average over the subdomain
        end
        Pi(i, iL) = acc/nf;
    end
end

figure('Name', 'Energy flux');
plot(Llist, Pi', 'o-', 'LineWidth', 1.2); hold on; yline(0, 'k:');
xlabel('L (cm)'); ylabel('\Pi(L)  (cm^2 s^{-3})');
legend(labels, 'Location', 'best');
title(sprintf('%s background, flux averaged over the analysis subdomain and %d frames', cfg, nf));
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
        th_all      = [th_all; th(S.mask)/180];                         % subdomain only, units of pi
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
title(sprintf('%s background, L = %g cm, analysis subdomain', cfg, Lpdf));

%% (3) fields at one frame: local flux, strain and stress eigenstructure
% columns: Pi^(L)(x,y) | rate of strain: lambda_max + extensional eigenvector
%          | sub-grid stress: lambda_max + eigenvector
% Full field shown, analysis subdomain dashed.  Eigenvectors are axial
% (sign-free), so they are drawn as short bars.  Image axes (y downwards),
% as in the raw frames of 3_photos/.
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
    imagesc(S.px(1,:), S.py(:,1), Q); axis ij equal tight; hold on
    cl = robust_max(Q); clim([-cl cl]); colormap(ax, cmap_div); colorbar
    draw_box(S.corner);
    title(sprintf('\\Pi^{(L)}(x,y)  (cm^2 s^{-3}),  box mean = %+.2e', mean(Q(S.mask), 'omitnan')));
    xlabel('x (cm)'); ylabel('y (cm)')

    ax = nexttile;                                          % rate of strain
    imagesc(S.px(1,:), S.py(:,1), lam_s); axis ij equal tight; hold on
    clim([0 robust_max(lam_s)]); colormap(ax, parula); colorbar
    draw_bars(S.px, S.py, ex_s, ey_s, 6, 0.8); draw_box(S.corner);
    title('rate of strain: \lambda_{max} (s^{-1}), extensional eigenvector'); xlabel('x (cm)')

    ax = nexttile;                                          % sub-grid stress
    imagesc(S.px(1,:), S.py(:,1), lam_t); axis ij equal tight; hold on
    clim([0 robust_max(lam_t)]); colormap(ax, parula); colorbar
    draw_bars(S.px, S.py, ex_t, ey_t, 6, 0.8); draw_box(S.corner);
    title('sub-grid stress: \lambda_{max} (cm^2 s^{-2}), eigenvector'); xlabel('x (cm)')
end

%% console summary
fprintf('Pi(L) [cm^2/s^3] over the analysis subdomain, columns = L, %s:\n', strjoin(tags, ', '));
disp([Llist(:) Pi']);
for i = 1:nset
    fprintf('%-18s  %4d grid points in subdomain   Pi(L=%g) = %+.3e   mean theta/pi = %.3f\n', ...
            tags{i}, nnz(D{i}.mask), Lpdf, Pi(i, Llist == Lpdf), mean(theta{i}));
end

%% local helpers
function mask = subdomain_mask(px, corner)
% logical mask of the analysis subdomain: grid points inside (or on) the
% quadrilateral with corners [x1 y1 ... x4 y4] (cm).  Same construction as
% the original analysis script (square grid: px(1,:) serves both axes).
    v = px(1, :);
    [X, Y] = meshgrid(v, v);
    P = reshape(corner, 2, 4)';
    [in, on] = inpolygon(X, Y, P(:,1), P(:,2));
    mask = in | on;
end

function draw_box(corner)
% outline of the analysis subdomain
    P = reshape(corner, 2, 4)';
    plot(P([1:4 1], 1), P([1:4 1], 2), 'k--', 'LineWidth', 1.2);
end

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
