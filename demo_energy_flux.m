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
% 20 cm x 20 cm field, but all averages and statistics are taken over an
% analysis subdomain (dashed box in the figures) that stays inside the zone
% swept by the reciprocating rod array.  The paper used a quadrilateral per
% set; here each is replaced by an axis-aligned rectangle with the same centre
% and area (a square where the rod array is rotated by 45 deg), listed below
% in cm (image axes).  The field plots are cropped to the neighbourhood of
% the subdomain.  The 21 archived frames replace the full-length records, so
% the numbers illustrate the method rather than reproduce the figure captions
% exactly.

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

% analysis subdomain of each set: rectangle [x_min x_max y_min y_max] in cm
% (same order as the rows of `tags`)
box = struct( ...
    'cellular', [ 6.5 12.9   7.8 11.5;      % set21
                  8.1 11.5   7.5 12.2;      % set27
                  7.1 12.4   7.9 11.1;      % set32
                  8.4 12.5   8.1 12.2], ... % set39 (array at 45 deg: square)
    'shear',    [ 8.7 14.8   7.8 11.7;      % set24
                  5.1  9.5   7.6 12.0;      % set30 (array at 45 deg: square)
                  5.1  9.5   7.6 12.0;      % set39 (array at 45 deg: square)
                  9.5 15.5   7.8 11.7]);    % set48
margin = 2;                               % field plots: box widened by this much on each side (cm)

Cfg  = [upper(cfg(1)) cfg(2:end)];
nset = numel(tags);

D = cell(1, nset);
for i = 1:nset
    file = sprintf('2_data/%s/%s_%s_set%d.mat', cfg, Cfg, tags{i}, setno.(cfg)(i));
    D{i} = load(file);
    D{i}.file   = file;
    D{i}.d0     = 20/size(D{i}.px, 1);   % grid spacing used by the filter (field width / points), cm
    D{i}.box    = box.(cfg)(i, :);
    D{i}.mask   = subdomain_mask(D{i}.px, D{i}.py, D{i}.box);   % logical n x n, true inside the subdomain
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
% Cropped to the analysis subdomain (dashed) plus `margin`.  Eigenvectors are
% axial (sign-free), so they are drawn as short bars.  Image axes (y
% downwards), as in the raw frames of 3_photos/.
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
    cl = robust_max(Q(S.mask)); clim([-cl cl]); colormap(ax, cmap_div); colorbar
    draw_box(S.box); crop_to(S.box, margin);
    title(sprintf('\\Pi^{(L)}(x,y)  (cm^2 s^{-3}),  box mean = %+.2e', mean(Q(S.mask), 'omitnan')));
    xlabel('x (cm)'); ylabel('y (cm)')

    ax = nexttile;                                          % rate of strain
    imagesc(S.px(1,:), S.py(:,1), lam_s); axis ij equal tight; hold on
    clim([0 robust_max(lam_s(S.mask))]); colormap(ax, parula); colorbar
    draw_bars(S.px, S.py, ex_s, ey_s, 3, 0.4); draw_box(S.box); crop_to(S.box, margin);
    title('rate of strain: \lambda_{max} (s^{-1}), extensional eigenvector'); xlabel('x (cm)')

    ax = nexttile;                                          % sub-grid stress
    imagesc(S.px(1,:), S.py(:,1), lam_t); axis ij equal tight; hold on
    clim([0 robust_max(lam_t(S.mask))]); colormap(ax, parula); colorbar
    draw_bars(S.px, S.py, ex_t, ey_t, 3, 0.4); draw_box(S.box); crop_to(S.box, margin);
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
function mask = subdomain_mask(px, py, b)
% logical mask of the analysis subdomain: grid points inside (or on the edge
% of) the rectangle b = [x_min x_max y_min y_max] (cm)
    mask = px >= b(1) & px <= b(2) & py >= b(3) & py <= b(4);
end

function draw_box(b)
% outline of the analysis subdomain
    plot(b([1 2 2 1 1]), b([3 3 4 4 3]), 'k--', 'LineWidth', 1.2);
end

function crop_to(b, margin)
% restrict the axes to the subdomain plus a margin (cm)
    xlim([b(1) - margin, b(2) + margin]); ylim([b(3) - margin, b(4) + margin]);
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
