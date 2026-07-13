    % run_STE_experiment.m
    % =====================================================================
    % INTEGRATED STE PIPELINE — Tasks 1 + 2 + 3.
    % Outputs saved to D:\Mtech\SEM3\MTP\STE_project\figs\images\after_lawnmower
    % =====================================================================
    
    clear; clc; close all;
    rng(42);
    
    figDir = 'D:\Mtech\SEM3\MTP\STE_project\figs\images\after_lawnmower';
    if ~exist(figDir, 'dir'); mkdir(figDir); end
    
    % ---- physics config (ground truth) ----------------------------------
    Nx = 60; Ny = 60;
    dx = 1.0; dy = 1.0;
    dt = 0.1;
    u  = 1.0; v = 0.0;
    D  = 0.5;
    Q  = 1.0;
    src_i = 10; src_j = 30;
    
    ad = AdvectionDiffusion2D(Nx, Ny, dx, dy, dt, u, v, D, src_i, src_j);
    ad.boundary_type = 'absorbing';   % zero-gradient walls -- plume flows out,
                                       % no artificial Dirichlet sink at the edge
    warmupSteps = 600;
    for n = 1:warmupSteps
        ad = ad.step(Q);   % must use step(), not step_loop() -- step_loop()
                            % has no boundary_type branch and leaves edge cells
                            % pinned at their initial value of zero regardless
    end
    
    % ---- TASK 2: 10x10 uniform lattice, [2,57] both axes, snake order ----
    margin = 2;
    xLat = linspace(margin, (Nx-1)*dx - margin, 10);
    yLat = linspace(margin, (Ny-1)*dy - margin, 10);
    sensorPath = zeros(100, 2);
    idx = 1;
    for jy = 1:10
        if mod(jy,2) == 1, xseq = 1:10; else, xseq = 10:-1:1; end
        for ix = xseq
            sensorPath(idx,:) = [xLat(ix), yLat(jy)];
            idx = idx + 1;
        end
    end
    nIter = size(sensorPath, 1);
    
    % ---- TASK 1: particle filter over (x_s, y_s, Q) ---------------------
    N       = 1000;
    theta_Q = 0.25;
    sigma_v = 0.01;
    pf = ParticleFilterSIR(N, 0, Nx, 0, Ny, theta_Q);
    
    % ---- TASK 3 setup: true level-set extracted ONCE --------------------
    lse      = LevelSetExtractor(dx, dy);
    Cth_true = 0.1 * max(ad.C(:));
    trueCont = lse.extract(ad.C, Cth_true);
    psTrue   = contourToPolyshape(trueCont);
    if psTrue.NumRegions == 0 || area(psTrue) <= 0
        error('True level-set empty/degenerate at Cth=%.4g.', Cth_true);
    end
    
    xPhys = (0:Nx-1) * dx;
    yPhys = (0:Ny-1) * dy;
    [XI, YJ] = ndgrid(xPhys, yPhys);
    
    ratioWarn = 1.0;   % flags any estimate larger than the true level-set area
    iouWarn   = 0.5;
    
    ESS_hist=zeros(nIter,1);  H_hist=zeros(nIter,1);
    IoU_hist=nan(nIter,1);     Coverage_hist=nan(nIter,1);    ratio_hist=nan(nIter,1);
    warn_hist=false(nIter,1); resamp_hist=false(nIter,1);
    theta_hist=nan(nIter,3);
    particles_hist = zeros(N, 3, nIter);
    weights_hist   = zeros(N, nIter);
    estCont_hist   = cell(nIter, 1);
    
    fprintf('%4s %8s %8s %8s %8s %8s %7s  %s\n', ...
        'iter','ESS','H','IoU','Coverage','ratio','resamp','warn');
    
    for k = 1:nIter
        sx = sensorPath(k,1);  sy = sensorPath(k,2);
        sensor = PointSensor(sx, sy, sigma_v);
        y_k = sensor.measure(ad);
        pf.update(y_k, sensor, u, D);
    
        ESS_hist(k) = pf.computeESS();
        H_hist(k)   = pf.computeEntropy();
    
        particles_hist(:,:,k) = pf.particles;
        weights_hist(:,k)     = pf.weights;
    
        [~, ord] = sort(pf.weights, 'descend');
        top = ord(1:5);
        w5  = pf.weights(top);  w5 = w5 / sum(w5);
        theta_hat = w5' * pf.particles(top, :);
        theta_hist(k,:) = theta_hat;
    
        estField = arrayfun(@(xx,yy) ...
            plumeConc(theta_hat(1), theta_hat(2), xx, yy, u, D, theta_hat(3)), ...
            XI, YJ);
        Cth_est = 0.1 * max(estField(:));

        if Cth_est <= 0 || ~isfinite(Cth_est)
            warn_hist(k) = true;
        else
            estCont = lse.extract(estField, Cth_est);
            estCont_hist{k} = estCont;
            psEst = contourToPolyshape(estCont);
            if psEst.NumRegions > 0 && area(psEst) > 0
                inter = area(intersect(psTrue, psEst));
                uni   = area(union(psTrue, psEst));
                IoU_hist(k)      = inter / uni;
                Coverage_hist(k) = inter / area(psTrue);
                ratio_hist(k) = area(psEst) / area(psTrue);
                warn_hist(k)  = (ratio_hist(k) > ratioWarn) && (IoU_hist(k) < iouWarn);
            else
                warn_hist(k) = true;
            end
        end
    
        if ESS_hist(k) < N/2
            pf.resample();  pf.roughen();  resamp_hist(k) = true;
        end
    
        fprintf('%4d %8.1f %8.4f %8.4f %8.4f %8.4f %7d  %s\n', ...
            k, ESS_hist(k), H_hist(k), IoU_hist(k), Coverage_hist(k), ratio_hist(k), ...
            resamp_hist(k), tern(warn_hist(k),'OVER-COVER',''));
    end
    
    pathColor = [1 0.5 0];   % orange -- confirmed legible, unchanged
    
    % ===================== FIGURE 1: field + path =========================
    figure;
    imagesc(xPhys, yPhys, ad.C'); set(gca,'YDir','normal');
    axis equal tight; colorbar; hold on;
    hSrc  = plot(src_i, src_j, 'r*', 'MarkerSize',14, 'LineWidth',2);
    hPath = plot(sensorPath(:,1), sensorPath(:,2), '.-', ...
        'Color', pathColor, 'MarkerSize',10, 'LineWidth',1.2);
    xlabel('x'); ylabel('y');
    legend([hSrc hPath], {'true source','sensor path'}, 'Location','northeastoutside');
    hold off;
    exportgraphics(gcf, fullfile(figDir,'field_and_path.pdf'), 'ContentType','vector');
    
    % ===================== FIGURE 2: IoU + Coverage + ratio history ========
    figure; tiledlayout(2,1);
    nexttile; hold on; grid on; ylim([0 1]);
    hU = plot(1:nIter, IoU_hist, 'm.-');
    hC = plot(1:nIter, Coverage_hist, 'b.-');
    ylabel('overlap metric');
    rIdx = find(resamp_hist);
    if isempty(rIdx)
        legend([hU hC], {'IoU','Coverage'}, 'Location','best');
    else
        hR = xline(rIdx, ':', 'Color',[.6 .6 .6]);
        legend([hU hC hR(1)], {'IoU','Coverage','resample'}, 'Location','best');
    end
    nexttile; hold on; grid on;
    hRt = plot(1:nIter, ratio_hist, 'k.-');
    hWl = yline(ratioWarn, 'r--');
    wIdx = find(warn_hist & isfinite(ratio_hist));
    hW = plot(wIdx, ratio_hist(wIdx), 'ro', 'MarkerSize',8, 'LineWidth',1.2);
    xlabel('iteration'); ylabel('area ratio  est/true');
    legend([hRt hWl hW], {'area ratio','warn threshold','over-coverage fired'}, ...
        'Location','best');
    exportgraphics(gcf, fullfile(figDir,'iou_ratio_history.pdf'), 'ContentType','vector');
    
    % ===================== OPTION C: spaghetti, fixed for legibility ======
    % FIX: early (unconverged) contours were visually dominating over the
    % converging trend. Now faded/thin early, bold/opaque late -- eye is
    % drawn to the convergence, not the early noise.
    figure; hold on;
    for k = 1:nIter
        ec = estCont_hist{k};
        if isempty(ec), continue; end
        frac  = k / nIter;                    % 0 (early) -> 1 (late)
        alpha = 0.08 + 0.75*frac;
        lw    = 0.5  + 2.0*frac;
        col   = [0.85 0.85 0.85] * (1-frac) + [0.85 0.1 0.1] * frac;  % grey->red
        for c = 1:numel(ec)
            plot(ec(c).x, ec(c).y, '-', 'Color', [col alpha], 'LineWidth', lw);
        end
    end
    hTrue = plot(psTrue, 'FaceColor','none', 'EdgeColor','b', 'LineWidth', 2.5);
    hSrc2 = plot(src_i, src_j, 'k*', 'MarkerSize',16, 'LineWidth',2);
    axis equal; xlim([0 Nx*dx]); ylim([0 Ny*dy]); xlabel('x'); ylabel('y');
    legend([hTrue hSrc2], ...
        {'true level-set','true source'}, 'Location','northeastoutside');
    % manual fade legend note via text, since colour here encodes iteration
    % progress rather than a single series
    text(2, 3, 'grey \rightarrow red : early \rightarrow late iteration', ...
        'FontSize', 9, 'BackgroundColor','w');
    hold off;
    exportgraphics(gcf, fullfile(figDir,'contour_spaghetti_convergence.pdf'), 'ContentType','vector');
    
    % ===================== GIF 1: level-set convergence only ==============
    gifFile1 = fullfile(figDir, 'levelset_convergence.gif');
    fig1 = figure('Position',[100 100 700 600]);
    for k = 1:nIter
        clf(fig1); hold on;
        imagesc(xPhys, yPhys, ad.C'); set(gca,'YDir','normal');
        axis equal tight; colormap(parula);
    
        hTrue = plot(psTrue, 'FaceColor','none', 'EdgeColor','b', 'LineWidth', 2);
        ec = estCont_hist{k};
        hEst = [];
        if ~isempty(ec)
            for c = 1:numel(ec)
                hEst = plot(ec(c).x, ec(c).y, 'r-', 'LineWidth', 2);
            end
        end
        hPath = plot(sensorPath(1:k,1), sensorPath(1:k,2), '.-', ...
            'Color', pathColor, 'MarkerSize', 10, 'LineWidth', 1.5);
        hCur = plot(sensorPath(k,1), sensorPath(k,2), 'o', ...
            'MarkerSize', 10, 'MarkerFaceColor', pathColor, 'MarkerEdgeColor','k');
        hSrc3 = plot(src_i, src_j, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
    
        xlim([0 Nx*dx]); ylim([0 Ny*dy]); xlabel('x'); ylabel('y');
        legH = [hTrue hSrc3 hPath hCur];
        legL = {'true level-set','true source','sensor path','current position'};
        if ~isempty(hEst)
            legH(end+1) = hEst; legL{end+1} = 'estimated level-set';
        end
        legend(legH, legL, 'Location','eastoutside', 'FontSize',8);
        text(2, Ny*dy-3, sprintf('iter %d/%d   cov=%.2f  ratio=%.2f', ...
            k, nIter, Coverage_hist(k), ratio_hist(k)), ...
            'Color','w', 'FontWeight','bold', 'BackgroundColor',[0 0 0 0.4]);
        hold off; drawnow;
    
        frame = getframe(fig1);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if k == 1
            imwrite(imind, cm, gifFile1, 'gif', 'Loopcount', inf, 'DelayTime', 0.15);
        else
            imwrite(imind, cm, gifFile1, 'gif', 'WriteMode', 'append', 'DelayTime', 0.15);
        end
    end
    close(fig1);
    
    % ===================== GIF 2: particle cloud only (size = weight) =====
    % Fixed colour, no colormap -- size is the only weight encoding, per
    % Dr. Nanavati's preference.
    gifFile2 = fullfile(figDir, 'particle_cloud_convergence.gif');
    fig2 = figure('Position',[100 100 700 600]);
    particleColor = [0.15 0.15 0.65];   % fixed, not weight-mapped
    for k = 1:nIter
        clf(fig2); hold on;
        imagesc(xPhys, yPhys, ad.C'); set(gca,'YDir','normal');
        axis equal tight; colormap(gray);   % neutral backdrop, doesn't compete
        caxis([0 max(ad.C(:))]);
    
        w = weights_hist(:,k);  wN = w / max(w);
        hPc = scatter(particles_hist(:,1,k), particles_hist(:,2,k), ...
            4 + 180*wN, particleColor, 'filled', 'MarkerFaceAlpha', 0.6);
    
        hPath = plot(sensorPath(1:k,1), sensorPath(1:k,2), '.-', ...
            'Color', pathColor, 'MarkerSize', 10, 'LineWidth', 1.5);
        hCur = plot(sensorPath(k,1), sensorPath(k,2), 'o', ...
            'MarkerSize', 10, 'MarkerFaceColor', pathColor, 'MarkerEdgeColor','k');
        hSrc4 = plot(src_i, src_j, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
    
        xlim([0 Nx*dx]); ylim([0 Ny*dy]); xlabel('x'); ylabel('y');
        legend([hPc hSrc4 hPath hCur], ...
            {'particle (size = weight)','true source','sensor path','current position'}, ...
            'Location','eastoutside', 'FontSize',8);
        text(2, Ny*dy-3, sprintf('iter %d/%d   ESS=%.0f', k, nIter, ESS_hist(k)), ...
            'Color','w', 'FontWeight','bold', 'BackgroundColor',[0 0 0 0.4]);
        hold off; drawnow;
    
        frame = getframe(fig2);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if k == 1
            imwrite(imind, cm, gifFile2, 'gif', 'Loopcount', inf, 'DelayTime', 0.15);
        else
            imwrite(imind, cm, gifFile2, 'gif', 'WriteMode', 'append', 'DelayTime', 0.15);
        end
    end
    close(fig2);
    
    fprintf('Level-set GIF:    %s\n', gifFile1);
    fprintf('Particle-cloud GIF: %s\n', gifFile2);
    
    % ===================== local helper functions =========================
    function ps = contourToPolyshape(contours)
        ps = polyshape();
        ws = warning('off','all');
        for c = 1:numel(contours)
            xc = contours(c).x(:);  yc = contours(c).y(:);
            if numel(xc) >= 3
                pc = polyshape(xc, yc, 'Simplify', true);
                if pc.NumRegions > 0, ps = union(ps, pc); end
            end
        end
        warning(ws);
    end
    
    function out = tern(c, a, b)
        if c, out = a; else, out = b; end
    end