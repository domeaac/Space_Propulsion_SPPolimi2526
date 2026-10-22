%% ================================================================
%  FINAL SUMMARY PLOTS — BARIA Bates Motor  (v2)
%  Run from main3.m after all sub-scripts have populated the workspace.
%
%  Figures produced:
%   1 — Experimental pressure traces (3×3 batch grid)
%   2 — Vieille law: true log-log axes + linear view
%   3 — c*: per-test scatter | violin+IQR+strip | histogram
%   4 — Motor simulation: p(t), Ab(t), rb(t), e(t),
%                         Ab–p phase portrait, grain cross-section
%% ================================================================

%— Shared palette —%
col_LP  = [0.2157, 0.4941, 0.7216];   % blue
col_MP  = [0.3020, 0.6863, 0.2902];   % green
col_HP  = [0.8941, 0.1020, 0.1098];   % red
col_fit = [0.10,   0.10,   0.10  ];   % near-black
col_bnd = [0.72,   0.72,   0.72  ];   % grey (uncertainty band)
col_g   = {col_LP, col_MP, col_HP};

fs = 11;   % base font size
lw = 1.8;  % base line width
ms = 60;   % scatter marker area [pt²]


%% ----------------------------------------------------------------
%  FIGURE 1 — Experimental pressure traces  (3 × 3 batch grid)
%% ----------------------------------------------------------------
figure('Name','Pressure Traces','Position',[30,30,1150,750]);
tl1 = tiledlayout(3,3,'TileSpacing','tight','Padding','compact');
title(tl1,'Experimental Pressure Traces — All 9 Batches', ...
     'FontSize',fs+2,'FontWeight','bold');
xlabel(tl1,'Time  t  [s]',     'FontSize',fs);
ylabel(tl1,'Pressure  p  [bar]','FontSize',fs);

for i = 1:9
    p_raw = data.(batch{i});
    t_i   = (0:size(p_raw,1)-1) * 1e-3;
    if i == 7                             % channel swap for batch 7
        p_L = p_raw(:,1); p_M = p_raw(:,3); p_H = p_raw(:,2);
    else
        p_L = p_raw(:,1); p_M = p_raw(:,2); p_H = p_raw(:,3);
    end

    nexttile;
    hold on;
    plot(t_i, p_L, 'Color', col_LP, 'LineWidth', 1.0);
    plot(t_i, p_M, 'Color', col_MP, 'LineWidth', 1.0);
    plot(t_i, p_H, 'Color', col_HP, 'LineWidth', 1.0);
    hold off;
    box on;  grid on;
    set(gca,'GridAlpha',0.18,'GridLineStyle',':','FontSize',fs-2);
    title(sprintf('Batch %d',i),'FontSize',fs-1,'FontWeight','bold');
    if i == 1
        legend({'LP','MP','HP'},'Location','northeast','FontSize',fs-3,'Box','on');
    end
end


%% ----------------------------------------------------------------
%  FIGURE 2 — Vieille Law: true log-log axes  +  linear scale
%  (using XScale/YScale log gives real pressure tick labels like
%   30, 40, 60 bar — far more readable than ln-transformed axes)
%% ----------------------------------------------------------------
p_pl   = logspace(log10(min(p_val)*0.84), log10(max(p_val)*1.18), 300);
rb_fit = a .* p_pl .^ n;
rb_hi  = a_new_up   .* p_pl .^ (n + Inc_n);
rb_lo  = a_new_down .* p_pl .^ (n - Inc_n);

figure('Name','Vieille Law','Position',[50,50,1050,490]);
tl2 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(tl2, sprintf('Vieille Burning-Rate Law   r_b = a\\cdotp^n     a = %.4f mm/s/bar^n,  n = %.4f,  R^2 = %.5f', ...
                    a, n, R2), 'FontSize',fs,'FontWeight','bold');

% ---- Left: log-log with true scales, each batch a distinct colour ----
nexttile;
hold on; box on; grid on;
set(gca,'XScale','log','YScale','log','GridAlpha',0.18,'GridLineStyle',':','FontSize',fs);

cmap9 = lines(9);
for i = 1:9
    scatter(p_eff_vec(i,1),rb_vec(i,1), ms*0.85, cmap9(i,:),'o','filled','MarkerFaceAlpha',0.80);
    scatter(p_eff_vec(i,2),rb_vec(i,2), ms*0.85, cmap9(i,:),'s','filled','MarkerFaceAlpha',0.80);
    scatter(p_eff_vec(i,3),rb_vec(i,3), ms*0.85, cmap9(i,:),'^','filled','MarkerFaceAlpha',0.80);
end

fill([p_pl,fliplr(p_pl)],[rb_hi,fliplr(rb_lo)], col_bnd, ...
     'FaceAlpha',0.28,'EdgeColor','none');
plot(p_pl, rb_fit, '-', 'Color',col_fit, 'LineWidth',lw+0.5);

% Symbol-only legend entries for pressure level
hLP = scatter(nan,nan,ms,col_LP,'o','filled','DisplayName','LP');
hMP = scatter(nan,nan,ms,col_MP,'s','filled','DisplayName','MP');
hHP = scatter(nan,nan,ms,col_HP,'^','filled','DisplayName','HP');
hFt = plot(nan,nan,'-','Color',col_fit,'LineWidth',lw,'DisplayName','Power-law fit');
hBd = fill(nan,nan,col_bnd,'FaceAlpha',0.28,'EdgeColor','none','DisplayName','\pm1\sigma band');
legend([hLP,hMP,hHP,hFt,hBd],'Location','northwest','FontSize',fs-2);

text(0.97, 0.04, ...
     sprintf('n = %.4f \\pm %.4f\na = %.4f \\pm %.4f\nR^2 = %.5f', n,Inc_n,a,Inc_a,R2), ...
     'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
     'FontSize',fs-1,'BackgroundColor','w','EdgeColor',[0.72 0.72 0.72]);

xlabel('p_{eff}  [bar]','FontSize',fs);
ylabel('r_b  [mm/s]',   'FontSize',fs);
title('Log-Log  (true scale)','FontSize',fs,'FontWeight','bold');

% ---- Right: linear scale with per-level mean ± std ----
nexttile;
hold on; box on; grid on;
set(gca,'GridAlpha',0.18,'GridLineStyle',':','FontSize',fs);

scatter(p_eff_vec(:,1),rb_vec(:,1), ms, col_LP,'o','filled','MarkerFaceAlpha',0.55);
scatter(p_eff_vec(:,2),rb_vec(:,2), ms, col_MP,'s','filled','MarkerFaceAlpha',0.55);
scatter(p_eff_vec(:,3),rb_vec(:,3), ms, col_HP,'^','filled','MarkerFaceAlpha',0.55);

fill([p_pl,fliplr(p_pl)],[rb_hi,fliplr(rb_lo)], col_bnd, ...
     'FaceAlpha',0.28,'EdgeColor','none');
plot(p_pl, rb_fit, '-', 'Color',col_fit, 'LineWidth',lw+0.5);

for k = 1:3
    errorbar(mean(p_eff_vec(:,k)), mean(rb_vec(:,k)), std(rb_vec(:,k)), ...
             'o','Color',col_g{k},'MarkerFaceColor','w', ...
             'MarkerSize',9,'LineWidth',lw,'CapSize',7);
end

xlabel('p_{eff}  [bar]','FontSize',fs);
ylabel('r_b  [mm/s]',   'FontSize',fs);
title('Linear  (circles = level mean ± std)','FontSize',fs,'FontWeight','bold');


%% ----------------------------------------------------------------
%  FIGURE 3 — c*  statistics
%  Three panels: per-test scatter | custom violin + strip | histogram
%  Violin uses ksdensity (Statistics & ML Toolbox)
%% ----------------------------------------------------------------
figure('Name','c* Statistics','Position',[70,70,1100,470]);
tl3 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
title(tl3, sprintf('Characteristic Velocity  c^*      mean c* = %.0f m/s,   \\sigma = %.0f m/s,   I_s = %.2f s', ...
                   c, Sigma_c, Is), 'FontSize',fs,'FontWeight','bold');

groups   = {c1, c2, c3};
labels_g = {'LP','MP','HP'};

% ---- Left: per-test scatter with ±σ band ----
nexttile;
hold on; box on; grid on;
set(gca,'GridAlpha',0.18,'GridLineStyle',':','FontSize',fs);

scatter(1:3:27, c_vec(1:3:27), ms, col_LP,'o','filled','MarkerFaceAlpha',0.88,'DisplayName','LP');
scatter(2:3:27, c_vec(2:3:27), ms, col_MP,'s','filled','MarkerFaceAlpha',0.88,'DisplayName','MP');
scatter(3:3:27, c_vec(3:3:27), ms, col_HP,'^','filled','MarkerFaceAlpha',0.88,'DisplayName','HP');

% Shaded ±1σ band
patch([0.5 27.5 27.5 0.5], [c-Sigma_c c-Sigma_c c+Sigma_c c+Sigma_c], ...
      col_fit, 'FaceAlpha',0.06,'EdgeColor','none','HandleVisibility','off');
yline(c,'-','Color',col_fit,'LineWidth',lw, ...
      'Label',sprintf(' mean c*= %.0f m/s',c), ...
      'LabelHorizontalAlignment','right','FontSize',fs-1);
yline(c+Sigma_c,'--','Color',col_bnd,'LineWidth',1.3, ...
      'Label','+\sigma','LabelHorizontalAlignment','left','FontSize',fs-2);
yline(c-Sigma_c,'--','Color',col_bnd,'LineWidth',1.3, ...
      'Label','-\sigma','LabelHorizontalAlignment','left','FontSize',fs-2);

legend('Location','south','FontSize',fs-2,'NumColumns',3);
xlabel('Test index','FontSize',fs);
ylabel('c^*  [m/s]','FontSize',fs);
title('c^* per Test','FontSize',fs,'FontWeight','bold');
xlim([0 28]);

% ---- Centre: custom violin + IQR box + jittered points ----
nexttile;
hold on; box on; grid on;
set(gca,'XTick',1:3,'XTickLabel',labels_g,'GridAlpha',0.18,'GridLineStyle',':','FontSize',fs);

rng(42);
vw = 0.38;   % violin half-width in x-axis units
for k = 1:3
    d   = groups{k};
    cx  = k;
    col = col_g{k};

    % KDE — let MATLAB choose bandwidth automatically (Silverman's rule)
    [fk, xik] = ksdensity(d, 'NumPoints',150);
    fn  = fk / max(fk) * vw;

    % Filled violin body
    fill([cx+fn, cx-fliplr(fn)], [xik, fliplr(xik)], col, ...
         'FaceAlpha',0.35, 'EdgeColor',col*0.60, 'LineWidth',1.2);

    % IQR box overlaid on the violin
    q25 = quantile(d, 0.25);
    q75 = quantile(d, 0.75);
    med = median(d);
    rectangle('Position',[cx-0.055, q25, 0.11, q75-q25], ...
              'FaceColor','w','EdgeColor',col*0.55,'LineWidth',2.0);
    plot([cx-0.06, cx+0.06],[med, med],'-','Color',col*0.45,'LineWidth',3.5);

    % Jittered strip chart
    jx = cx + (rand(9,1)-0.5)*0.16;
    scatter(jx, d, ms*0.85, col, 'filled', ...
            'MarkerFaceAlpha',0.80,'MarkerEdgeColor','w','LineWidth',0.6);
end

yline(c,'--','Color',col_fit,'LineWidth',1.3, ...
      'Label',sprintf(' mean c* = %.0f',c), ...
      'LabelHorizontalAlignment','right','FontSize',fs-1);

ylabel('c^*  [m/s]','FontSize',fs);
title('Violin  +  IQR  +  Data','FontSize',fs,'FontWeight','bold');
xlim([0.45 3.55]);

% ---- Right: stacked histogram + Gaussian overlay ----
nexttile;
hold on; box on; grid on;
set(gca,'GridAlpha',0.18,'GridLineStyle',':','FontSize',fs);

edges_h = linspace(min(c_vec)-45, max(c_vec)+45, 14);
for k = 1:3
    histogram(groups{k}, edges_h, ...
              'FaceColor',col_g{k}, 'FaceAlpha',0.48, ...
              'EdgeColor','w', 'LineWidth',0.7, 'DisplayName',labels_g{k});
end

xg  = linspace(min(c_vec)-90, max(c_vec)+90, 500);
scl = numel(c_vec) * (edges_h(2)-edges_h(1));
plot(xg, normpdf(xg,c,Sigma_c)*scl, '-', 'Color',col_fit, 'LineWidth',lw+0.3, ...
     'DisplayName',sprintf('\\mathcal{N}(%.0f, %.0f)', c, Sigma_c));

xl = xline(c, '-', ' mean c*', 'Color',col_fit, 'LineWidth',1.3, ...
           'LabelOrientation','horizontal','FontSize',fs-2);
xl.Alpha = 0.65;

legend('Location','northwest','FontSize',fs-2);
xlabel('c^*  [m/s]','FontSize',fs);
ylabel('Count','FontSize',fs);
title('Distribution by Level','FontSize',fs,'FontWeight','bold');
xlim([min(c_vec)-65, max(c_vec)+65]);


%% ----------------------------------------------------------------
%  FIGURE 4 — Motor simulation  (2×3 panel)
%  Panels: p(t) | Ab(t) | rb(t)
%          e(t) | Ab–p phase portrait (spanning 2 cols)
%% ----------------------------------------------------------------
Ab_t   = Ab(z_vec);
e_t    = a .* z_vec;           % web distance burned  [mm]
rb_t   = a .* (p_sol .^ n);   % instantaneous rb     [mm/s]

idx_bo = find(Ab_t > 1, 1, 'last');
if isempty(idx_bo), idx_bo = length(t_vec); end

figure('Name','Motor Simulation','Position',[100,100,1100,600]);
tl4 = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
title(tl4,'Motor Simulation — Mid Pressure (MP) Case', ...
     'FontSize',fs+2,'FontWeight','bold');

% ---- p(t) ----
nexttile;
hold on; box on; grid on;
set(gca,'GridAlpha',0.20,'GridLineStyle',':','FontSize',fs);
patch([t_vec(1),t_vec(idx_bo),t_vec(idx_bo),t_vec(1)], ...
      [0,0,max(p_sol)*1.10,max(p_sol)*1.10], ...
      col_HP,'FaceAlpha',0.06,'EdgeColor','none');
plot(t_vec, p_sol, 'Color',col_HP, 'LineWidth',lw);
yline(p_eff_MP,'--','Color',col_MP,'LineWidth',1.4, ...
      'Label',sprintf('\\bar{p}_{exp,MP} = %.1f bar',p_eff_MP), ...
      'LabelHorizontalAlignment','right','FontSize',fs-1);
xlabel('t  [s]','FontSize',fs); ylabel('p  [bar]','FontSize',fs);
title('Chamber Pressure','FontSize',fs,'FontWeight','bold');

% ---- Ab(t) ----
nexttile;
hold on; box on; grid on;
set(gca,'GridAlpha',0.20,'GridLineStyle',':','FontSize',fs);
plot(t_vec, Ab_t, 'Color',col_LP, 'LineWidth',lw);
yline(Ab0,'--k','LineWidth',1.2, ...
      'Label',sprintf('A_{b,0} = %.0f mm^2',Ab0), ...
      'LabelHorizontalAlignment','right','FontSize',fs-1);
xlabel('t  [s]','FontSize',fs); ylabel('A_b  [mm^2]','FontSize',fs);
title('Burning Area','FontSize',fs,'FontWeight','bold');

% ---- rb(t) ----
nexttile;
hold on; box on; grid on;
set(gca,'GridAlpha',0.20,'GridLineStyle',':','FontSize',fs);
plot(t_vec, rb_t, 'Color',col_fit, 'LineWidth',lw);
yline(rb_MP,'--','Color',col_MP,'LineWidth',1.4, ...
      'Label',sprintf('\\bar{r}_{b,exp,MP} = %.2f mm/s',rb_MP), ...
      'LabelHorizontalAlignment','right','FontSize',fs-1);
xlabel('t  [s]','FontSize',fs); ylabel('r_b  [mm/s]','FontSize',fs);
title('Burning Rate','FontSize',fs,'FontWeight','bold');

% ---- e(t) — web burn distance ----
nexttile;
hold on; box on; grid on;
set(gca,'GridAlpha',0.20,'GridLineStyle',':','FontSize',fs);
plot(t_vec, e_t, 'Color',col_MP, 'LineWidth',lw);
yline(wt,'--k','LineWidth',1.3, ...
      'Label',sprintf('Burnout  (e = %d mm)',wt), ...
      'LabelHorizontalAlignment','right','FontSize',fs-1);
xlabel('t  [s]','FontSize',fs); ylabel('e  [mm]','FontSize',fs);
title('Web Burn Distance  e = a\cdotz','FontSize',fs,'FontWeight','bold');

% ---- Ab–p phase portrait, trajectory coloured by time (spans 2 columns) ----
nexttile([1,2]);
hold on; box on; grid on;
set(gca,'GridAlpha',0.20,'GridLineStyle',':','FontSize',fs);

n_seg   = idx_bo - 1;
cmap_t  = parula(n_seg);
for k = 1:n_seg
    plot(p_sol(k:k+1), Ab_t(k:k+1), '-', 'Color',cmap_t(k,:), 'LineWidth',2.2);
end
colormap(gca, parula);
cb = colorbar;
cb.Label.String = 'Time  [s]';
cb.Label.FontSize = fs-1;
clim([0, t_vec(idx_bo)]);

scatter(p_sol(1),      Ab_t(1),      ms*1.3, [0.1 0.7 0.1], '^','filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7,'DisplayName','Ignition');
scatter(p_sol(idx_bo), Ab_t(idx_bo), ms*1.3, [0.8 0.1 0.1], 'v','filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7,'DisplayName','Burnout');
legend('Location','northeast','FontSize',fs-2);

xlabel('p  [bar]','FontSize',fs); ylabel('A_b  [mm^2]','FontSize',fs);
title('Phase Portrait  A_b(p)','FontSize',fs,'FontWeight','bold');
