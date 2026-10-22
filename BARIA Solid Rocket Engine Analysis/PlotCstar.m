%% c* per-test scatter — standalone figure
%  Run from main3.m (needs c_vec, c, Sigma_c, Is in workspace)

col_LP  = [0.2157, 0.4941, 0.7216];
col_MP  = [0.3020, 0.6863, 0.2902];
col_HP  = [0.8941, 0.1020, 0.1098];
col_fit = [0.10,   0.10,   0.10  ];
col_bnd = [0.72,   0.72,   0.72  ];

fs = 12;
lw = 1.8;
ms = 70;

figure('Name','c* per Test','Position',[200,200,700,430]);
hold on; box on; grid on;
set(gca,'GridAlpha',0.18,'GridLineStyle',':','FontSize',fs);

scatter(1:3:27, c_vec(1:3:27), ms, col_LP,'o','filled','MarkerFaceAlpha',0.88,'DisplayName','LP');
scatter(2:3:27, c_vec(2:3:27), ms, col_MP,'s','filled','MarkerFaceAlpha',0.88,'DisplayName','MP');
scatter(3:3:27, c_vec(3:3:27), ms, col_HP,'^','filled','MarkerFaceAlpha',0.88,'DisplayName','HP');

patch([0.5 27.5 27.5 0.5], [c-Sigma_c c-Sigma_c c+Sigma_c c+Sigma_c], ...
      col_fit,'FaceAlpha',0.06,'EdgeColor','none','HandleVisibility','off');
yline(c,'-','Color',col_fit,'LineWidth',lw, ...
      'Label',sprintf(' c* mean = %.0f m/s',c), ...
      'LabelHorizontalAlignment','right','FontSize',fs-1);
yline(c+Sigma_c,'--','Color',col_bnd,'LineWidth',1.3, ...
      'Label','+\sigma','LabelHorizontalAlignment','left','FontSize',fs-1);
yline(c-Sigma_c,'--','Color',col_bnd,'LineWidth',1.3, ...
      'Label','-\sigma','LabelHorizontalAlignment','left','FontSize',fs-1);

legend('Location','southeast','FontSize',fs-1,'NumColumns',3);
xlabel('Test index  [-]','FontSize',fs);
ylabel('c^*  [m/s]','FontSize',fs);
title(sprintf('Characteristic Velocity c*', ...
              c), 'FontSize',fs,'FontWeight','bold','Interpreter','latex');
xlim([0 28]);
