%   Dist: nSta x nSta, 站点间距离矩阵（km）
%   rhoMat: nSta x nSta, Pearson 相关系数矩阵
%   dc: 第一变点距离（例如 300 km）

nSta = size(Dist, 1);
dc = d_thr; 

sigma = dc / 2;  


W = zeros(nSta);
for i = 1:nSta
    for j = i+1:nSta
        d_ij = Dist(i,j);
        if d_ij <= dc
            weight = abs(Rho(i,j)) * exp(-d_ij^2 / (2*sigma^2));
            W(i,j) = weight;
            W(j,i) = weight;
        end
    end
end


% W = W ./ (max(W(:)) + eps);  % 归一化处理

D = diag(sum(W));                   
D_sqrt_inv = D^(-0.5);
D_sqrt_inv(isinf(D_sqrt_inv)) = 0;   
L_sym = eye(nSta) - D_sqrt_inv * W * D_sqrt_inv;

m = min(20, nSta-1);  %前20个

opt = struct('issym', 1, 'tol', 1e-6, 'maxit', 1000);
[V, P] = eigs(L_sym, m, 'smallestabs', opt);  

figure('Position', [100, 100, 800, 600]);
plot(1:m, eigvals, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 8);
hold on;
xq = 1:m;
yq = smoothdata(eigvals, 'movmean', 5);
plot(xq, yq, 'r--', 'LineWidth', 1.5);
dy = diff(yq);                          
ddy = diff(dy);                         
[~, elbow1] = max(dy);
suggest_k1 = elbow1;
thresh = mean(dy) + std(dy);
dy_min = min(dy);
dy_max = max(dy);
dy_range = dy_max - dy_min;
yyaxis right;
plot(1.5:m, dy, 'g-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'g', 'MarkerSize', 6);
hold on;
plot([1.5, m], [thresh, thresh], 'm--', 'LineWidth', 2, 'Color', [0.8 0 0.8]);
x_fill = [1.5, m, m, 1.5];
y_fill = [thresh, thresh, dy_max*1.1, dy_max*1.1];
fill(x_fill, y_fill, [0.9 0.9 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');


jumps = find(dy > thresh);
if ~isempty(jumps)
    for i = 1:length(jumps)
        k_idx = jumps(i);
        plot(k_idx+0.5, dy(k_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    suggest_k2 = jumps(1);
else
    suggest_k2 = suggest_k1;
end


ylabel('First difference (dy)', 'FontSize', 12);
yyaxis left;
xlabel('Eigenvalue index k', 'FontSize', 12);
ylabel('Smallest eigenvalue of the normalized Laplacian', 'FontSize', 12);
% title('拉普拉斯谱的"特征间隙"分析与阈值检测', 'FontSize', 14);
title('Laplacian Eigengap Analysis and Thresholding', 'FontSize', 14);
grid on;
suggested_k = suggest_k2;


% xline(suggested_k, '--r', sprintf('建议 k=%d', suggested_k), ...
%     'LineWidth', 2, 'Color', [0.9 0.2 0.2], 'FontSize', 11);


info_str = {
    sprintf('Threshold = %.4f (mean+std)', thresh);
    sprintf('Optimal number of clusters k = %d', suggested_k);
%     sprintf('超阈值点: %s', mat2str(jumps));
%     sprintf('最大差分点: %d', elbow1)
};
annotation('textbox', [0.15, 0.75, 0.2, 0.15], ...
    'String', info_str, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', [1 1 0.8], ...
    'EdgeColor', [0.5 0.5 0.5], ...
    'FontSize', 10);

% legend({'原始特征值', '平滑曲线', '阈值线', '阈值区域', '超阈值点', ''}, ...
%     'Location', 'best', 'FontSize', 10);
legend({'Original eigenvalue', 'Smoothed curve', 'First derivative', 'Threshold line', 'Points above threshold', ''}, ...
    'Location', 'best', 'FontSize', 10);


xlim([0.5, m+0.5]);
yyaxis left
ylim([0, max(eigvals)*1.1]);

disp('【增强版特征值间隙分析结果】');
fprintf('一阶差分统计信息:\n');
fprintf('  均值: %.6f\n', mean(dy));
fprintf('  标准差: %.6f\n', std(dy));
fprintf('  阈值 (mean+std): %.6f\n', thresh);
fprintf('  超阈值点索引: %s\n', mat2str(jumps));
fprintf('\n👉 最终建议聚类数 k = %d\n', suggested_k);


