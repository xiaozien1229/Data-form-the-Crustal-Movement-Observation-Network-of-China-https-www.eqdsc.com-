% 假设已知：
%   rhoMat: nSta x nSta 相关系数矩阵
%   Dist:   nSta x nSta 距离矩阵 (km)
%   dc:     第一变点距离（km）

dc = d_thr; % 例如 300 km
sigma = dc / 2;

% 构造相似性核
W = Rho .* exp(-Dist.^2 / (2*sigma^2));
W(Dist > dc) = 0;  % 超过变点距离置零
W = max(W, W');    % 对称化
W= abs(W);
% 谱聚类（需要指定簇数k，或用特征间隙估计）
k = 8; % 可根据中国大地构造块体初设（如青藏、华南、华北等）

% 使用 MATLAB Statistics and Machine Learning Toolbox
labels = spectralcluster(W, k, 'Distance', 'precomputed');

% 可视化
% scatterm(B, L, [], labels, 'filled');
% title('谱聚类分区结果');

%% 3. 绘图阶段 (优化版：分簇着色)
latlim = [15 55];
lonlim = [70 136];
figure('Color', 'w', 'Units','normalized','position',[0.15 0.15 0.5 0.7]);
m_proj('Equidistant Cylindrical','long',lonlim,'lat',latlim);

% A. 绘制陆地背景 (浅灰色填充)
m_coast('patch',[0.96 0.96 0.96],'edgecolor','none');
hold on;

% B. 绘制省界 (极细的浅灰色)
M_prov = m_shaperead('D:\Program Files\MATLAB\全国省行政区划\bou1_4l');
for i = 1:size(M_prov.ncst, 1)
    m_plot(M_prov.ncst{i,1}(:,1), M_prov.ncst{i,1}(:,2), 'color', [0.85 0.85 0.85], 'linewidth', 0.5);
end

% D. 核心修改：根据聚类标签 labels 绘制不同颜色的站点
% 定义颜色映射 (使用适合出版的离散颜色)
k = max(labels); 
colors = lines(k); % 或者使用 [1 0 0; 0 1 0; 0 0 1; 1 0 1] 等自定义 RGB
h_sta = zeros(k, 1); % 初始化句柄用于图例

for i = 1:k
    % 找到属于第 i 簇的索引
    idx = (labels == i);
    % 绘制该簇的站点
    h_sta(i) = m_plot(L(idx), B(idx), 'o', ...
                   'MarkerFaceColor', colors(i,:), ...
                   'MarkerEdgeColor', 'k', ... % 黑色描边让点更清晰
                   'MarkerSize', 6, ...
                   'LineWidth', 0.5);
end

% E. 设置坐标轴网格
m_grid('linestyle',':','tickdir','out','linewidth',1,...
       'fontsize',12,'fontname','Times New Roman',...
       'xtick',80:10:130,'ytick',20:10:50);

%% 4. 完善图例与标题
% 生成图例标签
cluster_names = arrayfun(@(x) sprintf('sub-region %d', x), 1:k, 'UniformOutput', false);
% 2. 执行绘制图例命令 (关键缺失步)
% h_sta 是你在循环中存储的句柄数组
l_obj = legend(h_sta, cluster_names, ...
    'FontSize', 10, ...
    'FontName', 'Microsoft YaHei', ...
    'Location', 'southwest', ... % 图例放在左下角
    'Box', 'on', ...             % 显示图例框
    'LineWidth', 0.8);

% 3. 给图例加个标题 (可选，增加专业感)
title(l_obj, 'Partition Details');

% % 4. 添加图像主标题
% title('中国 GNSS 站点谱聚类分区分布图', 'FontSize', 15, 'FontWeight', 'bold');

hold off;

