%   rhoMat: 相关系数矩阵
%   Dist:   距离矩阵 (km)
%   dc:     第一变点距离（km）

dc = d_thr; % 第一变点距离
sigma = dc / 2;


W = Rho .* exp(-Dist.^2 / (2*sigma^2));
W(Dist > dc) = 0;  
W = max(W, W');    
W= abs(W);

k = 4; 

labels = spectralcluster(W, k, 'Distance', 'precomputed');

latlim = [15 55];
lonlim = [70 136];
figure('Color', 'w', 'Units','normalized','position',[0.15 0.15 0.5 0.7]);
m_proj('Equidistant Cylindrical','long',lonlim,'lat',latlim);
m_coast('patch',[0.96 0.96 0.96],'edgecolor','none');
hold on;
M_prov = m_shaperead('D:\Program Files\MATLAB\全国省行政区划\bou1_4l');
for i = 1:size(M_prov.ncst, 1)
    m_plot(M_prov.ncst{i,1}(:,1), M_prov.ncst{i,1}(:,2), 'color', [0.85 0.85 0.85], 'linewidth', 0.5);
end

k = max(labels); 
colors = lines(k); 
h_sta = zeros(k, 1); 

for i = 1:k
    idx = (labels == i);
    h_sta(i) = m_plot(L(idx), B(idx), 'o', ...
                   'MarkerFaceColor', colors(i,:), ...
                   'MarkerEdgeColor', 'k', ... 
                   'MarkerSize', 6, ...
                   'LineWidth', 0.5);
end


m_grid('linestyle',':','tickdir','out','linewidth',1,...
       'fontsize',12,'fontname','Times New Roman',...
       'xtick',80:10:130,'ytick',20:10:50);


cluster_names = arrayfun(@(x) sprintf('sub-region %d', x), 1:k, 'UniformOutput', false);

l_obj = legend(h_sta, cluster_names, ...
    'FontSize', 10, ...
    'FontName', 'Microsoft YaHei', ...
    'Location', 'southwest', ... 
    'Box', 'on', ...             
    'LineWidth', 0.8);

title(l_obj, 'Partition Details');
% title('中国 GNSS 站点谱聚类分区分布图', 'FontSize', 15, 'FontWeight', 'bold');

hold off;

