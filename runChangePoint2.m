
B=Pre.B;
L=Pre.L;
N=Pre.pcainterN;
U=Pre.pcainterU;
E=Pre.pcainterE;

assert(isvector(B) && isvector(L), 'B 与 L 必须是一维向量');
assert(size(U,2) == numel(B), 'U 的列数必须等于站点数');
nSta = numel(B);
fprintf('站点数量 = %d，时间序列长度 = %d\n', nSta, size(U,1));

%%=====================================================================
%%  2. 计算站点间的距离和 Pearson 相关系数
%%=====================================================================
deg2rad = @(x) x*pi/180;
phi    = deg2rad(B(:));
lambda = deg2rad(L(:));
R = 6371.0088;                     


Dist = zeros(nSta);
for i = 1:nSta-1
    dphi    = phi(i+1:end)   - phi(i);
    dlambda = lambda(i+1:end)- lambda(i);
    a = sin(dphi/2).^2 + cos(phi(i))*cos(phi(i+1:end)).*sin(dlambda/2).^2;
    c = 2*atan2(sqrt(a), sqrt(max(0,1-a)));
    Dist(i,i+1:end) = R*c;
end
Dist = Dist + Dist.';               % 对称化


U0 = detrend(E,'constant');
U0 = U0 - mean(U0,1);
Rho = corr(U0,'type','Pearson');


idxUpper   = find(triu(ones(nSta),1));   
distVec    = Dist(idxUpper);
rhoVec     = Rho(idxUpper);

%%=====================================================================
%%  3. 形成距离递增序列 y(d) 并平滑
%%=====================================================================
[distVecSorted, sortIdx] = sort(distVec);
rhoVecSorted = rhoVec(sortIdx);
winSize = 7;                         % 移动平均窗口（奇数）
y = movmean(rhoVecSorted, winSize, 'Endpoints','shrink');

%%=====================================================================
%%  4. 变点检测（findchangepts / ischange）
%%=====================================================================
maxChanges = 3;                
statistic  = 'mean';

if exist('findchangepts','file')
%     if verLessThan('matlab','9.9')
        bkpts = findchangepts(y,'Statistic',statistic,'MaxNumChanges',maxChanges);
%     else
%         [bkpts,~] = findchangepts(y,'Statistic',statistic,...
%                                     'MaxNumChanges',maxChanges,'Penalty','BIC');
%     end
else
    changeLog = ischange(y,'mean','MaxNumChanges',maxChanges);
    idxRaw    = find(changeLog);
    diffIdx   = [true; diff(idxRaw)>1];
    bkpts     = idxRaw(diffIdx);
end

if isempty(bkpts)
    error('未检测到任何变点，请检查数据或放宽 maxChanges。');
end

fprintf('\n检测到的变点（按距离递增排序）:\n');
for k = 1:numel(bkpts)
    idx = bkpts(k);
    fprintf('  变点 %d : 索引 %d, 距离 = %.2f km, 相关系数 = %.3f\n',k,idx,distVecSorted(idx),y(idx));
end


firstIdx = bkpts(1);% 采用第一个变点作为阈值
d_thr    = distVecSorted(firstIdx);
rho_thr  = y(firstIdx);
fprintf('\n=== 采用的阈值 (第一个变点) ===\n');
fprintf('  距离阈值    d_thr = %.2f km\n', d_thr);
fprintf('  相关系数阈值 rho_thr = %.3f\n', rho_thr);

%%=====================================================================
%%  5. 置信区间估计 —— 对原始时间序列 U 进行 BlockBootstrap
%%=====================================================================
nBoot      = 500;                % 抽样次数
confLevel  = 0.95;
alpha      = (1-confLevel)/2;

blockLen   = 30;                 % 时间块长度
nTime      = size(U,1);
nBlocks    = ceil(nTime/blockLen);


boot_d_thr   = NaN(nBoot,1);
boot_rho_thr = NaN(nBoot,1);

fprintf('\n=== 开始基于时间序列的 BlockBootstrap（%d 次）===\n', nBoot);

parfor b = 1:nBoot

    U_boot = zeros(nTime,nSta);
    t = 1;
    while t <= nTime
        blkIdx   = randi(nBlocks);                    
        startT   = (blkIdx-1)*blockLen + 1;
        endT     = min(startT+blockLen-1, nTime);
        blockData = U(startT:endT,:);                  
        lenBlock  = size(blockData,1);
        if t+lenBlock-1 <= nTime
            U_boot(t:t+lenBlock-1,:) = blockData;
            t = t + lenBlock;
        else
            need = nTime - t + 1;
            U_boot(t:end,:) = blockData(1:need,:);
            t = nTime + 1;
        end
    end

    U0_boot = detrend(U_boot,'constant');
    U0_boot = U0_boot - mean(U0_boot,1);
    Rho_boot = corr(U0_boot,'type','Pearson');
    rhoVec_boot = Rho_boot(idxUpper);   

 
    rhoVecSorted_boot = rhoVec_boot(sortIdx);
    y_boot = movmean(rhoVecSorted_boot, winSize, 'Endpoints','shrink');

 
    if exist('findchangepts','file')
%         if verLessThan('matlab','9.9')
            bk = findchangepts(y_boot,'Statistic','mean','MaxNumChanges',maxChanges);
%         else
%             [bk,~] = findchangepts(y_boot,'Statistic','mean',...
%                                      'MaxNumChanges',maxChanges,'Penalty','BIC');
%         end
    else
        changeLog = ischange(y_boot,'mean','MaxNumChanges',maxChanges);
        idxRaw    = find(changeLog);
        diffIdx   = [true; diff(idxRaw)>1];
        bk        = idxRaw(diffIdx);
    end


    if ~isempty(bk)
        first = bk(1);

        if first > 5 && first < numel(y_boot)-5
            boot_d_thr(b)   = distVecSorted(first);
            boot_rho_thr(b) = y_boot(first);
        end
    end
end


boot_d_thr   = boot_d_thr(~isnan(boot_d_thr));
boot_rho_thr = boot_rho_thr(~isnan(boot_rho_thr));


ci_d   = quantile(boot_d_thr, [alpha, 1-alpha]);   
ci_rho = quantile(boot_rho_thr, [alpha, 1-alpha]);

fprintf('\nBlockBootstrap 完成（有效抽样 %d 次）\n', numel(boot_d_thr));
fprintf('距离阈值 95%% CI : [%.2f , %.2f] km\n', ci_d(1), ci_d(2));
fprintf('相关系数阈值 95%% CI : [%.3f , %.3f]\n', ci_rho(1), ci_rho(2));


%%=====================================================================
%%  6. 画图
%%=====================================================================

colScatter = [0.40 0.40 0.40];       
colSmooth  = [0 0.4470 0.7410];      
colThresh  = [0.8500 0.3250 0.0980]; 

figure('Name','Distance–Correlation with Change?Points & CI', ...
       'NumberTitle','off','Units','centimeters','Position',[5 5 16 10]);
hScatter = scatter(distVecSorted, rhoVecSorted, 12, colScatter, ...
                  'filled','Marker','o','MarkerEdgeColor','none'); hold on;
hSmooth = plot(distVecSorted, y, 'Color',colSmooth, ...
               'LineWidth',1.5,'LineStyle','-');

hChange = gobjects(numel(bkpts),1);
for k = 1:numel(bkpts)
    hChange(k) = xline(distVecSorted(bkpts(k)), '--', ...
                       'Color',colThresh,'LineWidth',1.2);
end


hThresh = xline(d_thr, '-', 'Color',colThresh, 'LineWidth',2);


dx = -0.8 * (max(distVecSorted)-min(distVecSorted));  
dy_up   = 0.01 * (max(y)-min(y));    
dy_down = 0.2 * (max(y)-min(y));   


txtDist = sprintf('d_{thr}=%.1f km\n95%% CI: [%.1f , %.1f] km', ...
                  d_thr, ci_d(1), ci_d(2));
text(d_thr-dx, max(y)-dy_up, txtDist, ...
     'FontName','Arial','FontSize',9,'FontWeight','bold', ...
     'Color',colThresh,'HorizontalAlignment','right','Interpreter','none');


txtRho = sprintf('\\rho_{thr}=%.3f\n95%% CI: [%.3f , %.3f]', ...
                 rho_thr, ci_rho(1), ci_rho(2));
text(d_thr-dx, min(y)+dy_down, txtRho, ...
     'FontName','Arial','FontSize',9,'FontWeight','bold', ...
     'Color',colThresh,'HorizontalAlignment','right','Interpreter','none');


dxTxt = 0.02*(max(distVecSorted)-min(distVecSorted));
dyTxt = 0.03*(max(rhoVecSorted)-min(rhoVecSorted));
for k = 1:numel(bkpts)
    idx = bkpts(k);
    dV  = distVecSorted(idx);
    rV  = y(idx);
    txt = sprintf('d=%.1f km\n\\rho=%.3f',dV,rV);
    text(dV+dxTxt, rV+dyTxt, txt, ...
         'FontName','Arial','FontSize',9,'FontWeight','bold', ...
         'Color',colThresh,'HorizontalAlignment','left', ...
         'VerticalAlignment','bottom');
end


xlabel('Station pair distance (km)',...
       'FontName','Arial','FontSize',10,'Interpreter','none');
ylabel('Pearson correlation coefficient',...
       'FontName','Arial','FontSize',10,'Interpreter','none');
title('Distance--Correlation with detected change-points and 95% CI',...
      'FontName','Arial','FontSize',12,'Interpreter','none');

legend([hScatter, hSmooth, hChange(1), hThresh], ...
       {'Raw scatter','Smoothed curve','Detected change-points','Threshold'}, ...
       'Location','best','FontName','Arial','FontSize',9,...
       'Interpreter','none','Box','off');


set(gca,'Box','off','LineWidth',1,'FontName','Arial','FontSize',9);
grid on;
xlim([0 4500]);    
ylim([-0.5 1]);    
hold off;

