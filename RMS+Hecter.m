clearvars -except  Pre labels fullpath proj_data path
proj_data.Llon=70;
proj_data.Rlon=136;
proj_data.Llat=15;
proj_data.Ulat=55;

[name,path]=uigetfile();
fullpath=fullfile(path,name);
%% %%%分区pca+计算振幅相位%%%%%%%%%%%%%%%%
labels=region_labels;%板块分区的站点区号
rms_out0.rms_FQ_afpc = cell(1, 8);
for k=1:8
    ooo_yyy_X999{k}=N(:,labels==k);%改方向
    idx = (labels == k);
    % 提取对应的 B 和 L
    lon_reordered = B(idx);
    lat_reordered = L(idx);
    [ooo_yyy_X2,ooo_yyy_D2,ooo_yyy_SRs_PC2]=ana_PCA(ooo_yyy_X999{k},1);
    [ooo_yyy_ISD3,ooo_yyy_SUMD3]=CCR(ooo_yyy_D2);%C Rate
    cme1=ooo_yyy_X2'*ooo_yyy_SRs_PC2;
    lon_all{k}=lon_reordered;
    lat_all{k}=lat_reordered;
    cme1_all{k} = cme1;
    ooo_yyy_final1{k}=ooo_yyy_X999{k}(:,:)-cme1;%after PCA on all sub_station
    % [output1,output2,output3]=ana_cme(ooo_yyy_X999,3,'ICA');
% [X_ZF1]=Amplitude_phase(lat_reordered',lon_reordered',cme1',fullpath,proj_data);
    %% %%%rms%%%%%%%%%%%%%%%%
    rms_out0.rms_FQ_afpc{k}=rms(ooo_yyy_final1{k});
    rms_out0.rms_mean_FQ_afpc(k)=mean(rms(ooo_yyy_final1{k}));
        
end
ooo_yyy_final1_matrix = horzcat(ooo_yyy_final1{:});
L_matrix=vertcat(lat_all{:});
B_matrix=vertcat(lon_all{:});
cme_matrix = horzcat(cme1_all{:});
rms1=horzcat(rms_out0.rms_FQ_afpc{:});
mean(rms_out0.rms_mean_FQ_afpc)

%% %%%整体pca+计算振幅相位%%%%%%%%%%%%%%%%
[ooo_yyy_X,ooo_yyy_D,ooo_yyy_SRs_PC]=ana_PCA(Pre.pcainterN,1);%改方向
[ooo_yyy_ISD,ooo_yyy_SUMD]=CCR(ooo_yyy_D);%C Rate
ppp3=ooo_yyy_X'*ooo_yyy_SRs_PC;
ooo_yyy_final=Pre.pcainterN-ppp3;%after PCA on all sub_station%改方向

rms_out0.rms_ooo_yyy_afpc=rms(ooo_yyy_final);
rms_out0.rms_ooo_yyy_mean_afpc=mean(rms(ooo_yyy_final));

rms_out0.rms_ooo_yyy_bf=rms(Pre.pcainterN);%改方向
rms_out0.rms_ooo_yyy_mean_bf=mean(rms(Pre.pcainterN));%改方向

[X_ZF1]=Amplitude_phase(L',B',ppp3',fullpath,proj_data);
%% %%%对比分析振幅相位与计算D值%%%%%%%%%%%%%%%%
% 调整L、B、CME

[~, index2] = ismember(L, L_matrix);% 返回的 index2 表示：L2(index2) 就能得到和 L1 一样的顺序
L2 = L_matrix(index2);% 2. 同步调整 L2（此时 L2_new 应该等于 L1）
% [~, index2] = ismember(B, B_matrix);
% B2 = B_matrix(index2);
cme_matrix_new = cme_matrix(:, index2);% 3. 同步调整 cme_matrix 的列
rms2=rms1(:, index2);
ooo_yyy_final1_matrix2=ooo_yyy_final1_matrix(:, index2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=NTOL_N+NTAL_N+HYDL_N;
x2=NTOL_U2+NTAL_U2+HYDL_U2;
x3=cme_matrix_new-x1;
x5=rms_out0.rms_FQ_afpc{1};
x4=ppp3-x1;
[X_ZF1]=Amplitude_phase(L',B',x3',fullpath,proj_data);
[X_ZF2]=Amplitude_phase(L',B',x1',fullpath,proj_data);
[X_ZF1]=Amplitude_phase(L',B',cme_matrix_new',fullpath,proj_data);
[X_ZF1]=Amplitude_phase(L_matrix',B_matrix',cme_matrix',fullpath,proj_data);
[X_ZF2]=Amplitude_phase(Lon,lat,x4',fullpath,proj_data);
% [X_ZF2]=Amplitude_phase(point_new2_lon',point_new2_lat',U_ppp2_matrix',fullpath,proj_data);%CME_ica1
% x2_1 = x1(:, unique_sorted_ids);%把负荷按照分区的数字排序

X1=X_ZF1(1,:);
X2=X_ZF1(2,:);
X3=X_ZF1(3,:);   
X4=X_ZF1(4,:);
X5=X_ZF2(1,:);
X6=X_ZF2(2,:);

amp1 = sqrt(X1.^2 + X2.^2);     % 振幅
pha1 = atan2(X1, X2);           % 相位（单位：弧度）
% 如需以度为单位输出：
pha_deg1 = mod(rad2deg(pha1), 360)';
amp2 = sqrt(X3.^2 + X4.^2);     % 振幅
pha2 = atan2(X3, X4);           % 相位（单位：弧度）
% 如需以度为单位输出：
pha_deg2 = mod(rad2deg(pha2), 360);

mean_amp=mean(amp1);
mean_deg=mean(pha_deg1);
D = value_D(X1,X5,X2,X6);
D = value_D(X1,X2,X5,X6);
%% %%%画RMS与噪声图%%%%%%%%%%%%%%%%
rms_final_subarea = zeros(size(labels)); 
lon_fq = zeros(213, 1);
lat_fq = zeros(213, 1);
for k = 1:8
    idx = (labels == k);
    rms_final_subarea(idx) = rms_out0.rms_FQ_afpc{k}(:);
    lon_fq(idx) = L(idx); 
    lat_fq(idx) = B(idx);
end

data_values = 1 - (rms2./ rms_out0.rms_ooo_yyy_bf);
data_values = 1 - (rms_out0.rms_ooo_yyy_afpc ./ rms_out0.rms_ooo_yyy_bf);
data_values =dpl{1};

figure('Units','normalized','position',[0.20 0.2 0.45 0.5])
m_proj('Equidistant Cylindrical','long',[70 136],'lat',[15 55]);
% m_contourf(ooo_L_yyy,ooo_B_yyy,ooo_yyy_SRs_PC,70,'linestyle','none');
m_grid('linestyle','none','tickdir','out','linewidth',1.6,'fontsize',18,'fontname','Times New Roman','fontweight','bold','xtick',[80:10:136],'ytick',[20:10:55]);
m_coast('patch',[0.86275,0.86275,0.86275],'edgecolor','none')
maa=shaperead(fullpath);
m_line([maa(:).X],[maa(:).Y],'color','k');
% colormap(winter(50));
colormap("parula");
hold on
% m_scatter(L, B, 100,(rms_out0.rms_final_subarea./rms_out0.rms_ooo_yyy_bf), 'filled');
% % m_scatter(L, B, 100,(rms_out0.rms_ooo_yyy_afpc./rms_out0.rms_ooo_yyy_bf), 'filled');
% 
% caxis([0 1]);%1：adp 2:pca 3:raw
m_scatter(L, B, 100, data_values, 'filled');
% m_scatter(L_matrix, B_matrix, 100, data_values, 'filled');
% m_scatter(filterL, filterB, 100, data_values, 'filled');
% m_scatter(ooo_L_yyy, ooo_B_yyy, 50, dwh{1}, 'filled');
% m_scatter(ooo_L_yyy, ooo_B_yyy, 50, dwh{1}, 'filled');
% m_scatter(L, B,100,'k','^','filled');%[0.8500 0.3250 0.0980]
caxis([0 ceil(max(data_values))]);%1：adp 2:pca 3:raw

hold on
%%%%%

