clc;clear
t1 = tic;
%% 参数
bathpath = 'D:\20231226_2\WF';
step = 1350;   % um
% % 分割阈值
% thr_GFP = 50;
% thr_RFP = 7;
% % 体积阈值
% T1_GFP = 470;
% T2_GFP = 800;
% T1_RFP = 510;
% T2_RFP = 1000;

thr_GFP = 50;
thr_RFP = 50;
% 体积阈值
T1_GFP = 470;
T2_GFP = 800;
T1_RFP = 470;
T2_RFP = 800;

%% 其余参数
% 三张图像路径及保存路径
temp = dir(fullfile(bathpath,'K_*'));
floder_name = {temp.name};

WF_folder = fullfile(bathpath,floder_name{1});
temp = dir(fullfile(WF_folder,'*.tif'));
file_name = {temp.name};
WFpath = fullfile(WF_folder,file_name{1});

GFP_folder = fullfile(bathpath,floder_name{2});
temp = dir(fullfile(GFP_folder,'*.tif'));
file_name = {temp.name};
GFPpath = fullfile(GFP_folder,file_name{1});

RFP_folder = fullfile(bathpath,floder_name{3});
temp = dir(fullfile(RFP_folder,'*.tif'));
file_name = {temp.name};
RFPpath = fullfile(RFP_folder,file_name{1});

% WFpath = fullfile(bathpath,'K_1','K_1_MMStack.ome.tif');
% GFPpath = fullfile(bathpath,'K_2','K_2_MMStack.ome.tif');
% RFPpath = fullfile(bathpath,'K_3','K_3_MMStack.ome.tif');
full_WF_savepath = fullfile(bathpath,'WF_Merge.tif');
full_GRFP_savepath = fullfile(bathpath,'GFP_RFP_Merge.tif');
full_GFP_savepath = fullfile(bathpath,'GFP.tif');
full_TopHatR_savepath = fullfile(bathpath,'TopHat_RFP.tif');
full_TopHatG_savepath = fullfile(bathpath,'TopHat_GFP.tif');
% 表格（共3页：所有孔坐标、高分细胞、所有细胞的总,GFP,RFP分数）
% Excel_savepath = fullfile(bathpath, 'WF(Y_mm,X_mm,X_index,Y_index,N1,N2,N3)-Cell(X_index,Y_index,score,GFP,RFP)-Allcell.xlsx');
Excel_savepath = fullfile(bathpath, 'result.xlsx');
% 表格每页矩阵
WF_POS = [];   % 1、用来存每个孔的坐标信息：(Y_pos,X_pos,X_index,Y_index,N1,N2,N3)   左上角为(1,1)
Cell_POS_Scores = [];   % 23、用来存每个细胞的坐标信息：(X_index,Y_index,score(,GFP,RFP))
for roi = 1:9   % 默认每个stack都是9个ROI
    t2 = tic;
    disp(['----------正在检测第',num2str(roi),'个ROI----------']);
    %% 读取图像
    WF = imread(WFpath,roi);
    GFP = imread(GFPpath,roi);
    RFP = imread(RFPpath,roi);
    % 最开始时定义，用来存放9个ROI结果。
    if roi == 1
        WFs = zeros([size(WF),9]);     WFs_BW = zeros([size(WF),9]);   
        GFPs = zeros([size(GFP),9]);   GFPs_BW = zeros([size(GFP),9]);   
        RFPs = zeros([size(RFP),9]);   RFPs_BW = zeros([size(RFP),9]);     
        Exact_Circles = zeros([size(WF),9]);   % 准确定位的⚪
        Template_Circles = zeros([size(WF),9]);   % 大模板⚪
        % TopHat之后的图像
        GFP_nobws = zeros([size(WF),9]);
        RFP_nobws = zeros([size(WF),9]);
    end

    %% 图像预处理，提取明场图像中孔信息，得到二值图（耗时0.15s）
    WF_BW = WF_processing(WF);

    %% 金字塔结构，遍历模板位置进行匹配（耗时5+1+2+1.5+2=11.5s）（0.65+0.70+1.20+0.80+1.15=4.5s）
    % 输入参数为：(WF,scale,x_best,x_range,y_best,y_range,theta_best,theta_range,theta_step)
    % 16倍下采样，默认遍历范围XY为-20~20，theta为-5~5°
    [~,~,x_best,y_best,theta_best] = WF_225Detect_1(WF_BW,16, 0,20, 0,20, 2, 0,5,1);    % 滑窗范围可以考虑
    % 8倍下采样
    [~,~,x_best,y_best,theta_best] = WF_225Detect_1(WF_BW,8 ,x_best*2,10,y_best*2,10, 2, theta_best,3,0.3); 
    % 4倍下采样
    [~,~,x_best,y_best,theta_best] = WF_225Detect_1(WF_BW,4 ,x_best*2,4,y_best*2,4, 1, theta_best,0.6,0.2);
    % 2倍下采样
    [~,~,x_best,y_best,theta_best] = WF_225Detect_1(WF_BW,2 ,x_best*2,2,y_best*2,2, 1, theta_best,0.3,0.15);
    % 1倍下采样（range=2实际对应5*5的滑窗）
    [template_center,template_circle,x_best,y_best,theta_best] = WF_225Detect_1(WF_BW,1 ,x_best*2,2,y_best*2,2, 1, theta_best,0.1,0.1);
    
    %% 利用template_circle中225个孔心坐标来crop孔的小图，并在每个crop小图中遍历得到225个准确孔心坐标（耗时3.3s）
    Hole_POS = WF_1Detect_1(WF_BW,template_center);   % 这里出来是整数

    %% 最后在每个准确的孔心坐标生成⚪，以此来验证匹配结果是否准确（以下全部耗时0.4s）
    hole_circle = zeros(size(WF));
    for i = 1 : size(Hole_POS,1)   % 找出所有孔心（理论是225个）
        hole_circle(Hole_POS(i,1),Hole_POS(i,2)) = 1;
    end
    hole_circle = bwdist(hole_circle) < 31/0.65;   % 生成⚪

    %% 将XY坐标换算成9个ROI大图中的坐标(pixel)，并保存到WF_POS中
    Hole_POS(:,5) = roi;   % N1
    % 根据roi算出后两列的XY_index
    Hole_POS(:,3) = Hole_POS(:,3) + mod((roi-1),3) * 15;      % X_index(1~15 → 1~45)
    Hole_POS(:,4) = Hole_POS(:,4) + floor((roi-1)/3) * 15;    % Y_index(1~15 → 1~45)

    %% 细胞评分
    [GFP_BW,GFP_POS_Scores_roi,GFP_nobw] = Cell_Score(1,GFP,Hole_POS, thr_GFP, T1_GFP, T2_GFP); 
    [RFP_BW,RFP_POS_Scores_roi,RFP_nobw] = Cell_Score(0,RFP,Hole_POS, thr_RFP, T1_RFP, T2_RFP); 
    % 保存TopHat之后的图像
    GFP_nobws(:,:,roi) = GFP_nobw;
    RFP_nobws(:,:,roi) = RFP_nobw;
    % 算一个加权总分（靶细胞有就行(0or1)，T细胞按分数来）
    weighted_score = GFP_POS_Scores_roi(:,3) .* RFP_POS_Scores_roi(:,3);
    Cell_POS_Scores_roi = [GFP_POS_Scores_roi(:,1:2),weighted_score,GFP_POS_Scores_roi(:,3),RFP_POS_Scores_roi(:,3)];

    %% 继续将XY坐标换算成9个ROI大图中的坐标(pixel)，并保存到WF_POS中
    % 根据roi算出在9张大图中的总坐标（到这里还是像素数）
    Hole_POS(:,2) = round(Hole_POS(:,2) + mod((roi-1),3) * (step/0.65));    % W(X)
    Hole_POS(:,1) = round(Hole_POS(:,1) + floor((roi-1)/3) * (step/0.65));  % H(Y)
    % 上下累加矩阵
    WF_POS = [WF_POS;Hole_POS];
    Cell_POS_Scores = [Cell_POS_Scores;Cell_POS_Scores_roi];

    %% 存到9个ROI的总和图像里
    % WF_Merge.tif
    WFs(:,:,roi) = uint16(WF);     % 输出预处理某一个阶段的图
    WFs_BW(:,:,roi) = uint16(WF_BW);
    Exact_Circles(:,:,roi) = uint16(hole_circle);
    Template_Circles(:,:,roi) = uint16(template_circle);
    % GFP1-9_Merge.tif
    GFPs(:,:,roi) = uint16(GFP);
    GFPs_BW(:,:,roi) = uint16(GFP_BW);
    % RFP1-9_Merge.tif
    RFPs(:,:,roi) = uint16(RFP);
    RFPs_BW(:,:,roi) = uint16(RFP_BW);

    toc(t2);
end

%% 拼接完整的芯片图像   位移台步长对应的像素数
% WF拼接
full_WF = Merge_1to9(WFs,step);   % 完整WF原图
full_WF_BW = Merge_1to9(WFs_BW,step);   % 完整WF_BW图
full_template_circle = Merge_1to9(Template_Circles,step);   % 完整大模板圆孔图
full_hole_circle = Merge_1to9(Exact_Circles,step);   % 完整定位圆孔图
% GFP拼接
full_GFP = Merge_1to9(GFPs,step);   % 完整GFP原图
full_GFP_BW = Merge_1to9(GFPs_BW,step);   % 完整GFP_BW图
% RFP拼接
full_RFP = Merge_1to9(RFPs,step);   % 完整GFP原图
full_RFP_BW = Merge_1to9(RFPs_BW,step);   % 完整GFP_BW图
% TopHat之后的图像拼接
full_GFP_nobws = Merge_1to9(GFP_nobws,step); 
full_RFP_nobws = Merge_1to9(RFP_nobws,step); 

%% 筛选并在GFP原图中框选出分数高的细胞
% 筛选前Cell_num个细胞
Cell_POS_Scores = sortrows(Cell_POS_Scores, 3, "descend");   % 全部细胞的评分，按分数从高到底排序
% % 随机打乱（若关闭功能则将thr_score改成无穷大）
% thr_score = 400;   % 高分阈值
% good = find(Cell_POS_Scores(:,3) > thr_score);   % 高于400分的随机打乱顺序
% if length(good) > Cell_num
%     Cell_good = Cell_POS_Scores(randperm(length(good)),:);   % 将超过400分的随机打乱
%     Cell_good = Cell_good(1:Cell_num,:);
% else
Cell_num = 0;
for i = 1:size(Cell_POS_Scores,1)
    if Cell_POS_Scores(i,3)>0
        Cell_num = Cell_num+1;
    else
        break
    end
end
Cell_good = Cell_POS_Scores(1:Cell_num,:);   % 筛选出分最高的一部分
% end
%% 标记在黑图中
WF_POS1 = WF_POS;   % 没啥用，debug时重复调用时方便，因为下面会把这个矩阵前两列换算成小数mm，再调用就会报错。调试好后把这行删掉下面1去掉
Cell_box = zeros(size(full_GFP));
for i = 1:Cell_num
    % 找出第i个细胞在WF_POS中对应行的Y_POS和X_POS
    H = WF_POS1(WF_POS1(:,3) == Cell_good(i,1) & WF_POS1(:,4) == Cell_good(i,2), 1);
    W = WF_POS1(WF_POS1(:,3) == Cell_good(i,1) & WF_POS1(:,4) == Cell_good(i,2), 2);
    % 在原图上标记出高分细胞
    H_min = max(H-64,1);     H_max = min(H+63,size(Cell_box,1));
    W_min = max(W-64,1);     W_max = min(W+63,size(Cell_box,2));
    % 画方框，宽度为7
    Cell_box(H_min:H_max,W_min:(W_min+6)) = 65535;
    Cell_box(H_min:H_max,(W_max-6):W_max) = 65535;
    Cell_box(H_min:(H_min+6),W_min:W_max) = 65535;
    Cell_box((H_max-6):H_max,W_min:W_max) = 65535;
end

%% 保存图像
% 保存为三通道图像WF_Merge.tif
% full_WF_merge = cat(3, uint16(full_hole_circle), uint16(full_WF_BW), uint16(full_template_circle), uint16(full_WF));
% imwrite(full_WF_merge(:,:,1)*65535,full_WF_savepath);   % 1、最终准确圆孔模板图（*65535是为了方便查看）
% imwrite(full_WF_merge(:,:,2)*65535,full_WF_savepath,'WriteMode','append');   % 2、WF_BW
% imwrite(full_WF_merge(:,:,3)*65535,full_WF_savepath,'WriteMode','append');   % 3、粗略圆孔模板图
% imwrite(full_WF_merge(:,:,4),full_WF_savepath,'WriteMode','append');   % 4、WF原图

full_WF_merge = cat(3, uint8(full_hole_circle), uint8(full_WF_BW),uint8(zeros(size(full_WF_BW))));
imwrite(full_WF_merge*255,full_WF_savepath);   % 1、最终准确圆孔模板图（*65535是为了方便查看）


% 保存为五通道图像GRFP_Merge.tif
% full_GRFP_merge = cat(3, uint16(Cell_box), uint16(full_GFP), uint16(full_GFP_BW),uint16(full_RFP),uint16(full_RFP_BW));
% imwrite(full_GRFP_merge(:,:,1),full_GRFP_savepath);   % 1、细胞检测筛选图（BoundingBox）
% imwrite(full_GRFP_merge(:,:,2),full_GRFP_savepath,'WriteMode','append');         % 2、GFP
% imwrite(full_GRFP_merge(:,:,3)*65535,full_GRFP_savepath,'WriteMode','append');   % 3、GFP_BW
% imwrite(full_GRFP_merge(:,:,4),full_GRFP_savepath,'WriteMode','append');         % 4、RFP
% imwrite(full_GRFP_merge(:,:,5)*65535,full_GRFP_savepath,'WriteMode','append');   % 5、RFP_BW

full_GRFP_merge = cat(3,uint8(full_RFP_BW), uint8(full_GFP_BW), uint8(Cell_box));
imwrite(full_GRFP_merge*255,full_GRFP_savepath);   % 1、细胞检测筛选图（BoundingBox）
imwrite(uint16(full_GFP),full_GFP_savepath);
% 保存TopHat之后的图
full_TopHat_merge = cat(3, uint16(full_GFP_nobws), uint16(full_RFP_nobws));
imwrite(full_TopHat_merge(:,:,1),full_TopHatG_savepath);  
imwrite(full_TopHat_merge(:,:,2),full_TopHatR_savepath);

%% 排序后保存xlsx表格
% 第一页（孔坐标）：按照N1(1~9)，N2(1~9)，N3(1~25)排序
WF_POS = sortrows(WF_POS, 7);   % 默认从低到高
WF_POS = sortrows(WF_POS, 6);
WF_POS = sortrows(WF_POS, 5);
% % 第二页（高分细胞）：按照X_index(1~45)，Y_index(1~45)排序
% Cell_good = sortrows(Cell_good, 2);
% Cell_good = sortrows(Cell_good, 1);
% 将前两列坐标换算成mm，注意！位移台的原点在左下角，因此y坐标要乘负数！
WF_POS(:,1) = WF_POS(:,1) * (-0.00065);   % pixel换算成mm（0.65um/pixel）
WF_POS(:,2) = WF_POS(:,2) * 0.00065;   % pixel换算成mm（0.65um/pixel）
writematrix(WF_POS,Excel_savepath,'Sheet',1,'Range','A1');   % Sheet1：所有孔坐标
writematrix(Cell_good,Excel_savepath,'Sheet',2,'Range','A1');   % Sheet2：高分细胞
writematrix(Cell_POS_Scores,Excel_savepath,'Sheet',3,'Range','A1');   % Sheet3：所有细胞评分

%% 如果检测出的细胞数量不够，则警告
% if(Cell_good(Cell_num,3) == 0)
%     warning(['检测出的有效细胞不足',num2str(Cell_num),'个！']);
% end
fprintf(['检测出细胞',num2str(Cell_num),'组\n'])

disp('运行结束，该程序总耗时：')
toc(t1);