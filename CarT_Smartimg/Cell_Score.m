%% 该函数用来做细胞检测
function [GFP_BW,POS_Scores,GFP_nobw] = Cell_Score(TorB,GFP,Hole_POS,threshold,T1,T2)
POS_Scores = zeros(size(Hole_POS,1),3);   % 3列：[X_index,Y_index,score]
%% 去背景（用于评分）
% se = strel('square',50);
% GFP1 = imerode(GFP,se);   % 在50*50的box内取最小的像素值代替原图中心点的像素值
% GFP_sub = GFP - GFP1;
%% 将GFP进行预处理
% 高斯滤波
f = fspecial('gaussian',5,4);
GFP = imfilter(GFP, f, "symmetric");
% TopHat滤波
se = strel('disk',15);
GFP_nobw = imtophat(GFP,se);
% 二值化（后续最好改成Otsu！）
GFP_BW = uint16(GFP_nobw > threshold);
% threshold = graythresh(im2double(GFP));
% GFP_BW = imbinarize(GFP,threshold);
% 改善图像质量
se1 = strel('square',5);
se2 = strel('square',7);
GFP_BW = imclose(imopen(GFP_BW,se1),se2);
% GFP_BW = imfill(GFP_BW,'holes');
%% 对每个Cell进行细胞检测
for i = 1 : size(Hole_POS,1)
    % 先记录坐标
    POS_Scores(i,1) = Hole_POS(i,3);   % X_index
    POS_Scores(i,2) = Hole_POS(i,4);   % Y_index
    %% 在GFP图中crop出第i个孔的细胞小图Cell(128*128)
    h0 = Hole_POS(i,1);   w0 = Hole_POS(i,2);
    % 鲁棒性，crop到边缘为止
    h_min = max(h0-64,1);   h_max = min(h0+63,size(GFP,1));
    w_min = max(w0-64,1);   w_max = min(w0+63,size(GFP,2));
    % crop
    Cell = GFP_nobw(h_min:h_max,w_min:w_max);   % 在去过背景的图上Crop
    Cell_BW = GFP_BW(h_min:h_max,w_min:w_max);  

    %% 连通域处理
    % 参数设置
    min_Area = 200;
    % 1、先将面积过小的连通域（杂质）像素值清零
    Cell_regions = bwconncomp(Cell_BW);   % 检测连通域
    for j = 1 : Cell_regions.NumObjects
        Cell_Area = length(cell2mat(Cell_regions.PixelIdxList(j)));   % 得到第j个连通域的面积
        if Cell_Area < min_Area
            Cell_BW(cell2mat(Cell_regions.PixelIdxList(j))) = 0;   % 将杂质连通域在原图上清0
        end
    end
    % 2、计算连通域（腐蚀后）个数，若不为1~2个，直接PASS该细胞
%     se = strel('square',5);
%     Cell_regions = bwconncomp(imerode(Cell_BW,se));   % 重新检测连通域
%     Cell_num = Cell_regions.NumObjects;
%     if Cell_num < 1 || Cell_num > 2
%         continue;
%     end
    % 2、计算连通域个数，若不为1~3个，直接PASS该细胞
    Cell_regions = bwconncomp(Cell_BW);   % 重新检测连通域
    if Cell_regions.NumObjects < 1 || Cell_regions.NumObjects > 3   % 只要1~3个细胞
        continue;
    end
    % 3、计算剩下细胞的总面积，判断面积区间，若面积过大则直接PASS该细胞（小于100的已筛除）
    All_Cell_Area = 0;
    for j = 1 : Cell_regions.NumObjects
        Cell_Area = length(cell2mat(Cell_regions.PixelIdxList(j)));   % 得到第j个连通域的面积
        All_Cell_Area = All_Cell_Area + Cell_Area;   % 顺便算出1~2个细胞的总面积
    end
    if All_Cell_Area > T2
        continue;
    elseif All_Cell_Area > T1 || Cell_regions.NumObjects ~= 1
        if TorB   % 1:T细胞   0:靶细胞
            All_Cell_Area = All_Cell_Area * 1.25;   % 相当于分数*0.8
            score = sum(sum(Cell .* Cell_BW)) / All_Cell_Area;
        else
            score = 0.9;
        end
    else
        if TorB
            score = sum(sum(Cell .* Cell_BW)) / All_Cell_Area;
        else
            score = 1;
        end
    end

    % 评分：总像素强度÷总像素面积=平均强度(与面积无关)
    % 记录分数
    POS_Scores(i,3) = score;   % 记录score

end

end