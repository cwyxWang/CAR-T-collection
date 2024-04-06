%% 该函数用来对原始明场图像进行预处理
function WF = WF_processing(WF)
%% 1、Sobel滤波提取边缘信息
sobel1 = [1 1 1;0 0 0;-1 -1 -1];
sobel2 = [1 0 -1;1 0 -1;1 0 -1];
edge1 = conv2(WF,sobel1,'same');   % same表示卷积后size不变
edge2 = conv2(WF,sobel2,'same');
WF = sqrt(edge1.^2 + edge2.^2);
WF(1:2,:) = 0;WF(end-1:end,:) = 0;WF(:,1:2) = 0;WF(:,end-1:end) = 0; % 将卷积的白边填0

%% 2、均值滤波平滑圆孔内信息
f = fspecial('average', [11,11]);
WF = imfilter(WF, f, "symmetric");

%% 3、自动阈值分割为二值图（百分比阈值或手动阈值）
% 将图像像素值排序
sorted_values = sort(WF(:),"descend");   % "descend"从高往低排序
% 根据百分位数计算阈值
percent_thr = 0.2;
threshold_index = round(percent_thr * numel(sorted_values));
threshold = sorted_values(threshold_index);
% 阈值分割
WF = uint8(WF > threshold);

end
