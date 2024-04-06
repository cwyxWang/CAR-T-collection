%% 该函数用来在当前x0,y0,theta位置创建模板，并与预处理后的原始明场图像进行模板匹配评分。
%% 返回参数：score：分数，template_center：当前模板中每个孔的中心坐标。
% scale:下采样比例   
% x0,y0:左上角孔的左上角相切坐标（像素数，索引从0开始，+1后才是坐标，最好整数）   theta:顺时针倾斜角度（deg°，支持小数）
function [score,template_center,template_circle] = WF_225Detect_2(WF,scale,x0,y0,theta)   
%% 芯片参数
pixel2um = 0.65 * scale;   % 10X明场下，0.65um=1格像素
R = 25 / pixel2um;   % 孔的半径25um，换算成像素数
d = 30 / pixel2um;   % 最小间隔30um，中等间隔60um，最大间隔120um，换算成像素数

%% 创建相同size的孔中心标准坐标模板（左上角孔的左上角相切处为原点）
template_center = zeros(size(WF));
% 得到标准的15个孔心坐标(像素数)
XY_std = zeros(15,1);
for i = 1:15 
    XY_std(i) = R + (i-1) * (2 * R + d) + floor((i-1) / 5) * d;   % [25,……,1205um] [38,1854pixel]
end
POS_std = ones(3,225);
for i = 1:225
    POS_std(1,i) = XY_std(ceil(i/15));   % H(Y)
    POS_std(2,i) = XY_std(mod((i-1),15) + 1);   % W(X)
end

%% 将标准坐标变换到x0,y0,theta(逆时针)
% 放射变换矩阵
theta = deg2rad(theta);   % 输入为角度，转换为弧度
A = [cos(theta), -sin(theta), x0;
     sin(theta),  cos(theta), y0;
        0,            0,       1];
hole_POS = A * POS_std;

% 生成·模板，理论上225个点按顺序像素强度为1~225，旋转和平移过程中可能会有超出边缘的点，如下所示
for i = 1:225
    h_pos = round(hole_POS(1,i));
    w_pos = round(hole_POS(2,i));
    % 索引超出范围的索引取边缘（！注意：极小概率会出现，例如在同一点处i=1,2，则最后template_circle中就只有224个点，缺少了1）
    h_pos = min(max(1,h_pos),size(WF,1));
    w_pos = min(max(1,w_pos),size(WF,2));
    template_center(h_pos,w_pos) = i;
end
% 生成⚪模板（R*1.24=31um是实际量出来的半径）
template_circle = uint8(bwdist(template_center) < R*1.24);   % 得到距离孔心的距离图，孔内为1，孔外为0

%% 将template_circle与原始明场图像做匹配评分
score = sum(sum(WF & template_circle));

end