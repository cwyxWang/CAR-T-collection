%% 该函数用来在每一个crop的小图上遍历孔心坐标，并与圆孔进行模板匹配
function [h_best,w_best] = WF_1Detect_2(hole_crop,scale,h0,h_range,w0,w_range)

if scale ~= 1
    % 将预处理后的原始明场图像进行下采样
    hole_crop = imresize(hole_crop,round(size(hole_crop) / scale),"nearest"); % 输出为uint8格式，最近邻插值不改变像素强度
end
best_score = -1;

%% 在crop小图中遍历孔心位置，进行滑窗
% 根据输入的范围滑窗
for temp_h0 = (h0 - h_range) : (h0 + h_range)
    for temp_w0 = (w0 - w_range) : (w0 + w_range)
        temp_circle = zeros(size(hole_crop));
        temp_circle(temp_h0,temp_w0) = 1;
        % 生成⚪模板（圆环、渐变权重、外圈负权重效果均不好）
        temp_circle = bwdist(temp_circle) < (31/0.65/scale);   % 31um是图上量出来的孔半径，10X明场下0.65um/pixel
        temp_score = sum(sum(hole_crop & uint8(temp_circle)));
        if temp_score > best_score
            best_score = temp_score;   h_best = temp_h0;   w_best = temp_w0;
        end
    end
end


end