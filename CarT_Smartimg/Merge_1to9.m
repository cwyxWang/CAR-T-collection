function [Merge_img] = Merge_1to9(img,step)

% 计算拼接时中间需要填充的部分
lack_H = round(step/0.65) - size(img,1);   % 高度补足
lack_W = round(step/0.65) - size(img,2);   % 宽度补足
zero_H = zeros(lack_H,size(img,2));   % 高度补足(29,2048)
zero_W = zeros(size(img,1),lack_W);   % 宽度补足(2048,29)
zero_B = zeros(lack_H,lack_W);   % 插入黑块(29,29)

% 拼接
Merge_img = [img(:,:,1),zero_W,img(:,:,2),zero_W,img(:,:,3);
             zero_H    ,zero_B,  zero_H  ,zero_B,    zero_H;
             img(:,:,4),zero_W,img(:,:,5),zero_W,img(:,:,6);
             zero_H    ,zero_B,  zero_H  ,zero_B,    zero_H;
             img(:,:,7),zero_W,img(:,:,8),zero_W,img(:,:,9)];



end