%% 该函数用来在预处理后的原始明场图中依次crop出225个孔的小图
function Hole_POS = WF_1Detect_1(WF,template_center)
%% 遍历225次，找出所有孔（理论也是225个）
Hole_POS = [];   % 用来存每个孔的坐标信息，该函数不更新N1。
for i = 1:225
    if length(find(template_center == i)) ~= 1    % 如果像素强度=i的点数量不为1（只可能是0或1），则i++
        continue;
    end
    temp_Hole_POS = zeros(1,7);   % (H,W,X_index,Y_index,N1,N2,N3) (H,W,1~15,1~15,1~9,1~9,1~25) 1~15主函数里换成1~45
    %% 将当前索引i(1~225)转换为N2(1~9),N3(1~25)
    % i(1~225) → x(1~15),y(1~15)
    x = mod((i-1),15) + 1;     y = ceil(i/15);
    temp_Hole_POS(3) = x;      temp_Hole_POS(4) = y;
    % x(1~15),y(1~15) → N2
    N2 = ceil(x/5) + (ceil(y/5) - 1) * 3;
    % x(1~15),y(1~15) → xx(1~5),yy(1~5) → N3
    xx = mod((x-1),5) + 1;     yy = mod((y-1),5) + 1;
    N3 = xx + (yy-1) * 5;
    temp_Hole_POS(6) = N2;
    temp_Hole_POS(7) = N3;

    %% 在孔心处crop小图（耗时0.004*225=0.9s）
    % 得到第i(1~225)个初始孔心在大图中的像素坐标(H,W)
    [h0, w0] = ind2sub(size(template_center), find(template_center == i));
    % 加强鲁棒性，当孔心在边缘时，crop超出范围的部分填黑
    WF_padding = padarray(WF,[64,64]);   % 给原图四周填64宽度的黑
    h_min = h0-64;   h_max = h0+63;   w_min = w0-64;   w_max = w0+63;   % 此处为在padding前的原图上的坐标
    % 在孔心处crop出128*128的小图
    hole_crop = WF_padding(h_min+64:h_max+64,w_min+64:w_max+64);

    %% 金字塔结构，在每个crop小图中遍历圆孔位置进行匹配（0.0020+0.0018+0.0018+0.0030=0.0086 * 225 = 2s）
    % 8倍下采样，图像为16*16，默认遍历范围为8*8,[8-3:8+4]
    [h_best,w_best] = WF_1Detect_2(hole_crop,8 ,   8    ,3,   8    ,4);   % 这里范围大了可能会报错！！！% 因为best可能在边缘，导致后面几层溢出图像索引
    % 4倍下采样
    [h_best,w_best] = WF_1Detect_2(hole_crop,4 ,h_best*2,3,w_best*2,3);
    % 2倍下采样
    [h_best,w_best] = WF_1Detect_2(hole_crop,2 ,h_best*2,2,w_best*2,2);
    % 1倍下采样（range=4实际对应9*9的滑窗）
    [h_best,w_best] = WF_1Detect_2(hole_crop,1 ,h_best*2,2,w_best*2,2);

    % (h_best,w_best)得到的是crop小图中的孔心坐标，小图的原点在原图中坐标为(h_min,w_min)，此处需变换为原图坐标
    % 鲁棒性：理论上，只要WF在阈值分割时不是大部分全黑，这里边缘的孔检测到的
    % h_best和w_best就不会是Padding的边缘，h1和w1就不会超索引
    h1 = h_best + h_min - 1;   % Y   减一是索引原因
    w1 = w_best + w_min - 1;   % X
    % 索引超出范围的索引取边缘
    temp_Hole_POS(1) = min(max(1,h1),size(WF,1));   % H(Y)
    temp_Hole_POS(2) = min(max(1,w1),size(WF,2));   % W(X)

    Hole_POS = [Hole_POS;round(temp_Hole_POS)];   % 坐标索引必须取整
end

end