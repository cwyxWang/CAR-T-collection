%% 该函数用来遍历x0,y0,theta的滑窗范围
function [template_center,template_circle,x_best,y_best,theta_best] = WF_225Detect_1(WF,scale,x0,x_range,y0,y_range,xy_step,theta,theta_range,theta_step) 

WF_downsampling = imresize(WF,round(size(WF) / scale),"nearest");   % 将预处理后的原始明场图像进行下采样
WF_downsampling = imbinarize(WF_downsampling);   % 下采样后再二值化，使用自适应阈值
best_score = -1;
for temp_theta = (theta - theta_range) : theta_step : (theta + theta_range)
    for temp_x0 = (x0 - x_range) : xy_step : (x0 + x_range)
        for temp_y0 = (y0 - y_range) : xy_step : (y0 + y_range)
            % 模板匹配得分
            [temp_score,temp_template_center,temp_template_circle] = WF_225Detect_2(WF_downsampling,scale,temp_x0,temp_y0,temp_theta);
            if temp_score > best_score
                best_score = temp_score;
                x_best = temp_x0; y_best = temp_y0; theta_best = temp_theta;
                template_center = temp_template_center;  template_circle = temp_template_circle;
            end
        end
    end
end

end