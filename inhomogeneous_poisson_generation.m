function [xx,yy,n_generated] = inhomogeneous_poisson_generation(r,xx0,yy0,lambda,M)

xx = [];
yy = [];
R = zeros(M,1);
for i = 1:M
    R(i) = sqrt(i/M)*r;
end
normalization = 1;
if lambda ~= 1
    sum = 0;
    for i=1:M
        sum = sum+1/i;
    end
    normalization = M/sum;
end

n_generated = zeros(M,1);
for i = 1:M
   lambda_i = lambda*normalization/i;
   A_i = pi*R(i)^2;
   if i > 1
        A_i = A_i - pi*R(i-1)^2;
   end
   numbPoints_i = poissrnd(A_i*lambda_i);
   n_generated(i) = numbPoints_i;
   theta = 2*pi*(rand(numbPoints_i,1));
   if i == 1
        rho = R(i)*sqrt(rand(numbPoints_i,1));
   else
        rho = R(i-1)+(R(i)-R(i-1))*sqrt(rand(numbPoints_i,1));
   end
   [xx_i,yy_i]=pol2cart(theta,rho);
   xx = [xx; xx_i];
   yy = [yy; yy_i];
end

xx=xx+xx0;
yy=yy+yy0;


% % --------- Plot of user distribution --------- 
% figure
% grid on
% hold on
% for i = 1:M
%     x_ = linspace(-R(i),R(i),100);
%     for j = 1:length(x_)
%         y_(j) = sqrt(R(i)^2-x_(j)^2);
%     end
%     plot(x_,y_)
%     plot(x_,-y_)
% end
% 
% scatter(xx,yy);
% xlabel('x');
% ylabel('y');
% axis square;




% % ---------  Plot user distribution on the Portello map --------- 
% img = imread('../Portello/Portello_esteso_con_fermate.png');
% 
% % Define the limits of the image
% [h, w, ~] = size(img);
% x_range = [- w/2, w/2];
% y_range = [-h/2, h/2];
% 
% % Plot the map
% figure
% imshow(img, 'XData', x_range, 'YData', y_range);
% hold on
% 
% % Define the center of the generation process in the cartesian coordinates
% scale_factor = 0.52; 
% offset_x = 65;
% offset_y = -30; 
% 
% % Plot center and circles
% for i = 0:M
%     if i ==0 
%         scatter(offset_x, offset_y, 100, 'g', 'filled'); % scenario a
%         % scatter(offset_x-scale_factor*150, offset_y+ scale_factor*125,100, 'g', 'filled'); % scenario b
%     else 
%         theta = linspace(0, 2*pi, 100);
%         x_ = offset_x + scale_factor * R(i) * cos(theta); % scenario a
%         y_ = offset_y + scale_factor * R(i) * sin(theta); % scenario a
%         % x_ = offset_x - scale_factor*150 + scale_factor * R(i) * cos(theta); % scenario b
%         % y_ = offset_y + scale_factor*125 + scale_factor * R(i) * sin(theta); % scenario b
%         plot(x_, y_, 'Color', [1 0.5 0], 'LineWidth', 1.5);
%     end
% end
% 
% % Plot users' position
% scatter(scale_factor*xx+offset_x, scale_factor*yy+offset_y, 'b', 'filled');
% 
% 
% axis image  
% hold off;
% 
% img_savepath = 'img\';
% saveas(gcf, fullfile(img_savepath, 'Scenario_a.eps'), 'epsc');
% % saveas(gcf, fullfile(img_savepath, 'Scenario_b.eps'), 'epsc');

end

