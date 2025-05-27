%% Impact boxplot
clear all

filename_path = 'workspace\';
filename_values = {'025','050','075','100','125','150','175','200',...
    '250','300','350','400'};
img_savepath = 'img\';
 
% titles = {'\textbf{Scenario} $a: \mu =\frac{1}{4}Q$',...
%     '\textbf{Scenario} $a: \mu =\frac{1}{2}Q$',...
%     '\textbf{Scenario} $a: \mu =\frac{3}{4}Q$',...
%     '\textbf{Scenario} $a: \mu =Q$',...
%     '\textbf{Scenario} $a: \mu =\frac{5}{4}Q$',...
%     '\textbf{Scenario} $a: \mu =\frac{3}{2}Q$',...
%     '\textbf{Scenario} $a: \mu =\frac{7}{4}Q$',...
%     '\textbf{Scenario} $a: \mu =2Q$',...
%     '\textbf{Scenario} $a: \mu =\frac{5}{2}Q$',...
%     '\textbf{Scenario} $a: \mu =3Q$',...
%     '\textbf{Scenario} $a: \mu =\frac{7}{2}Q$',...
%     '\textbf{Scenario} $a: \mu =4Q$'};
titles = {'\textbf{Scenario} $b: \mu =\frac{1}{4}Q$',...
    '\textbf{Scenario} $b: \mu =\frac{1}{2}Q$',...
    '\textbf{Scenario} $b: \mu =\frac{3}{4}Q$',...
    '\textbf{Scenario} $b: \mu =Q$',...
    '\textbf{Scenario} $b: \mu =\frac{5}{4}Q$',...
    '\textbf{Scenario} $b: \mu =\frac{3}{2}Q$',...
    '\textbf{Scenario} $b: \mu =\frac{7}{4}Q$',...
    '\textbf{Scenario} $b: \mu =2Q$',...
    '\textbf{Scenario} $b: \mu =\frac{5}{2}Q$',...
    '\textbf{Scenario} $b: \mu =3Q$',...
    '\textbf{Scenario} $b: \mu =\frac{7}{2}Q$',...
    '\textbf{Scenario} $b: \mu =4Q$'};

for h_=1:length(filename_values)
    % load(sprintf('%swrkspc_%sQ.mat',filename_path,filename_values{h_}));
    load(sprintf('%swrkspc_decentered_%sQ.mat',filename_path,filename_values{h_}));
    
    matrix_for_plot = zeros(size(Impact_matrix,1), 3*size(Impact_matrix,2));
    for i=1:size(Impact_matrix,2)
        matrix_for_plot(:, 3*i-2) = Impact_matrix(:, i);   
        matrix_for_plot(:, 3*i-1) = Impact_heuristic_matrix(:,i);
        matrix_for_plot(:, 3*i) = Impact_noBus_matrix(:, i);  
    end
    
    % boxplot
    figure()
    boxplot(matrix_for_plot, Colors="gbr", Symbol="+k", BoxStyle="filled");
    title([titles{h_}],"Interpreter","latex","FontSize", 18)
    xlabel('$\delta$', 'Interpreter','latex', "FontSize", 16)
    ylabel('Service Impact', "FontSize", 16)
    
    % label creation and positioning
    labels = cell(1, n_delta);
    for i = 0:n_delta-1
        labels{i+1} = sprintf('$\\frac{%d}{%d}\\delta_{max}$', i,n_delta-1);  
    end
    labels{1} = '0';
    labels{end} = '$\delta_{max}$';
    xticks(2:3:3*n_delta); 
    xticklabels(labels);   
    set(gca, 'TickLabelInterpreter', 'latex');
    ax = gca;
    ax.FontSize = 11;
    grid on
    
    % legend creation
    p1 = patch(NaN, NaN, [0 1 0], 'LineWidth', 1); % Green
    p2 = patch(NaN, NaN, [0 0 1], 'LineWidth', 1); % Blue
    p3 = patch(NaN, NaN, [1 0 0], 'LineWidth', 1); % Red
    legend([p1, p2, p3], {'MILP solution', 'heuristic solution', 'noBus solution'}, ...
        'Location', 'northwest', "FontSize",11);
    
    % saveas(gcf, sprintf('%sResults_%sQ.eps',img_savepath,filename_values{h_}), "epsc");
    saveas(gcf, sprintf('%sResults_decentered_%sQ.eps',img_savepath,filename_values{h_}), "epsc");
end


%% calculate the median Impact
clear all

% load the data
filename_path = 'workspace\';
filename_values = {'025','050','075','100','125','150','175','200',...
    '250','300','350','400'};
img_savepath = 'img\';

n_delta = 10;
median_Impact = zeros(n_delta,length(filename_values));
median_Impact_heuristic = zeros(n_delta,length(filename_values));
median_Impact_noBus = zeros(n_delta,length(filename_values));

for h_=1:length(filename_values)
    % load(sprintf('%swrkspc_%sQ.mat',filename_path,filename_values{h_}));
    load(sprintf('%swrkspc_decentered_%sQ.mat',filename_path,filename_values{h_}));

    median_Impact(:,h_) = median(Impact_matrix)';
    median_Impact_heuristic(:,h_) = median(Impact_heuristic_matrix)';
    median_Impact_noBus(:,h_) = median(Impact_noBus_matrix)';
end


%------------------ plot ------------------
figure()
x = [0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4]; %linspace(1,8,8);
plot(x,mean(median_Impact),"-xg", x,mean(median_Impact_heuristic),"-*b", x,mean(median_Impact_noBus),"-+r");
% title("\textbf{Scenario} $a$",'Interpreter','latex',"FontSize", 20)
title("\textbf{Scenario} $b$",'Interpreter','latex',"FontSize", 20)
xlabel('$\mu$', 'Interpreter','latex', "FontSize", 16)
ylabel('Service Impact', "FontSize", 16)
xlim([0, 4.25])
% ylim([0, 250000])
xline(1.5, LineWidth=2)
legend("MILP solution","heuristic solution","noBus solution",'Location', 'northwest',"FontSize",14)

labels = {'0','$\frac{1}{2}Q$','$Q$','$\frac{3}{2}Q$','$2Q$',...
    '$\frac{5}{2}Q$','$3Q$','$\frac{7}{2}Q$','$4Q$'};
xticklabels(labels);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = 11;
grid on

% saveas(gcf, sprintf('%sResults_mean_median.eps',img_savepath), "epsc");
saveas(gcf, sprintf('%sResults_mean_median_decentered.eps',img_savepath), "epsc");


%% Calculate the % of user served
clear all

filename_path = 'workspace\';
filename_values = {'025','050','075','100','125','150','175','200',...
    '250','300','350','400'};
img_savepath = 'img\';

n_delta = 10;
perc_user_served_MILP = zeros(n_delta,length(filename_values));
perc_user_served_heuristic = zeros(n_delta,length(filename_values));
ratio_MILP_heur = zeros(n_delta,length(filename_values));

for h_=1:length(filename_values)
    % load(sprintf('%swrkspc_%sQ.mat',filename_path,filename_values{h_}));
    load(sprintf('%swrkspc_decentered_%sQ.mat',filename_path,filename_values{h_}));

    perc_user_served_heuristic(:,h_) = mean(userServiceMatrix_heuristic(:,:,1)./...
        userServiceMatrix_heuristic(:,:,2))';
    perc_user_served_MILP(:,h_) = mean(userServiceMatrix_MILP(:,:,1)./userServiceMatrix_MILP(:,:,2))';
    ratio_MILP_heur(:,h_) = mean(userServiceMatrix_heuristic(:,:,1)./userServiceMatrix_MILP(:,:,1))';
end


%------------------ plot ------------------
figure()
x = [0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4];
hold on
plot(x,mean(perc_user_served_MILP),"-xg")
plot(x,mean(perc_user_served_heuristic),"-*b")
plot(x,mean(ratio_MILP_heur),"-+m", "LineWidth", 2);
hold off
% title("\textbf{Scenario} $a$",'Interpreter','latex',"FontSize", 20)
title("\textbf{Scenario} $b$",'Interpreter','latex',"FontSize", 20)
xlabel('$\mu$', 'Interpreter','latex', "FontSize", 16)
ylabel('%', "FontSize", 16)
xlim([0,4.25])
ylim([0,1.1])
xline(1.5, LineWidth=2)
legend("Served users: MILP","Served users: heuristic","$\rho$: heuristic/MILP",...
    'Interpreter','latex','Location','southwest',"FontSize",14)

labels = {'0','$\frac{1}{2}Q$','$Q$','$\frac{3}{2}Q$','$2Q$',...
    '$\frac{5}{2}Q$','$3Q$','$\frac{7}{2}Q$','$4Q$'};
xticklabels(labels);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = 11;
grid on

% saveas(gcf, sprintf('%sResults_mean_user_served.eps',img_savepath), "epsc");
saveas(gcf, sprintf('%sResults_mean_user_served_decentered.eps',img_savepath), "epsc");
