close all
clear all
clc

%% GLOBAL PARAMETERS
BETA = 1/10;
Max_loss = 0.85; % max relative loss the service can tolerate
Fc = 72; % fixed cost for the service in €
Kc = 1.22; % cost for the service in €/Km
fuel_price = 1.7; % €/ Liter
vehicle_efficiency = 6 / 100; % Liter / 100Km
cost_per_Km = fuel_price * vehicle_efficiency;
NMC = 10;


%% POLLUTION DIAGRAMS (presence of electric vehicles)
elec_values = [0 0.25 0.5 0.75 1];
elec_cases = length(elec_values);
lambda_values = [0.05, 0.1, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4];
lambda_cases = length(lambda_values);

% Legnaro simulations

pollL = zeros(elec_cases, lambda_cases, 2);
pollL = get_pollution_data(pollL,'./sim-ele-L/SL_b','workspace1L',...
    elec_cases,lambda_cases,elec_values,lambda_values);

ftsz = 20;
lw = 1.5;

figure % shuttle is not electric
hh = [];
for j = 1:elec_cases
    h = semilogy(lambda_values,pollL(j,:,1)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollL(j,:,1)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz)
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
ylabel('$\mathcal{I}^\star$','FontSize',ftsz,'Interpreter','Latex')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex'; 
yticks([0.5 1 2 3 4 5])
yticklabels([0.5 1 2 3 4 5])
ylim([0.3 5]) 
%title('Legnaro: Impact, normal bus')


figure % Users: shuttle is not electric
hh = [];
for j = 1:elec_cases
    h = plot(lambda_values,pollL(j,:,3)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollL(j,:,3)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz,'location','southeast')
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
yticks(0:10:100)
ylabel('$P(v^\star)$','FontSize',ftsz,'Interpreter','Latex')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';   
ylim([0 100])
%title('Legnaro: Users, normal bus')


figure % Distance ratio: shuttle is not electric
hh = [];
for j = 1:elec_cases
    h = semilogy(lambda_values,pollL(j,:,5)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollL(j,:,5)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz)
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
ylabel('shuttle distance / nonserved total travel')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';   
ylim([1e-3 1])
title('Legnaro: Distance ratio, normal bus')


figure % shuttle is electric
hh = [];
for j = 1:elec_cases
    h = semilogy(lambda_values,pollL(j,:,2)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollL(j,:,2)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz)
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
ylabel('$\mathcal{I}^\star$','FontSize',ftsz,'Interpreter','Latex')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';   
yticks([0.5 1 2 3 4 5])
yticklabels([0.5 1 2 3 4 5])
ylim([0.3 5])
%title('Legnaro: Impact, electric bus')


figure % Users: shuttle is electric
hh = [];
for j = 1:elec_cases
    h = plot(lambda_values,pollL(j,:,4)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollL(j,:,4)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz,'location','southeast')
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
yticks(0:10:100)
ylabel('$P(v^\star)$','FontSize',ftsz,'Interpreter','Latex')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';   
ylim([0 100])
%title('Legnaro: Users, electric bus')


figure % Distance ratio: shuttle is electric
hh = [];
for j = 1:elec_cases
    h = semilogy(lambda_values,pollL(j,:,6)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollL(j,:,6)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz)
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
ylabel('shuttle distance / nonserved total travel')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';   
ylim([1e-3 1])
title('Legnaro: Distance ratio, normal bus')


%%%%%%%%%%%%%%%%

% Padova simulations

pollP = zeros(elec_cases, lambda_cases, 2);
pollP = get_pollution_data(pollP,'./sim-ele-P/SP_b','workspace1L',...
    elec_cases,lambda_cases,elec_values,lambda_values);

%%
figure % shuttle is not electric
hh = [];
for j = 1:elec_cases
    h = semilogy(lambda_values,pollP(j,:,1)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollP(j,:,1)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz)
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
ylabel('$\mathcal{I}^\star$','FontSize',ftsz,'Interpreter','Latex')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';  
yticks([0.5 1 2 3 4 5])
yticklabels([0.5 1 2 3 4 5])
ylim([0.3 5])
%title('Padova: Impact, normal bus')
%%


figure % Users: shuttle is not electric
hh = [];
for j = 1:elec_cases
    h = plot(lambda_values,pollP(j,:,3)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollP(j,:,3)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz,'location','southeast')
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
yticks(0:10:100)
ylabel('$P(v^\star)$','FontSize',ftsz,'Interpreter','Latex')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';   
ylim([0 100])
%title('Padova: Users, normal bus')


figure % Distance ratio: shuttle is not electric
hh = [];
for j = 1:elec_cases
    h = semilogy(lambda_values,pollP(j,:,5)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollP(j,:,5)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz)
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
ylabel('shuttle distance / nonserved total travel')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';   
ylim([1e-3 1])
title('Padova: Distance ratio, normal bus')


figure % shuttle is electric
hh = [];
for j = 1:elec_cases
    h = semilogy(lambda_values,pollP(j,:,2)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollP(j,:,2)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz)
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
ylabel('$\mathcal{I}^\star$','FontSize',ftsz,'Interpreter','Latex')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';   
yticks([0.5 1 2 3 4 5])
yticklabels([0.5 1 2 3 4 5])
ylim([0.3 5])
%title('Padova: Impact, electric bus')


figure % Users: shuttle is electric
hh = [];
for j = 1:elec_cases
    h = plot(lambda_values,pollP(j,:,4)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollP(j,:,4)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz,'location','southeast')
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
yticks(0:10:100)
ylabel('$P(v^\star)$','FontSize',ftsz,'Interpreter','Latex')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';   
ylim([0 100])
%title('Padova: Users, electric bus')


figure % Distance ratio: shuttle is electric
hh = [];
for j = 1:elec_cases
    h = semilogy(lambda_values,pollP(j,:,6)','Color',clr(j),'LineWidth',lw);
    hold on 
    grid on
    hh = [hh h];
    scatter(lambda_values,pollP(j,:,6)',50,clr(j),"filled")
end
legend(hh,'$0$','$25\%$','$50\%$','$75\%$','$100\%$','Interpreter','Latex','FontSize',ftsz)
xlabel('$\mu$','Interpreter','Latex','FontSize',ftsz)
xticks([0 0.5 1 1.5 2 2.5 3 3.5 4])
xticklabels({'$0$','$\frac{1}{2} Q$','$Q$','$\frac{3}{2} Q$','$2Q$','$\frac{5}{2} Q$','$3Q$','$\frac{5}{2} Q$','$4 Q$'})
ylabel('shuttle distance / nonserved total travel')
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex';   
ylim([1e-3 1])
title('Padova: Distance ratio, normal bus')

%error('STOP')

%% PARETO DIAGRAMS

ticket_values = [0.5 1 2];
ticket_cases = length(ticket_values);

paretoL = zeros(ticket_cases, lambda_cases, 2);
[paretoL, pricesL] = get_pareto_data(paretoL,'L',ticket_cases,lambda_cases,ticket_values,lambda_values);

paretoP = zeros(ticket_cases, lambda_cases, 2);
[paretoP, pricesP] = get_pareto_data(paretoP,'P',ticket_cases,lambda_cases,ticket_values,lambda_values);

DD = 75;
ftsz = 1.4*ftsz;

% Legnaro simulations 
%%
figure % impact
grid on
hold on
hh = [];
for i = 1:lambda_cases
    h = plot(paretoL(:,i,1),paretoL(:,i,2),'Color',clr(i),'LineWidth',lw);
    hh = [hh h];
    for j = 1:ticket_cases
        sdim = DD*(1+j);
        scatter(paretoL(j,i,1),paretoL(j,i,2),sdim,clr(i),"filled")
        xscale("log")
    end
    %text(1,paretoL(end,i,2),[' $\quad\lambda = ' num2str(lambda_values(i)) '$'],...
    %    'Interpreter','Latex')
end
lgd = legend(hh,'$\lambda = 0.05$','$\lambda = 0.1$','$\lambda = 0.25$','$\lambda = 0.5$','$\lambda = 1$',...
    '$\lambda = 1.5$','$\lambda = 2$','$\lambda = 2.5$','$\lambda = 3$','$\lambda = 3.5$','$\lambda = 4$',...
    'Location','southwest','Interpreter','Latex','FontSize',ftsz);
lgd.NumColumns = 6;
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex'; 
xlim([0.1 12])
ylim([0.3 2])
xlabel('$\mathcal{L}^\star$','Interpreter','Latex','FontSize',ftsz)
ylabel('$\mathcal{I}^\star$','Interpreter','Latex','FontSize',ftsz)
%title('Legnaro: Pareto diagrams')
%%

figure % users
grid on
hold on
hh = [];
for i = 1:lambda_cases
    h = plot(paretoL(:,i,1),paretoL(:,i,3),'Color',clr(i),'LineWidth',lw);
    hh = [hh h];
    for j = 1:ticket_cases
        sdim = DD*(1+j);
        scatter(paretoL(j,i,1),paretoL(j,i,3),sdim,clr(i),"filled")
        xscale("log")
    end
    %text(1,paretoL(end,i,3),[' $\quad\lambda = ' num2str(lambda_values(i)) '$'],...
    %    'Interpreter','Latex')
end
legend(hh,'$\lambda = 0.05$','$\lambda = 0.1$','$\lambda = 0.25$','$\lambda = 0.5$','$\lambda = 1$',...
    '$\lambda = 1.5$','$\lambda = 2$','$\lambda = 2.5$','$\lambda = 3$','$\lambda = 3.5$','$\lambda = 4$',...
    'Location','northeast','Interpreter','Latex','FontSize',ftsz)
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex'; 
xlim([0.1 12])
ylim([0 100])
yticks(0:10:100)
xlabel('$\mathcal{L}^\star$','Interpreter','Latex','FontSize',ftsz)
ylabel('$P(v^\star)$','Interpreter','Latex','FontSize',ftsz)
%title('Legnaro: Pareto diagrams')
%%

figure % distance ratio
grid on
hold on
hh = [];
for i = 1:lambda_cases
    h = plot(paretoL(:,i,1),paretoL(:,i,4),'Color',clr(i));
    hh = [hh h];
    for j = 1:ticket_cases
        sdim = DD*(1+j);
        scatter(paretoL(j,i,1),paretoL(j,i,4),sdim,clr(i),"filled")
        xscale("log")
    end
    text(1,paretoL(end,i,4),[' $\quad\lambda = ' num2str(lambda_values(i)) '$'],...
        'Interpreter','Latex')
end
xlabel('loss []')
ylabel('shuttle distance / nonserved total travel')
title('Legnaro: Pareto diagrams')
%%

%%%%%%%%%%%%%%%

% Padova simulations
figure
grid on
hold on
hh = [];
for i = 1:lambda_cases
    h = plot(paretoP(:,i,1),paretoP(:,i,2),'Color',clr(i),'LineWidth',lw);
    hh = [hh h];
    for j = 1:ticket_cases
        sdim = DD*(1+j);
        scatter(paretoP(j,i,1),paretoP(j,i,2),sdim,clr(i),"filled")
        xscale("log")
    end
    %text(1,paretoP(end,i,2),[' $\quad\lambda = ' num2str(lambda_values(i)) '$'],...
    %    'Interpreter','Latex')
end
lgd = legend(hh,'$\lambda = 0.05$','$\lambda = 0.1$','$\lambda = 0.25$','$\lambda = 0.5$','$\lambda = 1$',...
    '$\lambda = 1.5$','$\lambda = 2$','$\lambda = 2.5$','$\lambda = 3$','$\lambda = 3.5$','$\lambda = 4$',...
    'Location','southwest','Interpreter','Latex','FontSize',ftsz);
lgd.NumColumns = 6;
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex'; 
xlim([0.1 12])
ylim([0.3 2])
xlabel('$\mathcal{L}^\star$','Interpreter','Latex','FontSize',ftsz)
ylabel('$\mathcal{I}^\star$','Interpreter','Latex','FontSize',ftsz)
%title('Padova: Pareto diagrams')
%%

figure % users
grid on
hold on
hh = [];
for i = 1:lambda_cases
    h = plot(paretoP(:,i,1),paretoP(:,i,3),'Color',clr(i),'LineWidth',lw);
    hh = [hh h];
    for j = 1:ticket_cases
        sdim = DD*(1+j);
        scatter(paretoP(j,i,1),paretoP(j,i,3),sdim,clr(i),"filled")
        xscale("log")
    end
    %text(1,paretoP(end,i,3),[' $\quad\lambda = ' num2str(lambda_values(i)) '$'],...
    %    'Interpreter','Latex')
end
legend(hh,'$\lambda = 0.05$','$\lambda = 0.1$','$\lambda = 0.25$','$\lambda = 0.5$','$\lambda = 1$',...
    '$\lambda = 1.5$','$\lambda = 2$','$\lambda = 2.5$','$\lambda = 3$','$\lambda = 3.5$','$\lambda = 4$',...
    'Location','northeast','Interpreter','Latex','FontSize',ftsz)
ax = gca;         
ax.FontSize = ftsz;       
ax.TickLabelInterpreter = 'latex'; 
xlim([0.1 12])
ylim([0 100])
yticks(0:10:100)
xlabel('$\mathcal{L}^\star$','Interpreter','Latex','FontSize',ftsz)
ylabel('$P(v^\star)$','Interpreter','Latex','FontSize',ftsz)
%title('Padova: Pareto diagrams')

%%

figure % distance ratio
grid on
hold on
for i = 1:lambda_cases
    plot(paretoP(:,i,1),paretoP(:,i,4),'Color',clr(i))
    for j = 1:ticket_cases
        sdim = DD*(1+j);
        scatter(paretoP(j,i,1),paretoP(j,i,4),sdim,clr(i),"filled")
        xscale("log")
    end
    text(1,paretoP(end,i,4),[' $\quad\lambda = ' num2str(lambda_values(i)) '$'],...
        'Interpreter','Latex')
end
xlabel('loss []')
ylabel('shuttle distance / nonserved total travel')
title('Padova: Pareto diagrams')


%% 
function diag = get_pollution_data(diag,path_1,foldername,...
    elec_cases,lambda_cases,elec_values,lambda_values)

for bb = 0:1
    path_2 = [path_1 num2str(bb) '_ev'];
    for jj = 1:elec_cases
        path_3 = [path_2 num2str(floor(100*elec_values(jj))) '/' ...
            foldername '/'];
        
        for ii = 1:lambda_cases

            zz = '';
            if ii == 1
                zz = '00';
            elseif ii == 2 || ii == 3 || ii == 4
                zz = '0';
            end

            filename = ['wrkspc_' foldername(end) '_' zz ... 
                num2str(floor(100*lambda_values(ii))) 'Q.mat'];
            path_4 = [path_3 filename];
            load(path_4,"add_info_MILP")
            load(path_4,"userServiceMatrix_MILP")

            % y-axis
            P_v_star = userServiceMatrix_MILP(:,1,1);
            beta = 162/718;
            if bb == 1
                beta = beta/0.25;
            end
            Dz0_stars = info_(P_v_star, add_info_MILP(:,1,1));
            M_stars = info_(P_v_star, add_info_MILP(:,1,2));
            M_0s = info_(P_v_star, add_info_MILP(:,1,3));
            Dv_stars = info_(P_v_star, add_info_MILP(:,1,15));
            % Impact ratio
            diag(jj,ii,bb+1) = mean((Dz0_stars+beta*M_stars)./(beta*M_0s));
            % Total passegners
            diag(jj,ii,bb+3) = mean(P_v_star);
            % Ratio: shuttle bus Km / sum(Km of unserved users vehicles)
            diag(jj,ii,bb+5) = mean(Dz0_stars ./ Dv_stars);

            clear add_info_MILP
            clear userServiceMatrix_MILP
        end
    end
end

end



function [diag,prices] = get_pareto_data(diag,scenario,...
    ticket_cases,lambda_cases,ticket_values,lambda_values)

    path_1 = './S';
    prices = zeros(lambda_cases,1);

    for jj = 1:ticket_cases
        ss = '05';
        switch ticket_values(jj)
            case 1
                ss = '1';
            case 2
                ss = '2';
        end
        path_2 = [path_1 scenario '_' ss '/workspace1L/'];
        for ii = 1:lambda_cases

            zz = '';
            if ii == 1
                zz = '00';
            elseif ii == 2 || ii == 3 || ii == 4
                zz = '0';
            end

            filename = ['wrkspc_L_' zz ... 
                num2str(floor(100*lambda_values(ii))) 'Q.mat'];
            path_3 = [path_2 filename];


            load(path_3,"add_info_MILP")
            load(path_3,"userServiceMatrix_MILP")
            load(path_3,"PRICE")
            load(path_3,"Fc")
            load(path_3,"Kc")
            

            if jj == 2
                load(path_3,"PRICE_REF")
                prices(ii) = PRICE_REF;
                lambda = lambda_values(ii);
                lambda_vs_PRICE_REF = [lambda PRICE_REF];
                lambda_vs_PRICE_REF
                clear PRICE_REF
            end

            % y-axis
            P_v_star = userServiceMatrix_MILP(:,1,1);
            Dz0_stars = info_(P_v_star, add_info_MILP(:,1,1));
            M_stars = info_(P_v_star, add_info_MILP(:,1,2));    
            M_0s = info_(P_v_star, add_info_MILP(:,1,3));
            Dv_stars = info_(P_v_star, add_info_MILP(:,1,15));
            beta = 162/718;
            
            diag(jj,ii,2) = mean((Dz0_stars+beta*M_stars)./(beta*M_0s));
            diag(jj,ii,3) = mean(P_v_star);
            diag(jj,ii,4) = mean(Dz0_stars ./ Dv_stars);

            % x-axis: loss
            Dz0_star = Dz0_stars/1000;
            P_v_star = info_(P_v_star, userServiceMatrix_MILP(:,1,1));
            one = ones(size(P_v_star));
            %PRICE = PRICE*ticket_values(jj);
            diag(jj,ii,1) = 0.3*mean( one - ...
                (P_v_star*PRICE-Fc*one-Kc*Dz0_star)./...
                (P_v_star*PRICE) );

            clear add_info_MILP
            clear userServiceMatrix_MILP
            clear PRICE
            clear Fc
            clear Kc
        end
    end



end



function v = info_(P_v_star,v)

counter = length(P_v_star);
while counter > 0
    if P_v_star(counter) == 0
        v(counter) = [];
    end
    counter = counter - 1;
end

end



function my_color = clr(n)

switch n
    case 1
        my_color = [220,20,60]; % crimson
    case 2
        my_color = [255,215,0]; % gold
    case 3
        my_color = [0,128,0]; % green
    case 4
        my_color = [30,144,255]; % dodger blue
    case 5
        my_color = [210,105,30]; % chocolate
    case 6
        my_color = [0,0,255]; % blue
    case 7
        my_color = [255,20,147]; % deep pink
    case 8
        my_color = [105,105,105]; % dim gray
    case 9
        my_color = [0,255,0]; % lime
    case 10
        my_color = [255,0,255]; % magenta
    case 11
        my_color = [0,0,0]; % black
end

my_color = my_color / 255;

end