% Description:  Reproduces Figure 1 of the letter, and computes the regression coefficients for each corresponding model. 
%       Outputs:
%           - Figure1.png in the results directory
%           - Regression coefficients for each model, printed in the log file.


%% Data generative model
N = 500;
z = ceil(5*rand(N,1));
x = 2*z + randn(N,1);
y = 4*z -x  + randn(N,1);


%% Produce two plots to visualize
figure
tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact'); 

nexttile
scatter(x,y,[],[0 0 0],'filled');
grid on;
grid minor;
title('','p(x,y)','FontSize',17);
xlabel('x','FontSize',15);
ylabel('y','FontSize',15);

nexttile
for i = 1:5
    scatter(x(z==i), y(z==i),'filled');
    hold on;
end
hold off;
grid on;
grid minor;
colororder(hsv(5))
caxis([0.5, 5.5])
title('','p(x,y|z)','FontSize',17);
xlabel('x','FontSize',15);
ylabel('y','FontSize',15);
legend('z=1','z=2','z=3','z=4',"z=5",'FontSize',15,'Location','best');

sgtitle("Simpson's Paradox",'FontSize',17);

%% Check the sign of the regression coefficients
a1 = x\y;
ab = [x,z]\y;
a2 = ab(1);

disp('Associated results for Figure 1: ')
disp('Linear model coefficients with and without confounder knowledge')
fprintf('If y = ax + eps, then a = %0.2f.\n',a1)
fprintf('If y = ax + bz + eps, then a = %0.2f.\n',a2)


%% Save the produced plot as output
saveas(gcf,'../results/Figure1.png')
