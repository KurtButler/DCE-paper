% Description:  Reproduces Figure 4 of the letter, and computes the correlation coefficient of z(t) vs estimated DCE
%       Outputs:
%           - Figure4.png in the results directory
%           - Slope of the best fit line of estimated DCE vs the signal z(t), output to the log file


N = 1000;
t = (1:N)';
y = zeros(N,1);
x = sin(t/2) + cos(sqrt(2)*t/6);
x = x + randn(size(x));
z = 1./(1+exp(15*sin(t/20)));
% x = normalize(x);
for n = 2:N
    y(n) = 0.9*y(n-1) + (0.1 + z(n-1))*x(n-1) + randn;
end

gp = fitrgp([x(1:end-1),y(1:end-1),z(1:end-1)], y(2:end),...
    'KernelFunction','ardsquaredexponential');
Mxyz = [x(1:end-1),y(1:end-1),z(1:end-1)];
a = gp.Alpha;
l = gp.KernelInformation.KernelParameters(1);
ll = gp.KernelInformation.KernelParameters(1:end-1);
sf= gp.KernelInformation.KernelParameters(end);
sg= gp.Sigma;
dkdx = (sf.^2)*exp(-0.5*(pdist2(Mxyz./ll',Mxyz./ll')).^2).*((x(1:end-1)'-x(1:end-1))/l.^2);
dFdx = dkdx*a;



%% Plot the results as a Figure
figure('Position',[80 80 590 826])

tiledlayout(3,2, 'Padding', 'compact', 'TileSpacing', 'compact'); 
nexttile(1,[1 2])
plot(t,x,t,y,'LineWidth',1);
grid on;
title('Signals','FontSize',15)
xlabel('t','FontSize',13)
grid minor;
legend('x(t)','y(t)','FontSize',13,'Location','best')

nexttile([1 2])
plot(t,z+0.1,'-',nan,nan,t(2:end),dFdx,'LineWidth',1.5);
grid on;
legend( '\partialF/\partialx = \beta+z(t)','',...
   'Estimated DCE','FontSize',13,'Location','best')
title('Modulation in time','FontSize',15)
grid minor;
xlabel('t','FontSize',13)


nexttile([1 2])
scatter(z(1:end-1), dFdx,'.')
grid on;
grid minor;
title('Differential causal effect vs z(t)','FontSize',15)
xlabel('z(t)','FontSize',13)
ylabel('Estimated \partialF/\partialx','FontSize',13)


%% Test output to the log
COVMAT = cov(z(1:end-1), dFdx);
LINEARCOEFF = COVMAT(1,2)/var(z(1:end-1));

disp(' ')
disp('Associated results for Figure 4: ')
fprintf('Slope of the estimated DCE vs the signal z(t) line: %0.3f\n',LINEARCOEFF)

%% Save the produced plot as output
saveas(gcf,'../results/Figure4.png')