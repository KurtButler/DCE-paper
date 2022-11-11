% Description:  Reproduces Figure 3 of the letter
%       Outputs:
%           - Figure3.png in the results directory


N = 500;
x= 5*rand(N,1);
z= 5*rand(N,1);
x= sort(x);
f = @(x) sin(x) + cos(2*x) + sin(3*x) + 0.1*x.^2;
df= @(x) cos(x) -2*sin(2*x)+3*cos(3*x) + 0.1*x;
y = f(x) + cos(z) + randn(N,1);

dFdx = 1+0*x;



%% SE kernel  GP
gp = fitrgp(x,y);
a = gp.Alpha;
l = gp.KernelInformation.KernelParameters(1);
sf= gp.KernelInformation.KernelParameters(end);
sg= gp.Sigma;
dkdx = (sf.^2)*exp(-0.5*(pdist2(x/l,x/l)).^2).*((x'-x)/l.^2);
dFdx = dkdx*a;

ACE = mean(dFdx);

%% Matern 52 kernel GP
gp2 = fitrgp(x,y,'KernelFunction','matern52','BasisFunction','none');
dFdx2 = [];
a   = gp2.Alpha;
sg  = gp2.Sigma;
sl  = gp2.KernelInformation.KernelParameters(1);
sf  = gp2.KernelInformation.KernelParameters(2);
r   = pdist2(x,x);
for k = 1:size(x,2)
dkdx = -(sf.^2)*((-5/3)*r./(sl^2) - 5*sqrt(5)*r.*r./(3*sl^3)).*exp(-sqrt(5)*r/sl).*(sign(x(:,k)'-x(:,k)));
dFdx2(:,k) = dkdx*a;
end

%% Matern 32 kernel GP
gp3 = fitrgp(x,y,'KernelFunction','matern32');
dFdx3 = [];
a   = gp3.Alpha;
sg  = gp3.Sigma;
sl  = gp3.KernelInformation.KernelParameters(1);
sf  = gp3.KernelInformation.KernelParameters(2);
r   = pdist2(x,x);
for k = 1:size(x,2)
dkdx = (sf.^2)*3.*exp(-sqrt(3)*r/sl).*((x(:,k)'-x(:,k))/sl^2);
dFdx3(:,k) = dkdx*a;
end

%% ARD-SE kernel GP
gp4 = fitrgp(x,y,'KernelFunction','ardsquaredexponential');
dFdx4 = [];
for k = 1:size(x,2)
a = gp4.Alpha;
ll = gp4.KernelInformation.KernelParameters(1:end-1);
sf= gp4.KernelInformation.KernelParameters(end);
sg= gp4.Sigma;
l = gp4.KernelInformation.KernelParameters(k);
dkdx = (sf.^2)*exp(-0.5*(pdist2(x./ll',x./ll')).^2).*((x(:,k)'-x(:,k))/l.^2);
dFdx4(:,k) = dkdx*a;
end




%% Evaluate GP predictions
xx = linspace(0,5,100)';
[yp,ys] = predict(gp,xx);
[yp2,~] = predict(gp2,xx);
[yp3,~] = predict(gp3,xx);
[yp4,~] = predict(gp4,xx);

%% Plot results
figure('Position',[60 60 735 740])
tiledlayout(3,1, 'Padding', 'compact', 'TileSpacing', 'compact'); 
nexttile
plot(x,y,'k.',xx,yp,'LineWidth',2)
patch([xx;flipud(xx)],[yp+ys;flipud(yp-ys)],[0.9 0.3 0],'EdgeColor','none','FaceAlpha',0.3);


grid on;
grid minor;
title('Gaussian process model (SE kernel)','FontSize',15)
legend('Observations','GP mean','GP std dev','FontSize',13,'Location','best')
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)

nexttile
plot(xx,df(xx),'k--',x,dFdx,xx,0*xx+ACE,'-.','LineWidth',2)
grid on;
grid minor;
title('DCE from SE Kernel','FontSize',15)
legend('\partialF/\partialx','GP estimate of DCE','ACE','FontSize',13,'Location','best')
xlabel('x','FontSize',13)

nexttile
plot(xx,df(xx),'k--',x,dFdx,xx,0*xx+ACE,'-.',...
    x,dFdx3,'',x,dFdx2,'','LineWidth',2)
grid on;
grid minor;
title('Comparison of Kernels','FontSize',15)
legend('\partialF/\partialx','SE Kernel','ACE',...
    'Mat3/2','Mat5/2','FontSize',13,'Location','best')
xlabel('x','FontSize',13)



%% Save the produced plot as output
saveas(gcf,'./results/Figure3.png')
