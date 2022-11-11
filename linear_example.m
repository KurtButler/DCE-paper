% Description:  Reproduces Figure 2 of the letter
%       Outputs:
%           - Figure2.png in the results directory



N = 200;

x = randn(N,4);
y = x*[1;2;3;4] + randn(N,1);


%% SE kernel
gp = fitrgp(x,y,'KernelFunction','squaredexponential');
dFdx = [];
for k = 1:size(x,2)
a = gp.Alpha;
ll = gp.KernelInformation.KernelParameters(1:end-1);
sf= gp.KernelInformation.KernelParameters(end);
sg= gp.Sigma;
l = gp.KernelInformation.KernelParameters(1);
dkdx = (sf.^2)*exp(-0.5*(pdist2(x./l,x./l)).^2).*((x(:,k)'-x(:,k))/l.^2);
dFdx(:,k) = dkdx*a;
end

beta = pinv([x,ones(size(x(:,1)))])*y;
beta(end) = [];


%% Other kernels

gp2 = fitrgp(x,y,'KernelFunction','rationalquadratic');
dFdx2 = [];
for k = 1:size(x,2)
a   = gp2.Alpha;
sg  = gp2.Sigma;
sl  = gp2.KernelInformation.KernelParameters(1);
arq = gp2.KernelInformation.KernelParameters(2);
sf  = gp2.KernelInformation.KernelParameters(3);
dkdx = (sf.^2)*((1 + (pdist2(x,x).^2)/(2*arq*sl^2)).^(-arq-1)).*((x(:,k)'-x(:,k))/l.^2);
dFdx2(:,k) = dkdx*a;
end

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


gp4 = fitrgp(x,y,'KernelFunction','ardsquaredexponential');
dFdx2 = [];
for k = 1:size(x,2)
a = gp4.Alpha;
ll = gp4.KernelInformation.KernelParameters(1:end-1);
sf= gp4.KernelInformation.KernelParameters(end);
sg= gp4.Sigma;
l = gp4.KernelInformation.KernelParameters(k);
dkdx = (sf.^2)*exp(-0.5*(pdist2(x./ll',x./ll')).^2).*((x(:,k)'-x(:,k))/l.^2);
dFdx2(:,k) = dkdx*a;
end

%% Kernel comparison
DOTSIZE = 80;
OLStarg = 'w+';
OLSSIZE = 1.5;

% Resize the fig
figure('Position',[40 40 568 377]);
tiledlayout(3,1, 'Padding', 'compact', 'TileSpacing', 'compact'); 
nexttile
edges = round(min(dFdx(:)),1):0.1:round(max(dFdx(:)),1);
for k = 1:size(x,2)
    histogram(dFdx(:,k),'EdgeAlpha',0.1,'FaceAlpha',0.5,'Normalization','pdf');
    hold on;
end
for k=1:4
    scatter(mean(dFdx(:,k)),0,DOTSIZE,'black','filled')
end
plot(beta,0*beta + 1e-2,OLStarg,'LineWidth',OLSSIZE);
hold off;
grid on;
title('Squared exponential kernel','FontSize',15);
COLORS = [        0.6350    0.0780    0.1840;
    0.9290    0.6940    0.1250  ;
        0.4660    0.6740    0.1880;
        0    0.4470    0.7410];
colororder(COLORS)





nexttile
edges = round(min(dFdx2(:)),1):0.1:round(max(dFdx2(:)),1);
for k = 1:size(x,2)
    histogram(dFdx2(:,k),'EdgeAlpha',0.1,'FaceAlpha',0.5,'Normalization','pdf');
    hold on;
end
for k=1:4
    scatter(mean(dFdx2(:,k)),0,DOTSIZE,'black','filled')
end
plot(beta,0*beta + 1e-2,OLStarg,'LineWidth',OLSSIZE);
hold off;
grid on;
title('ARD squared exponential kernel','FontSize',15);
COLORS = [        0.6350    0.0780    0.1840;
    0.9290    0.6940    0.1250  ;
        0.4660    0.6740    0.1880;
        0    0.4470    0.7410];
colororder(COLORS)



nexttile
edges = round(min(dFdx3(:)),1):0.1:round(max(dFdx3(:)),1);
for k = 1:size(x,2)
    histogram(dFdx3(:,k),'EdgeAlpha',0.1,'FaceAlpha',0.5,'Normalization','pdf');
    hold on;
end
for k=1:4
    scatter(mean(dFdx3(:,k)),0,DOTSIZE,'black','filled')
end
plot(beta,0*beta + 1e-2,OLStarg','LineWidth',OLSSIZE);
hold off;
grid on;
title('Matern 3/2 kernel','FontSize',15);
legend('\partialy/\partialx_1','\partialy/\partialx_2','\partialy/\partialx_3','\partialy/\partialx_4','FontSize',13,'Location','best')
COLORS = [        0.6350    0.0780    0.1840;
    0.9290    0.6940    0.1250  ;
        0.4660    0.6740    0.1880;
        0    0.4470    0.7410];
colororder(COLORS)



%% Save the produced plot as output
saveas(gcf,'./results/Figure2.png')

