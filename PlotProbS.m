
SB=0:1:50;
% m=100; %already per cell
est=7;
ProbS= 1-(1-1./est).^(SB);
figure
plot(SB,ProbS, 'b', 'LineWidth',5)
title('pine', 'fontsize',20, 'fontWeight','bold')
set(gca,'fontsize',20, 'fontWeight','bold');
xlabel('Number of seeds in the seed bank per cell'), ylabel('Probability of establishment')
set(gca,'fontsize',20, 'fontWeight','bold');
axis([0 50 0 1])
pause

SB=0:1:500;
% m=100; %already per cell
est=100;
ProbS= 1-(1-1./est).^(SB);
figure
plot(SB,ProbS,'r--.', 'LineWidth',5)
title('seeder', 'fontsize',20, 'fontWeight','bold')
set(gca,'fontsize',20, 'fontWeight','bold');
xlabel('Number of seeds in the seed bank per cell'), ylabel('Probability of establishment')
set(gca,'fontsize',20, 'fontWeight','bold');
axis([0 500 0 1])
pause

SB=0:1:10;
% m=100; %already per cell
est=2;
ProbS= 1-(1-1./est).^(SB);
figure
plot(SB,ProbS, 'g--*', 'LineWidth',5)
title('oak', 'fontsize',20, 'fontWeight','bold')
set(gca,'fontsize',20, 'fontWeight','bold');
xlabel('Number of seeds in the seed bank per cell'), ylabel('Probability of establishment')
set(gca,'fontsize',20, 'fontWeight','bold');
axis([0 10 0 1])



pause
