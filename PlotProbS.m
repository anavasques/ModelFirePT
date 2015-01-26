% SB=0:1:1000;
% % m=100;
% est=100;
% ProbS= 1-(1-1./est).^(SB);
% 
% figure
% plot(SB,ProbS)
% title('pines' )
% xlabel('Number of seeds in the seed bank per cell'), ylabel('Probability of establishment')
% axis([0 1000 0 1])

SB=0:1:10;
% m=100;
est=2;
ProbS= 1-(1-1./est).^(SB);

figure
plot(SB,ProbS)
title('oak' )
xlabel('Number of seeds in the seed bank per cell'), ylabel('Probability of establishment')
axis([0 10 0 1])
