
SB=0:1:50;
% m=100; %already per cell
est=7;
ProbS= 1-(1-1./est).^(SB);
figure
plot(SB,ProbS, 'b')
title('pine')
xlabel('Number of seeds in the seed bank per cell'), ylabel('Probability of establishment')
axis([0 50 0 1])
pause

SB=0:1:500;
% m=100; %already per cell
est=100;
ProbS= 1-(1-1./est).^(SB);
figure
plot(SB,ProbS,'r--.')
title('seeder')
xlabel('Number of seeds in the seed bank per cell'), ylabel('Probability of establishment')
axis([0 500 0 1])
pause

SB=0:1:10;
% m=100; %already per cell
est=2;
ProbS= 1-(1-1./est).^(SB);
figure
plot(SB,ProbS, 'g*')
title('oak' )
xlabel('Number of seeds in the seed bank per cell'), ylabel('Probability of establishment')
axis([0 10 0 1])
pause
