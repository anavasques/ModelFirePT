litter=0:0.1:6;
maxQ=0.9;
minQ=0.3
ProbQ=maxQ-(maxQ-minQ)*exp(-litter);
plot(litter,ProbQ)
xlabel('litter (cm)'), ylabel('Probability')
% axis([0 6 0 1])
