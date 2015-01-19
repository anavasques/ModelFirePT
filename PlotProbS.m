SB=0:1:10;
m=100;
est=100;
ProbS= 1-(1-1./est).^(SB/m/m);

plot(SB,ProbS)
xlabel('Number of seeds in the seed bank'), ylabel('Probability of establishment')
axis([0 100 0 1])
