                  
% OAK                  
                  Lit=0:0.1:6;
                  maxGQ=0.9;
                  minGQ=0.3
                  ProbGQ=maxGQ-(maxGQ-minGQ)*exp(-Lit);
                  plot(Lit,ProbGQ)
                  xlabel('litter (cm)'), ylabel('Probability of establishment of oak')
                  axis([0 6 0 1])