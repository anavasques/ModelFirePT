% PINE
                  
                  Lit=0:0.1:6;
                  maxGP=0.9;
                  minGP=0;
                  ProbPZeroL=0.7;
                  LitThreshP=3;
                  ampP=0.3;
                  ProbGP=(maxGP+minGP)/2+(maxGP-minGP)/2*tanh((LitThreshP-Lit)/ampP)-(maxGP-ProbPZeroL)*exp(-2/LitThreshP*exp(1)*Lit);
                  plot(Lit,ProbGP, 'b')
                  set(gca,'fontsize',20, 'fontWeight','bold');
                  title ('pine')
                  set(gca,'fontsize',20, 'fontWeight','bold');
                  xlabel('litter (cm)'), ylabel('Probability of establishment of pine')
                  set(gca,'fontsize',20, 'fontWeight','bold');
                  axis([0 6 0 1])
                  pause
%SEEDER
                  Lit=0:0.1:6;
                  maxGS=0.9;
                  minGS=0;
                  LitThreshS=2;
                  ampS=0.3;
                  ProbGS=(maxGS+minGS)/2+(maxGS-minGS)/2*tanh(LitThreshS-Lit/ampS);
                  plot(Lit,ProbGS, 'r--.')
                  set(gca,'fontsize',20, 'fontWeight','bold');
                  title('seeder')
                  set(gca,'fontsize',20, 'fontWeight','bold');
                  xlabel('litter (cm)'), ylabel('Probability of establishment of seeder')
                  set(gca,'fontsize',20, 'fontWeight','bold');
                  axis([0 6 0 1])
                  pause

% OAK                  
                  Lit=0:0.1:6;
                  maxGQ=0.9;
                  minGQ=0.3
                  ProbGQ=maxGQ-(maxGQ-minGQ)*exp(-Lit);
                  plot(Lit,ProbGQ, 'g*')
                  set(gca,'fontsize',20, 'fontWeight','bold');
                  title ('oak')
                  set(gca,'fontsize',20, 'fontWeight','bold');
                  xlabel('litter (cm)'), ylabel('Probability of establishment of oak')
                  set(gca,'fontsize',20, 'fontWeight','bold');
                  axis([0 6 0 1])
                  pause
                  
% ALL TOGETHER
plot(Lit,ProbGP, 'b', Lit,ProbGS, 'r--.', Lit,ProbGQ, 'g*')
legend('Pine','Seeder','Oak')