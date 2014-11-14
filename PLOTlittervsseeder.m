 % Seeder
                  Lit=0:0.1:6;
                  maxGS=0.9;
                  minGS=0;
                  LitThreshS=2;
                  ampS=0.3;
                  ProbGS=(maxGS+minGS)/2+(maxGS-minGS)/2*tanh(LitThreshS-Lit/ampS);
                  plot(Lit,ProbGS)
                  xlabel('litter (cm)'), ylabel('Probability of establishment of seeder')
                  axis([0 6 0 1])