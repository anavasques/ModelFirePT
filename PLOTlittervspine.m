                  Lit=0:0.1:6;
                  maxGP=0.9;
                  minGP=0;
                  ProbZeroL=0.7;
                  LitThreshP=3;
                  ampP=0.3;
                  ProbGP=(maxGP+minGP)/2+(maxGP-minGP)/2*tanh((LitThreshP-Lit)/ampP) ... 
                        -(maxGP-ProbPZeroL)*exp(-2/LitThreshP*exp(1)*Lit);
                  plot(Lit,ProbGP)
                  xlabel('litter (cm)'), ylabel('Probability of establishment of pine')
                  axis([0 6 0 1])

                  