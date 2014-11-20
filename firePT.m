%Model PT fire 
%CASCADE project
%Asynchronous CA
%Vasques et al. 2014 (model development)
% ------------------------------------------------------------------------

close all
clear all

%Parameter values
%-------------------------------------------------------------------------
%PER SPECIES

%%% PINE
AgeMP=10;                 % Age of maturity pine %start at 6 and regularly 10-15 % in Cronk and Fuller, 1995 [year] 
SeedFP=945;               % Seed production per pine mature tree number of seeds per
                          % cone (63)* cone per tree (15) Vega et al 2008
LSP= 100;                 % Life span of pine % in "practices centro pinus" [year]
canopyBank=0.5;           % Percent of the seeds that are stored in the canopy and released with fire
ReleaseSeeds=0;           % Pine seeds in the canopy that are released after the fire

%%% SEEDER
AgeMS=1;                  % Age of maturity seeder % field obs Calluna% [year]
                          % Cistus 3 years ref
SeedFS=100;               % Seed production per plant/occupied cell approx value ADJUST
                          % check in lit
LSS=30;                   % Life span calluna % in woodland education centre [year]

%%% OAK
AgeMO=50;                 % Age of maturity seeder % Kew % [year]
                          % !!! Pausas 1999 has maturity = 15!!!
                 
SeedFQ=10;                % Seed production oak per occupied cell - value that is not fundamented by literature assumes that the tree is dominated and small if under cover of pinus
                          % 120 acorns per tree refered in Martin?k et al. 2014% [n/m2/year] 
BirdSeedN=50;             % Annual seed input by birds - based on average values Q. suber Pons and Pausas 2007 - this value depends on surrounding populations
%BirdSeedN=200; 
LSO= 1000;                % Life span quercus robur % in forestar

%GENERAL
minG= [0 0 0.3];          % minumum germination first pine second seeder third oak
maxG= [0.9 0.9 0.9];      % maximum germination first pine second seeder third oak

ProbPZeroL=0.7;           % Germination probability for pine when litter=0 cm  

LitThreshP=3;             % Litter threshold for Pine above which ~no germination (cm)
LitThreshS=2;             % Litter threshold for seeders above which ~no germination (cm)
SeedLoss= [1 0.10 1];     % rate seed loss first pine second seeder third oak

lrate=0.42;               % rate of litter deposition [cm/year] Fernandes et al 2004
ProbG=[0 0 0];            % probability of germination first pine second seeder third oak
amp=[0.3 0.3 0];          % amplitude of curve interaction with litter

SBP1=0;                   % !SBP1 and SBP2 are only ways of initializing the seed bank every year
SBP2=0;
SBPC=0;                   % Inicialization seed bank pine canopy
SB= [0 0 0];              % seed bank first pine second seeder third oak
ProbS= [0 0 0];           % to calculate probability based on seed prod
est= [7 400];             % max number of seedlings per cell CCD field from which we infered a probability of establishment in one cell
%nrsp=3                   % number of species used in the model - to put in the prob expression

mort=1./[LSP,LSS,LSO];    % MORTALITY of pine, seeder, oak = 1/lifespan

AR= [0,0,0,1];            % ability to resprout: first element is fake (bare soil); pine=0, seeder=0, oak=1;

% Control constants
StartTime= 0;             % [year]
EndTime= 100;             % [year]
StoreStep = 1;            % [year]

% NOTE: put dt smaller than one year in a way that the probabilities are <1 but not too small otherwise the model runs slowly
% IMPROVE!!! dt needs to be smaller or change probab to max 3
dt=0.5;                 % improve!![year] very small dt for calculating probabilities
m= 100;                   % for size of lattice [meter]

% Control Variables
NrStore = 1;
StoreTime = StoreStep;
Time = StartTime;

D=0;                      % initialization only
Pine=0;                   % will count the number of cells with pine
Seeder=0;                 % will count the number of cells with pine
Oak=0;                    % will count the number of cells with pine

% Later ptr will be added
% experiments with management and fire will be done
% experiments with spatially explicit seed dispersal can be done
% experiments with microclimate adaptation at the germination level can be
% done

%%%Initialization of the matrices 
%TC - Type of Cover: 1- pine; 2- seeder; 3- oak; Age- plant age; SB- Seed Bank); Lit- Litter
% -------------------------------------------------------------------------
TC= zeros(m,m);           % Creates a matrix of size m*m filled with zeros       
Age= zeros(m,m);          % Creates a matrix of size m*m filled with zeros    
Lit= zeros(m,m);          % Creates a matrix of size m*m filled with zeros
z= 8;                     % Number of neighbours

% Fills the matrix TC with planted pines
%--------------------------------------------------------------------------
TC(4:4:m-4,4:4:m-4)= 1; % plants 1 pine every 4 meters - dense prodution stand excluding the borders

% Puts seeds in the matrix
%--------------------------------------------------------------------------
SB=[100 200 0+randi(BirdSeedN,1)]; %changing initial conditions !now same number of cells per seeder and oak

if SB(3)>0
    coordseed=randi(m,SB(3),2);
end
% Creates colormap
figure
white=[1 1 1];
green=[0 1 0];
red=[1 0 0];
blue=[0 0 1];
VegetationColormap=[white; green; red; blue];
% Plot image
h=subplot(1,1,1);
imagesc(TC)
set(h,'Clim',[-0.5 3.5]);
colormap(VegetationColormap);
colorbar


% colorbar; set(gco,'Clim',[1 4]);
%%%%%%%%%%%%%%%%%%%%%DYNAMIC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
while Time < EndTime
    for k=1:round(1/dt)
        % Creates litter in the neighborhod of pine (8 neighbors)+ the pine
        % site itself
        [x,y]=find(TC(2:end-1,2:end-1)==1);
        x=x+1;y=y+1;
        for i=1:length(x)
            Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)=Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)+lrate*dt;
% !!! Consider adding litter in the neighbourhood of oak when it is dominant
% and adult
        end
        
         % SEED BANK CALCULATION ONLY ONCE A YEAR
    SB(1)=SBP1+SBP2+SeedFP*(1-canopyBank)*sum(sum(TC(Age>AgeMP)==1))+ReleaseSeeds;% TWO YEARS OF SEED LIFE; 1-canopybank is doing the same as canopy bank, i.e. *0.5
%     SB(1)=SB(1)-SeedLoss(1)*SB(1);      %!!This term should be excluded
%     or corrected!!
    SB(2)=SB(2)+SeedFS*(sum(sum(TC==2)))-SeedLoss(2)*SB(2);           % LONG SEED LIFE
    SB(3)=SeedFQ*(sum(sum(TC(Age>AgeMO)==3)))+randi(BirdSeedN,1);     % NO MEMORY - !!!Before seed loss was 1 now there is no seed loss
    if SB(3)>0
   coordseed=randi(m,SB(3),2); % puts the seeds that arrive in random coordinates of the lattice
    end
    
    
   
    % accumulation of seeds in the canopy
    SBPC=SeedFP*canopyBank*sum(TC(Age>AgeMP)==1); %canoyBank=% of the seeds that stay in the
    % canopy and accumulate over time
        
        % Colonization of an empty cell & mortality for vegetated cells
        %------------------------------------------------------------------
        % At each time step tests who is going to colonize an empty cell in the lattice based on the
        % existent seed bank
        
        
        for i = 1 : m
            for j=1:m
                test=rand;
                if TC(i,j)==0 % colonization/germination
                    
                    % LITTER DEPENDENCE FOR COLONIZATION:
                    
% !! check this and maybe change it to sum total and divide by nrsp - should it be multiplying both terms of prob (ProbG and ProbS)? Understand what it does exactly
                    ProbG=ProbG*dt;                     %this is the trick to get probability small
                    ProbG(1)=(maxG(1)+minG(1))/2+(maxG(1)-minG(1))/2*tanh((LitThreshP-Lit(i,j))/amp(1)) ... 
                        -(maxG(1)-ProbPZeroL)*exp(-2/LitThreshP*exp(1)*Lit(i,j)); % PINE
                    ProbG(2)=(maxG(2)+minG(2))/2+(maxG(2)-minG(2))/2*tanh((LitThreshS-Lit(i,j))/amp(2)); % SEEDER ampS=0.3 max=.9 min=0.
                    ProbG(3)=maxG(3)-(maxG(3)-minG(3))*exp(-Lit(i,j)); % QUERCUS
                    
                    % Term to plot the sp relation with litter
                    % independently
                    % Pine
%                   Lit=0:0.1:6;
%                   maxGP=0.9;
%                   minGP=0;
%                   ProbZeroL=0.7;
%                   LitThreshP=3;
%                   ampP=0.3;
%                   ProbGP=(maxGP+minGP)/2+(maxGP-minGP)/2*tanh((LitThreshP-Lit)/ampP) ... 
%                         -(maxGP-ProbPZeroL)*exp(-2/LitThreshP*exp(1)*Lit);
%                   plot(Lit,ProbGP)
%                   xlabel('litter (cm)'), ylabel('Probability')
%                   axis([0 6 0 1])
%                     
                    % Seeder
%                   Lit=0:0.1:6;
%                   maxGS=0.9;
%                   minGS=0;
%                   LitThreshS=2;
%                   ampS=0.3;
%                   ProbGS=(maxGS+minGS)/2+(maxGS-minGS)/2*tanh(LitThreshS-Lit/ampS);
%                   plot(Lit,ProbGS)
%                   xlabel('litter (cm)'), ylabel('Probability')
%                   axis([0 6 0 1])
                    
                    % Oak
%                   Lit=0:0.1:6;
%                   maxGQ=0.9;
%                   minGQ=0.3
%                   ProbGQ=maxGQ-(maxGQ-minGQ)*exp(-Lit);
%                   plot(Lit,ProbGQ)
%                   xlabel('litter (cm)'), ylabel('Probability')
%                   axis([0 6 0 1])
                  

                    %%% Term for the relation between available seeds and Prob S.                      

                    ProbS(1:2)=1-(1-1./est(1:2)).^(SB(1:2)/m/m); % FOR PINE AND SEEDERS, SEEDS ARE EQUALLY SPREAD THROUGHOUT THE CELLS; this was taken in the paper: Cannas et al. 2003
                    ProbS(3)=0;                                   %!!!check what this part does exactly!!!                                
                    for ii=1:size(coordseed,1)
                        ProbS(3)=ProbS(3)+(coordseed(ii,1)==i&coordseed(ii,2)==j);
                    end
                    ProbS(3)=ProbS(3)>1; % !!!check what it does exactly!!!
                    ProbG=ProbG.*ProbS;
                    
% this step has the reference  of Alains' MSc thesis ->
%CAN BE IMPROVED
                    ProbG=ProbG*dt; %this is the trick to get probability small
                    if sum(ProbG)>1
                        'sum of probability higher than 1! Please decrease dt'
                        break
% !!! check how changes in dt affect (or not) the age etc

                    elseif test<ProbG(1)
                        TC(i,j)=1;
                        Age(i,j)=dt;
                        
                    elseif test<ProbG(1)+ProbG(2)
                        TC(i,j)=2;
                        Age(i,j)=dt;
                        
                    elseif test<ProbG(1)+ProbG(2)+ProbG(3)
                        TC(i,j)=3;
                        Age(i,j)=dt;
                    end
                else
                    Age(i,j)=Age(i,j)+dt;
                    
                    if test< mort(TC(i,j))*dt
                        TC(i,j)=0;
                        Age(i,j)=0;
% !! term for litter when a pine dies (same cell) - same or not?
                    end
                end
            end
            SBP2=SBP1; % PINE SEED BANK OF TWO YEARS BEFORE
            SBP1=SB(1);% PINE SEED BANK OF 1 YEAR BEFORE
        end
      Time= Time+dt  
    end
    
  
    %%% DISTURBANCE
    
    if Time>=12   % initial time for plant development before disturbance - we let pine establish
    D=randi(10,1);%%% !!!! CHANGE THIS TO MAKE IT MORE INTUITIVE AND REALISTIC!!! Now we let pine reproduce once without disturbance
    end
    % if Time/10 is an integer there is a probability of 1/10 of fire every year and this does not depend from previous events
    if D == 1 
        'fire'
        Lit(:,:)=0;
        for i=1:m
            for j=1:m
                TC(i,j)= TC(i,j)*AR(TC(i,j)+1);  
                Age(i,j)= Age(i,j)*AR(TC(i,j)+1); 
                ReleaseSeeds= sum(SBPC); %the production of seeds when there is a fire is the total of the canopy seeds produced until that moment
                SBP1=0;SBP2=0;
                
    % !! check!! the canopy seed bank should be released when there is
    % disturbance - fire or tree harvesting
            end
        end
  
                
end
    
    imagesc(TC)
    set(h,'Clim',[-0.5 3.5]);
    colormap(VegetationColormap);
    colorbar
 
    drawnow;%pause
    
    %%% update abundance of different species in the lattice
    
    Pine=sum(sum(TC==1));
    Seeder=sum(sum(TC==2));
    Oak=sum(sum(TC==3));
    
    
    %%%%%%%%%%%%%%%% STORING AND VISUALIZATION %%%%%%%%%%%%%%%%%
    StoreTime = StoreTime - Time;
    if StoreTime <= 0
        StorePine(NrStore,:) = [Pine];
        StoreSeeder(NrStore,:) = [Seeder];
        StoreOak(NrStore,:) = [Oak];
        VectorTime(NrStore,:)= [Time];
        NrStore = NrStore+1;
        StoreTime = StoreStep;
    end %if StoreTime <= 0  
end
    
%%%Plotting over time
figure
plot(VectorTime,StorePine/m/m*100,VectorTime,StoreSeeder/m/m*100,VectorTime,StoreOak/m/m*100)
legend('Pine','Seeder','Resprouter')
set(gca,'fontsize',14);
set(gcf,'Position',[560         582        800         366],'PaperPositionMode','auto');
% saveas(gcf,'figureTime.png','png')

% Creates movie - not working yet
%     imagesc(TC);
%     movie(Frame)=getframe;
%     Frame=Frame+1;
%writerObj = VideoWriter(imageTC) %constructs a VideoWriter object to write video data to an AVI file with Motion JPEG compression.
    