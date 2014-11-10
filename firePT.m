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
LSO= 1000;                % Life span quercus robur % in forestar

%GENERAL
minG= [0 0 0.3];          % minumum germination first pine second seeder third oak
maxG= [0.1 0.9 0.9];      % maximum germination first pine second seeder third oak
SeedLoss= [0.50 0.05 1];  % rate seed loss first pine second seeder third oak

lrate=0.42;               % rate of litter deposition [cm/year] Fernandes et al 2004
ProbG=[0 0 0];            % probability of germination first pine second seeder third oak
amp=[0 0.3 0];            % amplitude of curve interaction with litter

SBP1=0;                   % !SBP1 and SBP2 are only ways of initializing the seed bank every year
SBP2=0;
SBPC=0;                   % Inicialization seed bank pine canopy
SB= [0 0 0];              % seed bank first pine second seeder third oak
ProbS= [0 0 0];           % to calculate probability based on seed prod
est= [7 400];             % max number of seedlings per cell CCD field

mort=1./[LSP,LSS,LSO];    % MORTALITY of pine, seeder, oak = 1/lifespan




AR= [0,0,0,1];            % ability to resprout: first element is fake (bare soil); pine=0, seeder=0, oak=1;

% Control constants
StartTime= 0;             % [year]
EndTime= 200;             % [year]
PlotStep = 1;             % [year]
StoreStep = 1;            % [year]

% NOTE: put dt smaller than one year in a way that the probabilities are <1 but not too small otherwise the model runs slowly
dt= 1;                    % [year]
m= 100;                   % for size of lattice [meter]

% Control Variables
NrStore = 1;
PlotTime = PlotStep;
StoreTime = StoreStep;
Time = StartTime;

D= 0;                     % initialization only
Pine=0;                   % will count the number of cells with pine
Seeder=0;                 % will count the number of cells with pine
Oak=0;                    % will count the number of cells with pine

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
SB=[100 1000 0+randi(BirdSeedN,1)];

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
        end
        % Colonization of an empty cell & mortality for vegetated cells
        %------------------------------------------------------------------
        % At each time step tests who is going to colonize an empty cell in the lattice based on the
        % existent seed bank
        
        for i = 1 : m
            for j=1:m
                test=rand;
                if TC(i,j)==0 % colonization/germination
                    % LITTER DEPENDENCE:
                    ProbG(1)=(maxG(1)+minG(1))/2;%(maxG(1)+minG(1))/2+(maxG(1)-minG(1))/2*tanh((Lit(i,j)-2)/amp(1)); % PINE
                    ProbG(2)=(maxG(2)+minG(2))/2+(maxG(2)-minG(2))/2*tanh((2-Lit(i,j))/amp(2)); % SEEDER ampS=0.3 max=.9 min=0.
                    ProbG(3)=maxG(3)-(maxG(3)-minG(3))*exp(-Lit(i,j)); % QUERCUS
                    
                    ProbG=ProbG*dt; %this is the trick to get probability small
                    
                    %%% Term for the relation between available seeds and Prob S.                      
                    ProbS(1:2)=1-(1-1./est(1:2)).^(SB(1:2)/m/m); % FOR PINE AND SEEDERS, SEEDS ARE EQUALLY PSREAD THROUGHOUT THE CELLS
                    ProbS(3)=1;                                  % Check this again; it was zero that's why there were never quercus
                    for ii=1:length(coordseed)
                        ProbS(3)=ProbS(3)+(coordseed(ii,1)==i&coordseed(ii,2)==j);
                    end
                    ProbS(3)=ProbS(3)>1;
                    
                    ProbG=ProbG.*ProbS;
                    %this step has the reference  of Alains' MSc thesis ->
                    %CAN BE IMPROVED
                    if sum(ProbG)*dt>1
                        'sum of probability higher than 1! Please decrease dt'
                        break
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
                    end
                end
            end
        end
      Time= Time+dt  
    end
    
    % SEED BANK CALCULATION ONLY ONCE A YEAR
    SB(1)=SBP1+SBP2+SeedFP*(1-canopyBank)*sum(sum(TC(Age>AgeMP)==1)); % TWO YEARS OF SEED LIFE
    SB(1)=SB(1)-SeedLoss(1)*SB(1);
    SB(2)=SB(2)+SeedFS*(sum(sum(TC==2)))-SeedLoss(2)*SB(2);           % LONG SEED LIFE
    SB(3)=SeedFQ*(sum(sum(TC(Age>AgeMO)==3)))+randi(BirdSeedN,1);     % NO MEMORY
    SB(3)=SB(3)-SeedLoss(3)*SB(3);
    if SB(3)>0
        coordseed=randi(m,SB(3),2);
    end
    
    SBP2=SBP1; % PINE SEED BANK OF TWO YEARS BEFORE
    SBP1=SB(1);% PINE SEED BANK OF 1 YEAR BEFORE
   
    % accumulation of seeds in the canopy
    SBPC=SeedFP*canopyBank*sum(TC(Age>AgeMP)==1); %canoyBank=% of the seeds that stay in the
    % canopy and accumulate over time
    
    %%% DISTURBANCE
    
    D=randi(10,1);%%% !!!! CHANGE THIS TO MAKE IT MORE INTUITIVE AND REALISTIC!!!
    % if Time/10 is an integer there is a probability of 1/10 of fire every year and this does not depend from previous events
    if D == 1 
        'fire'
        Lit(:,:)=0;
        for i=1:m
            for j=1:m
                TC(i,j)= TC(i,j)*AR(TC(i,j)+1);  
                Age(i,j)= Age(i,j)*AR(TC(i,j)+1); 
                SB(1)=SBPC;
                SBP1=0;SBP2=0;
            end
        end
    end
    
    %%% update abundance of different species in the lattice
    
    Pine=sum(sum(TC==1));
    Seeder=sum(sum(TC==2));
    Oak=sum(sum(TC==3));
    
    imagesc(TC)
    set(h,'Clim',[-0.5 3.5]);
    colormap(VegetationColormap);
    colorbar
 
    drawnow;%pause
    
    %%%%%%%%%%%%%%%% STORING AND VISUALIZATION %%%%%%%%%%%%%%%%%
    StoreTime = StoreTime - Time;
    if StoreTime <= 0
    StorePine(NrStore,:) = [Time Pine]; 
    StoreSeeder(NrStore,:) = [Time Seeder];
    StoreOak(NrStore,:) = [Time Oak]; 
    VectorTime(NrStore,:)= [Time];
    NrStore = NrStore+1;
    StoreTime = StoreStep;
    end %if StoreTime <= 0
  
    
end
    %%% Improve this part to get the final plot working
    PlotTime = PlotTime-dt;
    if PlotTime <= 0
  
    x = VectorTime;
    y1 = StorePine;
    y2 = StoreSeeder;
    y3= StoreOak;
    figure
    plot(x,y1,x,y2,x,y3)
    legend('Pine','Seeder','Resprouter');

    PlotTime = PlotStep;
    end
  
%     Creates movie
%     showimagesc(TC);
%     movie(Frame)=getframe;
%     Frame=Frame+1;

                
                