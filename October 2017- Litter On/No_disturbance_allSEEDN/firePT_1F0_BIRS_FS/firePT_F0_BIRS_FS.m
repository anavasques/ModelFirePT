%Model PT fire
%CASCADE project
%Asynchronous CA
%Vasques et al. 2014-2016 (model development)
% ------------------------------------------------------------------------

close all
clear all

%VARIABLES
%-------------------------------------------------------------------------
%PER SPECIES

%%% PINE
AgeMP=10;                 % Age of maturity pine [year]
%SeedFP=1000;             % Number of seeds per cone (63)* cone per tree(15) Vega et al 2008 [n/m2/year]
SeedFP=100;               % Reduced in 10 times as in Quercus [n/m2/year]
LSP= 100;                 % Life span of pine % in "practices centro pinus" [year]
canopyBank=0.5;           % Percent of the seeds that are stored in the canopy (indicative)

%%% SEEDER
AgeMS=3;                  % Age of maturity seeder % field obs Calluna 1 [year]; Cistus 3 years [year]
SeedFS=1000;              % Seed production general value [n/m2/year]
LSS=30;                   % Life span Calluna % in woodland education centre [year]

%%% OAK
%AgeMO=50;                % Age of maturity [year] %Kew % !in Pausas 1999 maturity = 15!
AgeMO=20;                 % Age of maturity [year] % According to Ramon an oak can produce acorns after 15-20 years
SeedFQ=10;                % Seed production per occupied cell [n/m2/year]- 120 acorns per tree refered in Martin et al. 2014 [n/m2/year]
BirdSeedN=2;              % Annual seed input by birds in the whole laticce [n/year] % average values Q. suber Pons and Pausas 2007 50 seeds per hectar - this value depends on surrounding populations
RespAge=1;                % Age at which the oaks can resprout [year]

%LSO= 1000;               % Life span [year] Quercus robur % in forestar;
LSO=500;                  % Life span [year] it was reduced by half to be more realistic (also in the calculations of mortality)

%LITTER RELATED
minG= [0 0 0.3];          % Minimum probability of germination/establishment to generate the litter function(first pine second seeder third oak)
maxG= [0.9 0.9 0.9];      % Maximum probability of germination/establishment to generate the litter function(first pine second seeder third oak)
ProbPZeroL=0.7;           % Germination probability for pine when litter=0 cm
LitThreshP=3;             % Litter threshold for Pine above which there is no germination (cm)
LitThreshS=2;             % Litter threshold for seeders above which there is no germination (cm)
lrate=0.8;                % These values were approximated to have a curve close to the sigmoid- indication from literature: rate of litter deposition [cm/year] Fernandes et al 2004 have 0.42; Indication from Ramon: after 20-30 litter stabilizes
eflit=1-0.4;              % Effective litter (e.g. if 0.90 then 0.10 of the total litter is decomposed yearly - estimation that is not based in literature)
amp=[0.3 0.3 0];          % Amplitude of curve for the function litter vs. establishment(first pine second seeder third oak)


LitOn=1;                  % Switch for the experiments with litter (on=1; off=0)
lconv=0.5*ones(3,3);      % Convolution matrix for the deposition of litter around pine (distant neighbours)
lconv(2,2)=1;             % Convolution matrix for the deposition of litter around pine (closest neighbours)
z= 8;                     % Number of neighbours used in the convolution matrix


litLS=4;                  % [cm] This is the value that is used in the function of low severity

% GENERAL

est= [10 100 2];          % ProbS - max number of seedlings per cell CCD field from which we inferred a probability of establishment in one cell (type of carrying capacity of cell)
SeedLoss= [0 0.10 1];     % Rate of seed loss in the soil per year (1 pine 2 seeder 3 oak)
nrsp=length(SeedLoss);    % Number of species used in the model - variable that is used in the probability function
mort=1./[LSP,LSS,LSO];    % Yearly mortality (1/lifespan) (1 pine 2 seeder 3 oak)

% DISTURBANCE
AR= [0,0,1,1];            % Ability to resprout (0= not able; 1 = able) (first element is fake (bare soil); 2 pine, 3 seeder, 4 oak)

%D=0;                     % initialization of disturbance


% CONTROL VARIABLES (for the runs)
StartTime= 0;             % [year]
EndTime= 3000;            % [year]
StoreTime = 1;            % [year]

dt=1;                     % [year]
m= 100;                   % number of cells in lattice


%%%%% CODE FOR MULTIRUNS %%%%%%%

% here we can vary the values that will be used for the repetitions of the
% runs 

BirdSeedNv=2; 
%maxseedSeed=100:100000:1000000; % makes runs changing the parameter of SB (2) between the three values determined and keeping all other values fixed

%pd=3:10:43
b=142:1:151; % repetitions for the fire series (each number is a random number that corresponds to a fire sequence)

nruns=20;                % number of repetitions of each run (exactly the same conditions)


%%%% VECTOR FOR MULTIPLE RUNS WITH DIFFEREENT NUMBER OF SEEDS 
% oak
for k=BirdSeedNv
    BirdSeedN=BirdSeedNv;
% initial pine density
% for k=1:length(pd)
%     TC(pd:pd:m-pd,pd:pd:m-pd)=1;
%     SB=[0 100*m*m 0+randi(BirdSeedN,1)];
% seeds of seeder
% for k=1:length(maxseedSeed)
%     SB=[0 maxseedSeed(k) 0+randi(BirdSeedN,1)];
%     TC(3:3:m-3,3:3:m-3)=1;

%%%VECTOR OF FIRE OCCURRENCE FOR DIFFERENT FIRE SERIES - b

for b=122:141;
D=0*[StartTime:dt:EndTime];%#ok<NBRAK>
tf=4000;                   % time without fires
fireret=7;                % average interval between fires (fires are still stochastic)
rng(b); % INITIALISE THE RAND COMMANDS TO OBTAIN A SPECIF FIRE SEQUENCE         

% Multiruns
      for irun=1:nruns;
        k,irun

        %%%------------------------------------------------------------------------
        %%%INICIALIZATION OF THE MATRICES
        
        %TC - Type of Cover: 1- pine; 2- seeder; 3- oak; Age- plant age; SB- Seed Bank); Lit- Litter
        % -------------------------------------------------------------------------
        TC= zeros(m,m);           % Creates a matrix of size m*m filled with zeros
        Age= zeros(m,m);          % Creates a matrix of size m*m filled with zeros
        Lit= zeros(m+2,m+2);      % Creates a matrix of size (m+2)*(m+2) filled with zeros - the borders should be excluded in stats done with litter
        PosQSeed=zeros(m,m);      % Number of acorns per cell
        Seeder=0;                 % Initializes the number of cells with seeder
        Oak=0;                    % Initializes the number of cells with oak
        Pine=0;                   % Initializes the number of cells with pine
        
        %PLANT PINES aND INITIALIZING THE SEED BANK
        TC(3:3:m-3,3:3:m-3)=1;
        SB=[0 100*m*m BirdSeedN]
        SBP1=0;                   % [number of seeds] SBP1 and SBP2 are the only ways of initializing the seed bank every year
        SBP2=0;
        SBPC=0;                   % Initialization seed bank pine canopy
       
        %FUNCTIONS
        
        ProbL=[0 0 0];            % Probability of germination due to litter (first pine second seeder third oak)
        Litter=0;                 % Initializes the sum of the number of cells with litter
        indAge=0;                 % Parameter used in litter deposition
        ProbG=[0 0 0];            % Probability of germination first pine second seeder third oak
        ProbS= [0 0 0];           % To calculate probability based on seed production
        
        %STORING AND PLOTTING
        Time = StartTime;
        NrStore = 1;
        StoreStep = StoreTime;
        StorePine=zeros(EndTime,1);
        StoreSeeder=zeros(EndTime,1);
        StoreOak=zeros(EndTime,1);
        StoreLitter=zeros(EndTime,1);
        VectorTime=zeros(EndTime,1);
        MatPine=0;
        StoreMatPine=zeros(EndTime,1);
          
   %%%%%%%%%%%%%%%%%%%% DYNAMIC LOOP FOR FIRES

while tf<EndTime% EndTime can be substituted for the time when disturbance should stop
   
    tf=tf-fireret*log(rand(1,1)); % function stochastic fire recurrence mentioned in Baudena et al 2010
    D(round(tf))=1;
  
end %while loop tf
        
        %--------------------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%% MAIN DYNAMIC LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %--------------------------------------------------------------------------
        
        while Time < EndTime
            Time= Time+dt;
            %%%% Creates LITTER in the neighborhod of pine (8 neighbors)+ pine plce itself
            [x,y]=find(TC==1); %finds cells with pine in the whole matrix
            x=x+1;y=y+1;
            for i=1:length(x)
                %%%%%% Closer to a sigmoid curve
                indAge=Age(x(i)-1,y(i)-1)>AgeMP;
                Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)=Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)+lrate*lconv*dt*(indAge+(1-indAge)*.2); % .2 IS THE 20% OF THE MAXIMUM VALUE OF LITTER DEPOSITION RATE FOR PINE<AGEMP
            end %litter deposition
            
            % SEED BANK CALCULATION (ONCE A YEAR)
            
            %PINE soil seeds (SB1) and canopy seeds (SBPC)
            SB(1)=SBP1+SBP2+SeedFP*(1-canopyBank)*sum(sum(TC(Age>AgeMP)==1));% Pine TWO YEARS OF SEED LIFE; 1-canopybank is doing the same as canopy bank, i.e. *0.5
            SBPC=SBPC+SeedFP*canopyBank*sum(TC(Age>AgeMP)==1);
            %SEEDER
            SB(2)=SB(2)+SeedFS*(sum(sum(TC==2)))-SeedLoss(2)*SB(2);           % LONG SEED LIFE
            %OAK
            SB(3)=SeedFQ*(sum(sum(TC(Age>AgeMO)==3)))+BirdSeedN;
            
            %%%%%% Small loop only to spread the acorns over the lattice
            for kk=1:SB(3)      %only happens if SB3 is bigger than 1
                cc=randi(m*m);  % MARA May17 2016: I CHANGED THIS TO GENERATE ONE RAND NUMBER ONLY BTW 1 AND m*m AND USING THE FACT THAT MATLAB CAN READ MATRICES AS LONG VECTOR (SO YOU CAN GIVE ONE COORDINATE ONLY FOR EACH CELL)
                    PosQSeed(cc)=PosQSeed(cc)+1;
            end %end of the loop seeds oak
            
            %------------------------------------------------------------------
            % COLONIZATION OF AN EMPTY CELL AND MORTALITY FOR VEGETATED CELLS
            %------------------------------------------------------------------
            % At each time step tests who is going to colonize an empty cell in the lattice based on the
            % availiable seeds (seed bank) and on litter
            for i = 1 : m
                for j=1:m           % for all the cells of the lattice
                    test=rand*nrsp; % tests for a random number BETWEEN 0 AND THE NUMBER OF SPECIES (LENGTH(PROBg=3))
                    if TC(i,j)==0   % empty cells
                        
                        %%% TERM FOR PROBABILITY OF COLONIZATION vs. NUMBER
                        %%% OF SEEDS AND ESTABLISHMENT (Prob S)
                        
                        ProbS(1)=1-(1-1/est(1))^(SB(1)/m/m);
                        ProbS(2)=1-(1-1/est(2))^(SB(2)/m/m); % FOR PINE AND SEEDERS, SEEDS ARE EQUALLY SPREAD THROUGHOUT THE CELLS; this was taken in the paper: Cannas et al. 2003
                        % only for oak the seeds are spread over the lattice
                        ProbS(3)=1-(1-1/est(3))^PosQSeed(i,j); %if ProsQSeed=0 the whole term goes to zero
                        
                        %%%COLONIZATION vs LITTER (Prob L)
                        
                        ProbL(1)=(maxG(1)+minG(1))/2+(maxG(1)-minG(1))/2*tanh((LitThreshP-Lit(i+1,j+1))/amp(1)) ...
                            -(maxG(1)-ProbPZeroL)*exp(-2/LitThreshP*exp(1)*Lit(i+1,j+1)); % PINE
                        ProbL(2)=(maxG(2)+minG(2))/2+(maxG(2)-minG(2))/2*tanh((LitThreshS-Lit(i+1,j+1))/amp(2)); % SEEDER ampS=0.3 max=.9 min=0.
                        ProbL(3)=maxG(3)-(maxG(3)-minG(3))*exp(-Lit(i+1,j+1)); % QUERCUS
                        ProbL=ProbL*LitOn+(1-LitOn)*[1 1 1];% term for test IN ABSENCE of litter
                        
                        % COMBINING PROBABILITIES OF ESTABLISHMENT (Prob G)
                        % version modified from Alain's i.e. prob only needs to
                        % be <than 3 and we don't need to multiply it by dt anymore)
                        
                        ProbG=ProbL.*ProbS; %overall term for colonization test
                        
                        if test>ProbG(1)+ProbG(2)+ProbG(3)
                            TC(i,j)=0; %this was modified to make the code faster so if the cell is not colonized it stops running here
                        elseif test>ProbG(1)+ProbG(2)
                            TC(i,j)=3;
                            Age(i,j)=dt;
                        elseif test>ProbG(1)
                            TC(i,j)=2;
                            Age(i,j)=dt;
                            SB(2)=SB(2)-100; % for every cell that is colonized by pine there are 100 less shrub seeds in SB
                        else
                            TC(i,j)=1;
                            Age(i,j)=dt;
                            SB(1)=SB(1)-10; % for every cell that is colonized by pine there are 10 less pine seeds in SB
                        end
                    else
                        Age(i,j)=Age(i,j)+dt;
                        
                        if test< mort(TC(i,j))*dt %determines if a plant in a cell dies, mortality is defined according to life span
                            TC(i,j)=0;
                            Age(i,j)=0;
                        end
                    end % end of if empty cell
                end % end of lookinf for columns
            end % end of looking for rows
            
%%%% DISTURBANCE
 % NOTE - the rest of the disturbance term is at the beggining of the code
    D1=D(Time);
    if D1== 1
        'fire'; % can be written when a fire occurs

        %HIGH SEVERITY
        Lit(:,:)=0;
        
        %LOW SEVERITY
        %Lit(:,:)=litLS*sum(sum(TC(Age>AgeMP)==1))/m/m*8; %leaves from the canopy fall creating a litter
        %layer - for simplification purposes the litter in the soil is
        %maintained - it is porportional to the canopy cover anyway - canopy effect on 8
        %neightbours
        for i=1:m
            for j=1:m
                % Function that determines that a seeder dies and a
                % resprouter continues to occupy the cell
                % Kind IS A TRICK TO AVOID IF STRUCTURE (IF AGE>RespAge THEN RESPROUT)
                kind=floor(AR(TC(i,j)+1)*Age(i,j)/RespAge);
                kind=round(kind/(kind+eps));
                TC(i,j)= TC(i,j)*kind;
                Age(i,j)= Age(i,j)*kind;
            end
        end
        
        % IN ANY CASE (HS or LS) PINE CANOPY SEEDS FALL INTO SEEDBANK
        SB(1)= SBPC; %the production of seeds when there is a fire is the total of the canopy seeds produced until that moment
        SBP1=0;SBPC=0; %SBP2=0; Not needed to put sbp2 to zero
    else
        Lit=eflit*Lit; % effective litter, i.e. litter that is not degraded and remain for the years after - should be updated here to also occur in the empty cells
    end %end of fire occurence D=1
            
        % RE- INITIALIZATION OF SEED BANKS
        
        % UPDATES OF THE PINE SEEDBANK
        SBP2=0.5*SBP1; % PINE SEED BANK OF TWO YEARS BEFORE IS 50%
        SBP1=SB(1);% PINE SEED BANK OF 1 YEAR BEFORE
        % RESETS NUMBER OF QUERCUS SEEDS TO ZERO FOR NEXT YEAR, REINITIALISE
        % EVERY YEAR THE NUMBER OF SEEDS COMING IN
        PosQSeed=0*PosQSeed;      % NUMBER OF QUERCUS SEEDS PER CELL
        SB(3)=BirdSeedN; %This is the term to get a new random number between 1-5 every year
            
        % STORES THE VARIABLES FOR PLOTTING
        Pine=sum(sum(TC==1));
        MatPine=sum(sum(TC==1&Age>AgeMP));
        Seeder=sum(sum(TC==2));
        Oak=sum(sum(TC==3));
        Litter=sum(sum(Lit(2:m+1,2:m+1)));
        %Litter=max(max(Lit(2:m+1,2:m+1)));%
        %AgeMTX=mean(mean(Age));
        StorePine(NrStore) = Pine; % NOTICE THESE VECTORS ARE AS LONG AS ENDTIMES, and as wide as 1 (vector not matrices)
%       StoreMatPine(NrStore)= MatPine;
        StoreSeeder(NrStore) = Seeder;
        StoreOak(NrStore) = Oak;
        StoreLitter(NrStore)= Litter;
        %StoreAge(NrStore)=AgeMTX; %correct this
        VectorTime(NrStore)= Time;
        NrStore = NrStore+1;
        StoreTime = StoreStep;
        
        %%% Plots final figure of vegetation (also the pattern)
        % figure (pattern)
        % h=subplot(1,1,1);
        % imagesc(TC)
        % set(h,'Clim',[-0.5 3.5]);
        % set(gca,'FontSize',20,'fontWeight','bold')
        % colormap(VegetationColormap);
        % colorbar
        % drawnow; %pause
        %%%%% figure over time
        % figure
        % plot(VectorTime,mPineMR/m/m*100,'b', VectorTime,mSeederMR/m/m*100, 'r--.', VectorTime,mOakMR/m/m*100, 'g*')%, VectorTime,StoreLitter/m/m, 'k.')%, VectorTime,StoreAge,'gr')
        %%%% figure for Litter depth over time
        % figure
        % set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
        % plot(VectorTime,StoreLitter/m/m)
        % xlabel('Time (year)');
        % ylabel ('Mean litter depth (cm)'); 
            
        end % end of the main dynamic loop
        
        %%%%%% MULTIRUNS CODE
        filename=strcat(['firePT_F0_BIRS_FS',num2str(BirdSeedNv),'_',num2str(irun),'_',num2str(b),'.mat' ]);
        save(filename,'StorePine','StoreSeeder','StoreOak','StoreLitter','VectorTime')

    end % stochastic runs
    
%     %%%%%% FIGURES OF THE MULTIRUNS(INCLUDED EIN SEPARATE FILES)
%     figure
%     plot(VectorTime,mPineMR/m/m*100,'b', VectorTime,mSeederMR/m/m*100, 'r--.', VectorTime,mOakMR/m/m*100, 'g*')%, VectorTime,StoreLitter/m/m, 'k.')%, VectorTime,StoreAge,'gr')

end % of different fire series
end %of multiruns number of seeds
    
