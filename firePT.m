%Model PT fire
%CASCADE project
%Asynchronous CA
%Vasques et al. 2014-2015 (model development)
%PROBABILITY MARA'S - MODIFIED FROM ALAIN'S
% ------------------------------------------------------------------------

close all
clear all

%Parameter values
%-------------------------------------------------------------------------
%PER SPECIES

%%% PINE
AgeMP=10;                 % Age of maturity pine %start at 6 and regularly 10-15 % in Cronk and Fuller, 1995 [year]
SeedFP=945;               % Seed production per pine mature tree number of seeds per cone (63)* cone per tree (15) Vega et al 2008
LSP= 100;                 % Life span of pine % in "practices centro pinus" [year]
canopyBank=0.5;           % Percent of the seeds that are stored in the canopy maybe reduce it too many pines
ReleaseSeeds=0;           % Pine seeds in the canopy that are released after the fire
SBP1=0;                   % !SBP1 and SBP2 are only ways of initializing the seed bank every year
SBP2=0;
SBPC=0;                   % Initialization seed bank pine canopy
Pine=0;                   % will count the number of cells with pine

%%% SEEDER
AgeMS=3;                  % Age of maturity seeder % field obs Calluna 1 [year]; Cistus 3 years ref
%SeedFS=400;
SeedFS=1000;              % Seed production per plant/occupied cell approx value ADJUST
LSS=30;                   % Life span calluna % in woodland education centre [year]
Seeder=0;                 % will count the number of cells with seeder

%%% OAK
AgeMO=50;                 % Age of maturity % Kew % [year]% !Pausas 1999 has maturity = 15!
SeedFQ=12;                % Seed production oak per occupied cell - 120 acorns per tree refered in Martin?k et al. 2014% [n/m2/year]
BirdSeedN=5;              % Annual seed input by birds - based on average values Q. suber Pons and Pausas 2007 - this value depends on surrounding populations
RespAge=10;                % ONLY OAKS OLDER THAN THIS AGE CAN RESPROUT

%LSO= 1000;               % Life span quercus robur % in forestar;
LSO=500;                  % this was reduced by half to be more realistic (also in the calculations of mortality)
Oak=0;                    % will count the number of cells with oak

%LITTER RELATED
minG= [0 0 0.3];          % minimum germination first pine second seeder third oak
maxG= [0.9 0.9 0.9];      % maximum germination first pine second seeder third oak
ProbPZeroL=0.7;           % Germination probability for pine when litter=0 cm
LitThreshP=3;             % Litter threshold for Pine above which ~no germination (cm)
LitThreshS=2;             % Litter threshold for seeders above which ~no germination (cm)
lrate=0.42;               % rate of litter deposition [cm/year] Fernandes et al 2004
eflit=0.90;               % effective litter: if 0.90 then 0.10 of the total litter is decomposed - estimated value not from literature
ProbL=[0 0 0];            % probability of germination due to litter (first pine second seeder third oak)
amp=[0.3 0.3 0];          % amplitude of curve interaction with litter
Litter=0;                 % to sum the sumber of cells with litter

% GENERAL
ProbG=[0 0 0];            % probability of germination first pine second seeder third oak
ProbS= [0 0 0];           % to calculate probability based on seed prod
est= [7 100 1.5];         % max number of seedlings per cell CCD field from which we inferred a probability of establishment in one cell
%est= [7 400 7];          % infered probability of establishment for ProbS

SeedLoss= [0 0.10 1];     % rate seed loss soil seed bank 1 pine 2 seeder 3 oak
%nrsp=3                   % number of species used in the model - to put in the prob expression
mort=1./[LSP,LSS,LSO];    % MORTALITY of pine, seeder, oak = 1/lifespan
AR= [0,0,0,1];            % Ability to resprout: first element is fake (bare soil); pine=0, seeder=0, oak=1;

D=0;                      % initialization of disturbance

% CONTROL CONSTANTS AND VARIABLES
StartTime= 0;             % [year]
EndTime= 500;             % [year]
StoreTime = 1;            % [year]

dt=1;                     % [year]
m= 100;                   % for size of lattice [meter]

Time = StartTime;
NrStore = 1;
StoreStep = StoreTime;
StorePine=zeros(EndTime,1);
StoreSeeder=zeros(EndTime,1);
StoreOak=zeros(EndTime,1);
StoreLitter=zeros(EndTime,1);
VectorTime=zeros(EndTime,1);

% FOR MULTIRUNS
nruns=10;

%%%INICIALIZATION OF THE MATRICES

%TC - Type of Cover: 1- pine; 2- seeder; 3- oak; Age- plant age; SB- Seed Bank); Lit- Litter
% -------------------------------------------------------------------------
TC= zeros(m,m);           % Creates a matrix of size m*m filled with zeros
Age= zeros(m,m);          % Creates a matrix of size m*m filled with zeros
Lit= zeros(m,m);          % Creates a matrix of size m*m filled with zeros
z= 8;                     % Number of neighbours

TC(4:4:m-4,4:4:m-4)= 1;   % plants 1 pine every 4 meters - dense prodution stand excluding the borders
%TC(40:40:m-40,40:40:m-40)= 1; % plants 1 pine every 40 meters - dense prodution stand excluding the borders

PosQSeed=zeros(m,m);      % NUMBER OF QUERCUS SEEDS PER CELL

SB=[0 1*m*m 0+randi(BirdSeedN,1)]; %initial seed bank
%SB=[0 1000 0+randi(BirdSeedN,1)]; %previous number of seeds changing initial conditions for seeder and oak, pine is planted but can also be seeded randomly

%%%VECTOR OF FIRE OCCURRENCE
D=0*[StartTime:dt:EndTime];%#ok<NBRAK>
tf=10;                     %time without fires
fireret=5;                 %intervale between fires - fire return
rand('state',120)
while tf<EndTime
    tf=tf-fireret*log(rand(1,1)); %stochastic fire recurrence Baudena et al 2010
    D(round(tf))=1;
end

%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%DYNAMIC LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

while Time < EndTime
    Time= Time+dt

    %%%Creates LITTER in the neighb of pine (8 neighbors)+ pine site
    
    %if TC(Age>AgeMP)==10 % pine does not create litter in the
    %first years
        [x,y]=find(TC(2:end-1,2:end-1)==1); %consider including here litter deposition only after a certain age
        x=x+1;y=y+1;
        for i=1:length(x)
        Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)=Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)+lrate*dt;
        end
 
    % SEED BANK CALCULATION (ONCE A YEAR)
    %%%%%%%%%%%%%%Check this with Mara%%%%%%%%%
    %%%%%%%%%%%%%%this was why the seeds of pine did not come back%
    SB(1)=SB(1)+ SBP1+SBP2+SeedFP*(1-canopyBank)*sum(sum(TC(Age>AgeMP)==1));% TWO YEARS OF SEED LIFE; 1-canopybank is doing the same as canopy bank, i.e. *0.5
    %SB(1)=SB(1)-SeedLoss(1)*SB(1); % this ter is not needed as pine seeds only last 2 years and then die
    SB(2)=SB(2)+SeedFS*(sum(sum(TC==2)))-SeedLoss(2)*SB(2);           % LONG SEED LIFE
    SB(3)=SeedFQ*(sum(sum(TC(Age>AgeMO)==3)))+randi(BirdSeedN,1);     % NO MEMORY
    % puts the seeds that arrive in random coordinates of the lattice
    % OLD WAY
    %     if SB(3)>0
    %     coordseed=randi(m,SB(3),2);
    %     end
    % NEW WAY
    for kk=1:SB(3)
        c1=randi(m,1,1);c2=randi(m,1,1);
        PosQSeed(c1,c2)=PosQSeed(c1,c2)+1;
    end
    
    % accumulation of pine seeds in the canopy
    SBPC=SBPC+SeedFP*canopyBank*sum(TC(Age>AgeMP)==1); %canopyBank=% of the seeds that stay in the canopy and accumulate over time
    
    %------------------------------------------------------------------
    % COLONIZATION OF AN EMPTY CELL AND MORTALITY FOR VEGETATED CELLS
    %------------------------------------------------------------------
    % At each time step tests who is going to colonize an empty cell in the lattice based on the
    % availiable seeds (seed production and seed bank)
    for i = 1 : m
        for j=1:m
            test=rand*length(ProbL); %RANDOM NUMBER BETWEEN 0 AND THE NUMBER OF SPECIES (LENGTH(PROBg=3))
            if TC(i,j)==0 % colonization/germination
                
                
   %%% TERM FOR PROBABILITY OF COLONIZATION vs. NUMBER OF SEEDS AND ESTABLISHMENT
                
                ProbS(1)=1-(1-1/est(1))^(SB(1)/m/m); 
                ProbS(2)=1-(1-1/est(2))^(SB(2)/m/m); % FOR PINE AND SEEDERS, SEEDS ARE EQUALLY SPREAD THROUGHOUT THE CELLS; this was taken in the paper: Cannas et al. 2003
                % to oak
% %             mm=find(sum(ismember(coordseed(:,1:2),[i,j]),2)>=2);
% %             SAME AS LINE ABOVE BUT WITH A FASTER ALGOORITHM (FOUND ON THE
% %             INTERNET):
%               mm=find(sum(builtin('_ismemberoneoutput',coordseed(:,1:2),[i,j]),2)>=2);
%               lengthmm=length(mm);
                % EXPRESSION FOR PROB QUERCUS WITH A THRESHOLD FROM ONE SEED
                % UP, PROBS=0.8:
                %ProbS(3)=0;
                %ProbS(3)=ProbS(3)+round(lengthmm/(lengthmm+eps))*.8; %
                % EXPRESSION SIMILAR TO SEEDERS AND PINES: 
%                 ProbS(3)=1-(1-1./est(3)).^lengthmm;
                ProbS(3)=1-(1-1/est(3))^PosQSeed(i,j);
                
                %%%% VERSION 1 (discarded)
                %for ii=1:size(coordseed,1)
                %ProbS(3)=ProbS(3)+(coordseed(ii,1)==i&coordseed(ii,2)==j);
                %end
                %ProbS(3)=ProbS(3)>=1; % if there is one seed or more -> prob=1
                
% COLONIZATION vs LITTER
%               ProbL=[1 1 1];% term for test IN ABSENCE of litter
                ProbL(1)=(maxG(1)+minG(1))/2+(maxG(1)-minG(1))/2*tanh((LitThreshP-Lit(i,j))/amp(1)) ...
                    -(maxG(1)-ProbPZeroL)*exp(-2/LitThreshP*exp(1)*Lit(i,j)); % PINE
                ProbL(2)=(maxG(2)+minG(2))/2+(maxG(2)-minG(2))/2*tanh((LitThreshS-Lit(i,j))/amp(2)); % SEEDER ampS=0.3 max=.9 min=0.
                ProbL(3)=maxG(3)-(maxG(3)-minG(3))*exp(-Lit(i,j)); % QUERCUS
                
% COMBINING PROBABILITIES OF ESTABLISHMENT DUE TO SEED NUMBERS AND LITTER
                
                ProbG=ProbL.*ProbS; %term for test with litter
                
                % this step has the improved version of Alains' (done by
                % Mara i.e. prob only needs to be <than 3 and we don't multiply by dt anymore)
                
                if test>ProbG(1)+ProbG(2)+ProbG(3)
                    TC(i,j)=0;
                elseif test>ProbG(1)+ProbG(2)
                    TC(i,j)=3;
                    Age(i,j)=dt;
                elseif test>ProbG(1)
                    TC(i,j)=2;
                    Age(i,j)=dt;
                else
                    TC(i,j)=1;
                    Age(i,j)=dt;
                end
            else
                Age(i,j)=Age(i,j)+dt;
                Lit(i,j)=eflit*Lit(i,j); % effective litter, i.e. litter that is not degraded and remain for the years after
                
                if test< mort(TC(i,j))*dt
                    TC(i,j)=0;
                    Age(i,j)=0;
                end
            end
        end
        SBP2=SBP1; % PINE SEED BANK OF TWO YEARS BEFORE
        SBP1=SB(1);% PINE SEED BANK OF 1 YEAR BEFORE
    end
    
    %%%% DISTURBANCE
    %%%D=randi(10,1)*(Time>=12); %this term is no longer needed
    D1=D(Time);
    if D1== 1
        'fire';
        Lit(:,:)=0;
        for i=1:m
            for j=1:m
                % PINES AND SEEDERS DIE; QUERCUS RESPROUTS IF OLDER OR
                % EQUAL THEN Respage
                % Kind IS A TRICK TO AVOID IF STRUCTURE (IF AGE>RespAge THEN RESPROUT)
                kind=floor(AR(TC(i,j)+1)*Age(i,j)/RespAge);
                kind=round(kind/(kind+eps));
                TC(i,j)= TC(i,j)*kind;
                Age(i,j)= Age(i,j)*kind;
                
                % PINE CANOPY SEEDS FALL INTO SEEDBANK 
                SB(1)= SB(1)+sum(SBPC); %the production of seeds when there is a fire is the total of the canopy seeds produced until that moment
                SBP1=0;SBP2=0;SBPC=0;
            end
        end
    end
    
    %%% updates abundance of different species in the lattice
    
    Pine=sum(sum(TC==1));
    Seeder=sum(sum(TC==2));
    Oak=sum(sum(TC==3));
    Litter=sum(sum(Lit>0.1));
    
    % RESET NUMBER OF QUERCUS SEEDS TO ZERO FOR NEXT YEAR
    PosQSeed=zeros(m,m);      % NUMBER OF QUERCUS SEEDS PER CELL
    
    %%%%%%%%%%%%%%%% STORING AND VISUALIZATION %%%%%%%%%%%%%%%%%
    %     StoreTime = StoreTime - Time; % (Mara) I COMMENTED THIS BECAUSE YOU WANT TO
    %     STORE EVERY TIME STEP SO IT' NOT USEFUL TO HAVE THIS EXTRA IF. TO BE
    %     RESTORED (AND CHANGED DIMENSIONS OF THE VECTORS BELOW) IF DT<1 YEAR
    %     OR YOU WANT TO SAVE EVERY E.G. 10 YEARS
    %     if StoreTime <= 0
    StorePine(NrStore) = Pine; % NOTICE THESE VECTORS ARE AS LONG AS ENDTIMES, and as wide as 1 (vector not matrices) NOT M BY M AS YOU DEFINED THEM.. -> REDIFINING ABOVE SHOULD TAKE LESS TIME -> LET ME KNOW!
    StoreSeeder(NrStore) = Seeder;
    StoreOak(NrStore) = Oak;
    StoreLitter(NrStore)= Litter;
    VectorTime(NrStore)= Time;
    NrStore = NrStore+1;
    StoreTime = StoreStep;
    % end %if StoreTime <= 0
    
    %StoreSpecies=[StorePine StoreSeeder StoreOak];
    % xlswrite('Sp abundance pine,seeder,oak',StoreSpecies)
    end
    
 %%% Plots final figure

figure
white=[1 1 1];
blue=[0 0 1];
red=[1 0 0];
green=[0 1 0];
VegetationColormap=[white; blue;red;green];
h=subplot(1,1,1);
imagesc(TC)
set(h,'Clim',[-0.5 3.5]);
colormap(VegetationColormap);
colorbar
drawnow;pause

%%%Plotting over time
figure
plot(VectorTime,StorePine/m/m*100,'b', VectorTime,StoreSeeder/m/m*100, 'r--.', VectorTime,StoreOak/m/m*100, 'g*', VectorTime,StoreLitter/m/m*100, 'kx')
legend('Pine','Seeder','Oak', 'Litter')
set(gca,'fontsize',14);
set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
xlabel('Time (year)');
ylabel ('Cover (%)');
drawnow;pause
% saveas(gcf,'figureTime.png','png')

% Creates movie - not working yet
%     imagesc(TC);
%     movie(Frame)=getframe;
%     Frame=Frame+1;
%     writerObj = VideoWriter(imageTC) %constructs a VideoWriter object to write video data to an AVI file with Motion JPEG compression.
