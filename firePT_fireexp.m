%Model PT fire
%CASCADE project
%Asynchronous CA
%Vasques et al. 2014 (model development)
%PROBABILITY MARA'S - MODIFIED FROM ALAIN'S
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
canopyBank=0.5;           % Percent of the seeds that are stored in the canopy maybe reduce it too many pines
ReleaseSeeds=0;           % Pine seeds in the canopy that are released after the fire

%%% SEEDER
AgeMS=2;                  % Age of maturity seeder % field obs Calluna 1 [year]
% Cistus 3 years ref
SeedFS=400;               % Seed production per plant/occupied cell approx value ADJUST
% check in lit
LSS=30;                   % Life span calluna % in woodland education centre [year]

%%% OAK
AgeMO=50;                 % Age of maturity seeder % Kew % [year]% !!! Pausas 1999 has maturity = 15!!!

SeedFQ=12;                % Seed production oak per occupied cell - 120 acorns per tree refered in Martin?k et al. 2014% [n/m2/year]
BirdSeedN=5;              % Annual seed input by birds - based on average values Q. suber Pons and Pausas 2007 - this value depends on surrounding populations

LSO= 1000;                % Life span quercus robur % in forestar

% FIRE
fireret=5;                  %year

% GENERAL
minG= [0 0 0.3];          % minimum germination first pine second seeder third oak
maxG= [0.9 0.9 0.9];      % maximum germination first pine second seeder third oak

ProbPZeroL=0.7;           % Germination probability for pine when litter=0 cm

LitThreshP=3;             % Litter threshold for Pine above which ~no germination (cm)
LitThreshS=2;             % Litter threshold for seeders above which ~no germination (cm)


lrate=0.42;               % rate of litter deposition [cm/year] Fernandes et al 2004
eflit=0.90;               % effective litter: if 0.90 then 0.10 of the total litter is decomposed - estimated value not from literature
ProbG=[0 0 0];            % probability of germination first pine second seeder third oak
amp=[0.3 0.3 0];          % amplitude of curve interaction with litter

SBP1=0;                   % !SBP1 and SBP2 are only ways of initializing the seed bank every year
SBP2=0;
SBPC=0;                   % Initialization seed bank pine canopy
ProbS= [0 0 0];           % to calculate probability based on seed prod
est= [7 400];             % max number of seedlings per cell CCD field from which we inferred a probability of establishment in one cell
SeedLoss= [0 0.10 1];     % rate seed loss first pine second seeder third oak
%nrsp=3                   % number of species used in the model - to put in the prob expression

mort=1./[LSP,LSS,LSO];    % MORTALITY of pine, seeder, oak = 1/lifespan

AR= [0,0,0,1];            % ability to resprout: first element is fake (bare soil); pine=0, seeder=0, oak=1;

% Control constants
StartTime= 0;             % [year]
EndTime= 500;             % [year]
StoreTime = 1;            % [year]

dt=1;                    % [year] 
m= 100;                   % for size of lattice [meter]

% Control Variables
NrStore = 1;
StoreStep = StoreTime;
Time = StartTime;

D=0;                      % initialization only
Pine=0;                   % will count the number of cells with pine
Seeder=0;                 % will count the number of cells with pine
Oak=0;                    % will count the number of cells with pine
StorePine=zeros(EndTime,1);
StoreSeeder=zeros(EndTime,1);
StoreOak=zeros(EndTime,1);
VectorTime=zeros(EndTime,1);
nruns=10;


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
SB=[0 1000 0+randi(BirdSeedN,1)]; %changing initial conditions for seeder and oak, pine is planted but can also be seeded randomly

% VECTOR OF FIRE OCCURRENCE
D=0*[StartTime:dt:EndTime]; %#ok<NBRAK>
tf=12;
while tf<EndTime
    tf=tf-fireret*log(rand(1,1));
    D(round(tf))=1;
end

% %Creates colormap
% figure
% white=[1 1 1];
% green=[0 1 0];
% red=[1 0 0];
% blue=[0 0 1];
% VegetationColormap=[white; green; red; blue];
% % Plot image
% h=subplot(1,1,1);
% imagesc(TC)
% set(h,'Clim',[-0.5 3.5]);
% colormap(VegetationColormap);
% colorbar; set(gco,'Clim',[1 4]);
%%%%%%%%%%%%%%%%%%%%%DYNAMIC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

while Time < EndTime
    Time= Time+dt
%     Creates LITTER in the neighborhod of pine (8 neighbors)+ the pine
%     site itself
    if TC(Age>AgeMP)==1 %because in the first years pine do not create litter
    [x,y]=find(TC(2:end-1,2:end-1)==1); %consider including here litter deposition only after a certain age
    x=x+1;y=y+1;
    for i=1:length(x)
        Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)=Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)+lrate*dt;
%         !!! maybe consider adding litter in the neighbourhood of oak when it is dominant
%         and adult
    end
    end
    
    % SEED BANK CALCULATION ONLY ONCE A YEAR
    SB(1)=SBP1+SBP2+SeedFP*(1-canopyBank)*sum(sum(TC(Age>AgeMP)==1))+ReleaseSeeds;% TWO YEARS OF SEED LIFE; 1-canopybank is doing the same as canopy bank, i.e. *0.5
    %SB(1)=SB(1)-SeedLoss(1)*SB(1);      %not needed as pine seeds only
    %last 2 years and then die
    SB(2)=SB(2)+SeedFS*(sum(sum(TC==2)))-SeedLoss(2)*SB(2);           % LONG SEED LIFE
    SB(3)=SeedFQ*(sum(sum(TC(Age>AgeMO)==3)))+randi(BirdSeedN,1);     % NO MEMORY - !!!Before seed loss was 1 now there is no seed loss
    if SB(3)>0
        coordseed=randi(m,SB(3),2); % puts the seeds that arrive in random coordinates of the lattice
    end
    % accumulation of seeds in the canopy
    SBPC=SBPC+SeedFP*canopyBank*sum(TC(Age>AgeMP)==1); %canopyBank=% of the seeds that stay in the
    % canopy and accumulate over time
    
    % COLONIZATION OF AN EMPTY CELL AND MORTALITY FOR VEGETATED CELLS
    %------------------------------------------------------------------
    % At each time step tests who is going to colonize an empty cell in the lattice based on the
    % availiable seeds (seed production and seed bank)
   
%     % %!!!!!!!!!!!THIS TERM IS for test without litter
%     ProbG=[1 1 1]; 
%     
    for i = 1 : m
        for j=1:m
            test=rand*length(ProbG); %RANDOM NUMBER BETWEEN 0 AND THE NUMBER OF SPECIES (LENGTH(PROBg=3))
            if TC(i,j)==0 % colonization/germination
                
%                 COLONIZATION vs LITTER
%                 

                
                ProbG(1)=(maxG(1)+minG(1))/2+(maxG(1)-minG(1))/2*tanh((LitThreshP-Lit(i,j))/amp(1)) ...
                    -(maxG(1)-ProbPZeroL)*exp(-2/LitThreshP*exp(1)*Lit(i,j)); % PINE
                ProbG(2)=(maxG(2)+minG(2))/2+(maxG(2)-minG(2))/2*tanh((LitThreshS-Lit(i,j))/amp(2)); % SEEDER ampS=0.3 max=.9 min=0.
                ProbG(3)=maxG(3)-(maxG(3)-minG(3))*exp(-Lit(i,j)); % QUERCUS
                
                
                %%% TERM FOR PROBABILITY OF COLONIZATION vs. NUMBER OF SEEDS AND ESTABLISHMENT
                
                ProbS(1:2)=1-(1-1./est(1:2)).^(SB(1:2)/m/m); % FOR PINE AND SEEDERS, SEEDS ARE EQUALLY SPREAD THROUGHOUT THE CELLS; this was taken in the paper: Cannas et al. 2003
                % QUERCUS HAS PROB 1 IF THERE IS A SEED IN THAT CELL
                % VERSION 2
               
                ProbS(3)=0; 
                mm=find(sum(ismember(coordseed(:,1:2),[i,j]),2)>=2);
                ProbS(3)=ProbS(3)+round(length(mm)/(length(mm)+eps))*.8;
               
% %                 % VERSION 1
%                 for ii=1:size(coordseed,1)
%                     ProbS(3)=ProbS(3)+(coordseed(ii,1)==i&coordseed(ii,2)==j);
%                 end
%                 ProbS(3)=ProbS(3)>=1; % if there is one seed or more -> prob=1
                
               
                
                % combining probabilities of establishment due to litter and seed numbers
                
              
              
                
                ProbG=ProbG.*ProbS;  %
                
                % this step has the improved version (Mara's) of Alains' MSc thesis trick
                
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
              Lit(i,j)=eflit.*Lit(i,j); %10% of the accumulated litter is degraded each year and 90% remains
                
                if test< mort(TC(i,j))*dt
                    TC(i,j)=0;
                    Age(i,j)=0;
                end
            end
        end
        SBP2=SBP1; % PINE SEED BANK OF TWO YEARS BEFORE
        SBP1=SB(1);% PINE SEED BANK OF 1 YEAR BEFORE
    end

    
    
    %%% DISTURBANCE
    
   D=randi(10,1)*(Time>=12);%allows pine to reproduce before the first disturbance
 
    D1=D(Time);
    if D1== 1
        'fire';
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
    
%     imagesc(TC)
%     set(h,'Clim',[-0.5 3.5]);
%     colormap(VegetationColormap);
%     colorbar
%     
%     drawnow;pause
    
    %%% update abundance of different species in the lattice
    
    Pine=sum(sum(TC==1));
    Seeder=sum(sum(TC==2));
    Oak=sum(sum(TC==3));
    

    %%%%%%%%%%%%%%%% STORING AND VISUALIZATION %%%%%%%%%%%%%%%%%
%     StoreTime = StoreTime - Time; % I COMMENTED THIS BECAUSE YOU WANT TO
%     STORE EVERY TIME STEP SO IT' NOT USEFUL TO HAVE THIS EXTRA IF. TO BE
%     RESTORED (AND CHANGED DIMENSIONS OF THE VECTORS BELOW) IF DT<1 YEAR
%     OR YOU WANT TO SAVE EVERY E.G. 10 YEARS 
%     if StoreTime <= 0
        StorePine(NrStore) = Pine; % NOTICE THESE VECTORS ARE AS LONG AS ENDTIMES, and as wide as 1 (vector not matrices) NOT M BY M AS YOU DEFINED THEM.. -> REDIFINING ABOVE SHOULD TAKE LESS TIME -> LET ME KNOW!
        StoreSeeder(NrStore) = Seeder;
        StoreOak(NrStore) = Oak;
        VectorTime(NrStore)= Time;
        NrStore = NrStore+1;
%         StoreTime = StoreStep;
%     end %if StoreTime <= 0

end

StoreSpecies=[StorePine StoreSeeder StoreOak];
% xlswrite('Sp abundance pine,seeder,oak',StoreSpecies)

imagesc(TC)
set(h,'Clim',[-0.5 3.5]);
colormap(VegetationColormap);
colorbar
%%%Plotting over time
figure
plot(VectorTime,StorePine/m/m*100,VectorTime,StoreSeeder/m/m*100,VectorTime,StoreOak/m/m*100)
legend('Pine','Seeder','Resprouter')
set(gca,'fontsize',14);
set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
xlabel('Time (year)');
ylabel('Cover (%)');

% saveas(gcf,'figureTime.png','png')

% Creates movie - not working yet
%     imagesc(TC);
%     movie(Frame)=getframe;
%     Frame=Frame+1;
%writerObj = VideoWriter(imageTC) %constructs a VideoWriter object to write video data to an AVI file with Motion JPEG compression.
