%Model PT fire
%CASCADE project
%Asynchronous CA
%Vasques et al. 2014-2015 (model development)
% ------------------------------------------------------------------------

close all
clear all

%Parameter values
%-------------------------------------------------------------------------
%PER SPECIES

%%% PINE
AgeMP=10;                 % Age of maturity pine 
%SeedFP=1000;             % Seed production per pine mature tree number of seeds per cone (63)* cone per tree (15) Vega et al 2008
SeedFP=100;               % Reduced in 10 times as in Quercus
LSP= 100;                 % Life span of pine % in "practices centro pinus" [year]
canopyBank=0.5;           % Percent of the seeds that are stored in the canopy maybe reduce it too many pines
% ReleaseSeeds=0;           % Pine seeds in the canopy that are released after the fire% not used in the current version


%%% SEEDER
AgeMS=3;                  % Age of maturity seeder % field obs Calluna 1 [year]; Cistus 3 years ref
SeedFS=1000;
%SeedFS=2;                % Seed production per plant/occupied cell approx value ADJUST
LSS=30;                   % Life span calluna % in woodland education centre [year]


%%% OAK
%AgeMO=50;                % Age of maturity % Kew % [year]% !Pausas 1999 has maturity = 15!
AgeMO=20;                 % According to Ramon an oak can produce acorns after 15-20 years
%SeedFQ=12;               % Seed production oak per occupied cell - 120 acorns per tree refered in Martin?k et al. 2014% [n/m2/year]
SeedFQ=10;
%SeedFQ=0;                % If FQ=0 oak does not produce seeds, it creates a reserve of saplings in the understory
%BirdSeedN=0;
BirdSeedN=5;              % Annual seed input by birds - average values Q. suber Pons and Pausas 2007 50seeds per hectar - this value depends on surrounding populations
%BirdSeedN=500;           % to experiment
RespAge=1;                % ONLY OAKS OLDER THAN THIS AGE CAN RESPROUT
%RespAge=10;              % RESPROUT ABILITY at X years - for experiments

%LSO= 1000;               % Life span quercus robur % in forestar;
LSO=500;                  % this was reduced by half to be more realistic (also in the calculations of mortality)


%LITTER RELATED
minG= [0 0 0.3];          % minimum germination first pine second seeder third oak
maxG= [0.9 0.9 0.9];      % maximum germination first pine second seeder third oak
ProbPZeroL=0.7;           % Germination probability for pine when litter=0 cm
LitThreshP=3;             % Litter threshold for Pine above which ~no germination (cm)
LitThreshS=2;             % Litter threshold for seeders above which ~no germination (cm)
lrate=0.8;                % These values were approximated to have a curve close to the sigmoid- indication from literature: rate of litter deposition [cm/year] Fernandes et al 2004 have 0.42; Indication from Ramon: after 20-30 litter stabilizes
eflit=1-0.4;              % effective litter: if 0.90 then 0.10 of the total litter is decomposed - estimated value not from literature
%lrate=0.05

amp=[0.3 0.3 0];          % amplitude of curve interaction with litter


LitOn=1;                  % switch for litter on/off
lconv=0.5*ones(3,3);lconv(2,2)=1; %convolution matrix for the litter around pine
%%expanded matrix for litter
%lconv=0.25*ones(5,5);lconv(2,2)=0.5;lconv(2,3)=0.5;lconv(2,4)=0.5;lconv(3,2)=0.5;lconv(3,3)=1;lconv(3,4)=0.5;lconv(4,2)=0.5;lconv(4,3)=0.5;lconv(4,4)=0.5;
litLS=4;
% GENERAL
est= [10 100 2];          % max number of seedlings per cell CCD field from which we inferred a probability of establishment in one cell
%est= [7 400 7];          % infered probability of establishment for ProbS

SeedLoss= [0 0.10 1];     % rate seed loss soil seed bank 1 pine 2 seeder 3 oak
nrsp=length(SeedLoss);    % number of species used in the model - to put in the prob expression
mort=1./[LSP,LSS,LSO];    % MORTALITY of pine, seeder, oak = 1/lifespan
AR= [0,0,0,1];            % Ability to resprout: first element is fake (bare soil); pine=0, seeder=0, oak=1;
z= 8;                     % Number of neighbours
%D=0;                      % initialization of disturbance

% % FIGURE
% white=[1 1 1];
% blue=[0 0 1];
% red=[1 0 0];
% green=[0 1 0];
% VegetationColormap=[white; blue;red;green];

% CONTROL CONSTANTS AND VARIABLES
StartTime= 0;             % [year]
EndTime= 10;             % [year]
StoreTime = 1;            % [year]

dt=1;                     % [year]
m= 100;                   % for size of lattice [meter]

nruns=3;
maxseedSeed=100:100000:1000000; % makes runs changing the parameter of SB (2) between the three values determined and keeping all other values fixed
% BirdSeedNv=1:5:50; 
% pd=3:10:43

%%%%% CODE FOR MULTIRUNS %%%%%%%
%%%%BirdSeed
% for k=1:length(BirdSeedNv)
%     BirdSeedN=BirdSeedNv(k);
%     SB=[0 100*m*m 0+randi(BirdSeedN,1)];
%%%%%Initial pine density
% for k=1:length(pd)
%     TC(pd:pd:m-pd,pd:pd:m-pd)=1
%%%%MaxSeedSeeder
        %%%VECTOR OF FIRE OCCURRENCE
%%%VECTOR OF FIRE OCCURRENCE
tf=100;                     %time without fires
%tf=40000    %no disturbance

fireret=15;                %interval between fires - fire return

for b=121:131
D=0*[StartTime:dt:EndTime];%#ok<NBRAK>
rand('state',b)           % if you move this to l.107 (right before the loop over the runs) you will have also different fire sequences in the runs

while tf<EndTime% EndTime can be substituted for the time when disturbance should stop
    %%%%%%%%%% CHECK WITH MARA %%%%%%%%%%%%%%%%%%%%%
    %%Multiruns for fire return
    % varfireret= 4000,30,15,7;
    % for r=1:length(varfireret)
    % for l=1:10
    % tf=tf-varfireret(l,:)*log(rand(1,1));
    % D(round(tf))=1
    %end
    %end
    tf=tf-fireret*log(rand(1,1)); %stochastic fire recurrence Baudena et al 2010
    D(round(tf))=1;
end

for k=1:length(maxseedSeed)
    SB=[0 maxseedSeed(k) 0+randi(BirdSeedN,1)];
      % for  l=1:length(BirdSeedN)
    for irun=1:nruns
        k,irun
      
        %     filename=strcat(['par',num2str(k),'.mat']);
        %     save(filename)
        
        %%%------------------------------------------------------------------------
        %%%INICIALIZATION OF THE MATRICES
        
        %TC - Type of Cover: 1- pine; 2- seeder; 3- oak; Age- plant age; SB- Seed Bank); Lit- Litter
        % -------------------------------------------------------------------------
        TC= zeros(m,m);           % Creates a matrix of size m*m filled with zeros
        Age= zeros(m,m);          % Creates a matrix of size m*m filled with zeros
        Lit= zeros(m+2,m+2);      % Creates a matrix of size (m+2)*(m+2) filled with zeros - the borders should be excluded in stats done with litter
        PosQSeed=zeros(m,m);      % NUMBER OF QUERCUS SEEDS PER CELL
        %     Seeder=0;                 % will count the number of cells with seeder
        %     Oak=0;                    % will count the number of cells with oak
        %     Pine=0;                   % will count the number of cells with pine
        %
        ProbL=[0 0 0];            % probability of germination due to litter (first pine second seeder third oak)
        %     Litter=0;                 % to sum the number of cells with litter
        ProbG=[0 0 0];            % probability of germination first pine second seeder third oak
        ProbS= [0 0 0];           % to calculate probability based on seed prod
        
        %PLANT PINES
        %TC(2:2:m-2,2:2:m-2)= 1;   % plants 1 pine every X meters - dense prodution stand excluding the borders
        TC(3:3:m-3,3:3:m-3)=1;  % pine is planted every 3 meters, there is no gap
        %between pines - homogeneous when canopy closes
        %TC(40:40:m-40,40:40:m-40)= 1; % plants 1 pine every X meters - for
        %experiments
        
        Time = StartTime;
        NrStore = 1;
        StoreStep = StoreTime;
        StorePine=zeros(EndTime,1);
        StoreSeeder=zeros(EndTime,1);
        StoreOak=zeros(EndTime,1);
        StoreLitter=zeros(EndTime,1);
        VectorTime=zeros(EndTime,1);
        StoreMatPine=zeros(EndTime,1);
        
        SBP1=0;                   % !SBP1 and SBP2 are only ways of initializing the seed bank every year
        SBP2=0;
        SBPC=0;                   % Initialization seed bank pine canopy
        
        %     MatPine=0;
        
        
        %Puts cover of seeder or oak randomly in the lattice
        % s=1000;%puts a number of cells occupied with seeder or oak(in this case seeders in a random manner)
        % rp=randperm(m*m,s);
        % TC(rp)=2;
        
        %SB=[0 100*m*m 0+randi(BirdSeedN,1)]; %NOT for pine!! initial seed bank %comment this on the multiruns
        %SB=[0 0 0]; %starting seeds of oak and seeder changing to analyse one
        %species at a time
        %initial conditions for seeder and oak, pine is planted but can also be seeded randomly
        %SBP1=100*m*m %to start the seeds of pine
        
        
% %%%The Fire part was here -> MB: SINCE FOR THESE EXPERIMENTS YOU WANT TO
% HAVE THE SAME FIRE SERIES, YOU NEED TO CREAD D ONLY ONCE (SEE ABOVE). IF
% YOU CREATE IT HERE WITH ONLY ONE INITIALISATION (UP) EACH RUN WILL HAVE
% DIFFERENT FIRE SEQUENCES -> THAT' WHY I MOVED IT UP
        
        %--------------------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%%%DYNAMIC LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %--------------------------------------------------------------------------
        
        while Time < EndTime
            Time= Time+dt;
            %%%% Creates LITTER in the neighborhod of pine (8 neighbors)+ pine
            [x,y]=find(TC==1); %finds cells =1 in the whole matrix - already has if the AgeMP is not needed anymore because the accumulation only starts after 10 y anyway
            x=x+1;y=y+1;
            for i=1:length(x)
                %%%%%% Closer to a sigmoid curve
                indAge=Age(x(i)-1,y(i)-1)>AgeMP;
                Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)=Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)+lrate*lconv*dt*(indAge+(1-indAge)*.2); % .2 IS THE 20% OF THE MAXIMUM VALUE OF LITTER DEPOSITION RATE FOR PINE<AGEMP
            end
            
            % SEED BANK CALCULATION (ONCE A YEAR)
            
            %PINE soil seeds (SB1) and canopy seeds (SBPC)
            SB(1)=SBP1+SBP2+SeedFP*(1-canopyBank)*sum(sum(TC(Age>AgeMP)==1));% Pine TWO YEARS OF SEED LIFE; 1-canopybank is doing the same as canopy bank, i.e. *0.5
            SBPC=SBPC+SeedFP*canopyBank*sum(TC(Age>AgeMP)==1);
            %SEEDER
            SB(2)=SB(2)+SeedFS*(sum(sum(TC==2)))-SeedLoss(2)*SB(2);           % LONG SEED LIFE
            %OAK
            SB(3)=SeedFQ*(sum(sum(TC(Age>AgeMO)==3)))+randi(BirdSeedN,1);
            %%%%%UNCOMENT TO ELIMINATE OAKS
            %SB(3)=0; %no oaks in the lattice
            %%%%%% NEW WAY
            for kk=1:SB(3) %only happens if SB3 is bigger than 1
                cc=randi(m,1,2);%c2=randi(m,1,1);
                PosQSeed(cc(1),cc(2))=PosQSeed(cc(1),cc(2))+1;
            end
            
            %------------------------------------------------------------------
            % COLONIZATION OF AN EMPTY CELL AND MORTALITY FOR VEGETATED CELLS
            %------------------------------------------------------------------
            % At each time step tests who is going to colonize an empty cell in the lattice based on the
            % availiable seeds (seed bank) and on litter
            for i = 1 : m
                for j=1:m
                    test=rand*nrsp; %RANDOM NUMBER BETWEEN 0 AND THE NUMBER OF SPECIES (LENGTH(PROBg=3))
                    if TC(i,j)==0 % colonization/germination
                        
                        %%% TERM FOR PROBABILITY OF COLONIZATION vs. NUMBER OF SEEDS AND ESTABLISHMENT
                        
                        ProbS(1)=1-(1-1/est(1))^(SB(1)/m/m);
                        ProbS(2)=1-(1-1/est(2))^(SB(2)/m/m); % FOR PINE AND SEEDERS, SEEDS ARE EQUALLY SPREAD THROUGHOUT THE CELLS; this was taken in the paper: Cannas et al. 2003
                        % only for oak the seeds are spread over the lattice
                        ProbS(3)=1-(1-1/est(3))^PosQSeed(i,j); %if ProsQSeed=0 the whole term goes to zero
                        
                        %%%COLONIZATION vs LITTER
                        
                        ProbL(1)=(maxG(1)+minG(1))/2+(maxG(1)-minG(1))/2*tanh((LitThreshP-Lit(i+1,j+1))/amp(1)) ...
                            -(maxG(1)-ProbPZeroL)*exp(-2/LitThreshP*exp(1)*Lit(i+1,j+1)); % PINE
                        ProbL(2)=(maxG(2)+minG(2))/2+(maxG(2)-minG(2))/2*tanh((LitThreshS-Lit(i+1,j+1))/amp(2)); % SEEDER ampS=0.3 max=.9 min=0.
                        ProbL(3)=maxG(3)-(maxG(3)-minG(3))*exp(-Lit(i+1,j+1)); % QUERCUS
                        ProbL=ProbL*LitOn+(1-LitOn)*[1 1 1];% term for test IN ABSENCE of litter
                        
                        % COMBINING PROBABILITIES OF ESTABLISHMENT DUE TO SEED
                        % NUMBERS AND LITTER % version Mara modified from Alain's i.e. prob only needs to
                        % be <than 3 and we don't need to multiply it by dt anymore)
                        
                        ProbG=ProbL.*ProbS; %term for test with litter
                        
                        if test>ProbG(1)+ProbG(2)+ProbG(3)
                            TC(i,j)=0; %this was modified to make the code faster so if the cell is not colonized it stops running here
                        elseif test>ProbG(1)+ProbG(2)
                            TC(i,j)=3;
                            Age(i,j)=dt;
                        elseif test>ProbG(1)
                            TC(i,j)=2;
                            Age(i,j)=dt;
                            SB(2)=SB(2)-100; % FOR EVERY ADULT SHRUB SEEDER OCCUPYING A CELL 100 (?) SEEDS ARE LOST
                        else
                            TC(i,j)=1;
                            Age(i,j)=dt;
                            SB(1)=SB(1)-10;
                        end
                    else
                        Age(i,j)=Age(i,j)+dt;
                        
                        if test< mort(TC(i,j))*dt %determines if a cell dies, mort is defined according to life span
                            TC(i,j)=0;
                            Age(i,j)=0;
                        end
                    end
                end
            end
            
 %%%% DISTURBANCE
 % the rest of the disturbance term is at the beggining of the code
    D1=D(Time);
    if D1== 1
        'fire';
%         %HIGH SEVERITY
%         Lit(:,:)=0;
        %LOW SEVERITY CHANGE PAR NAME PLEASE
        Lit(:,:)=litLS*sum(sum(TC(Age>AgeMP)==1))/m/m*8; %leaves from the canopy fall creating a litter
        %layer - for simplification purposes the litter in the soil is
        %maintained - it is porportional to the canopy cover anyway - canopy effect on 8
        %neightbours
        for i=1:m
            for j=1:m
                % PINES AND SEEDERS DIE; QUERCUS RESPROUTS IF OLDER OR
                % EQUAL THEN Respage
                % Kind IS A TRICK TO AVOID IF STRUCTURE (IF AGE>RespAge THEN RESPROUT)
                kind=floor(AR(TC(i,j)+1)*Age(i,j)/RespAge);
                kind=round(kind/(kind+eps));
                TC(i,j)= TC(i,j)*kind;
                Age(i,j)= Age(i,j)*kind;
            end
        end
               % PINE CANOPY SEEDS FALL INTO SEEDBANK
        SB(1)= SBPC; %the production of seeds when there is a fire is the total of the canopy seeds produced until that moment
        SBP1=0;SBPC=0; %SBP2=0; Not needed to put sbp2 to zero
    else
        Lit=eflit*Lit; % effective litter, i.e. litter that is not degraded and remain for the years after - should be updated here to also occur in the empty cells
    end
            
            % UPDATE OF THE PINE SEEDBANK
            SBP2=0.5*SBP1; % PINE SEED BANK OF TWO YEARS BEFORE IS 50%
            SBP1=SB(1);% PINE SEED BANK OF 1 YEAR BEFORE
            
            % RESETS NUMBER OF QUERCUS SEEDS TO ZERO FOR NEXT YEAR, REINITIALISE
            % EVERY YEAR THE NUMBER OF SEEDS COMING IN
            PosQSeed=0*PosQSeed;      % NUMBER OF QUERCUS SEEDS PER CELL
            
            %%% UNCOMMENT TO HAVE OAKS
            SB(3)=randi(BirdSeedN,1); %This is the term to get a new random number between 1-5 every year
            
            % Store variables for plotting
            Pine=sum(sum(TC==1));
            MatPine=sum(sum(TC==1&Age>AgeMP));
            Seeder=sum(sum(TC==2));
            Oak=sum(sum(TC==3));
            Litter=sum(sum(Lit(2:m+1,2:m+1)));
            %Litter=max(max(Lit(2:m+1,2:m+1)));%
            
            %AgeMTX=mean(mean(Age));
            
            StorePine(NrStore) = Pine; % NOTICE THESE VECTORS ARE AS LONG AS ENDTIMES, and as wide as 1 (vector not matrices)
%             StoreMatPine(NrStore)= MatPine;
            StoreSeeder(NrStore) = Seeder;
            StoreOak(NrStore) = Oak;
            StoreLitter(NrStore)= Litter;
            %StoreAge(NrStore)=AgeMTX; %correct this
            VectorTime(NrStore)= Time;
            NrStore = NrStore+1;
            StoreTime = StoreStep;
            
        end
        %%%%%% MULTIRUNS CODE
        filename=strcat(['nofire_maxseed',num2str(maxseedSeed(k)),'_',num2str(irun),'.mat' ]);
        save(filename,'StorePine','StoreSeeder','StoreOak','StoreLitter','VectorTime')
        %     filename=strcat(['seedstart',num2str(maxseedSeed(k)),'.mat' ]);
        %     save(filename,'matr','-ascii') %save a text file

        %%% Plots final figure of vegetation
% figure
% h=subplot(1,1,1);
% imagesc(TC)
% set(h,'Clim',[-0.5 3.5]);
% set(gca,'FontSize',20,'fontWeight','bold')
% colormap(VegetationColormap);
% colorbar
% drawnow; %pause

% if LitOn==1
    %%%figure for Litter depth over time
%     figure
%     set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
%     plot(VectorTime,StoreLitter/m/m)
%     xlabel('Time (year)');
%     ylabel ('Mean litter depth (cm)');
%     
% end %if Lit
% saveas(gcf,'figureTime.png','png')
    end % stochastic runs
%     % %%%%% MULTIRUNS Figures
%         filename=strcat(['nofire_maxseed',num2str(maxseedSeed(k)),'_',num2str(irun),'.mat' ]);
%         save(filename,'VectorTime','mPineMR','mSeederMR','mOakMR')
%         %     filename=strcat(['seedstart',num2str(maxseedSeed(k)),'.mat' ]);
%         %     save(filename,'matr','-ascii') %save a text file
%     %%%Plotting over time
%     figure
%     plot(VectorTime,mPineMR/m/m*100,'b', VectorTime,mSeederMR/m/m*100, 'r--.', VectorTime,mOakMR/m/m*100, 'g*')%, VectorTime,StoreLitter/m/m, 'k.')%, VectorTime,StoreAge,'gr')
%     %%%%ERROR BARS NOT WORKING YET
% %     errorbar(VectorTime,mPineMR/m/m*100,EP)
% %     errorbar(VectorTime,mSeederMR/m/m*100,ES)
% %     errorbar(VectorTime,mOakMR/m/m*100,EO)
%     legend('Pine','Seeder','Oak')%, 'Litter mean depth')%, 'Average age')
%     set(gca,'fontsize',14, 'fontWeight','bold');
%     set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
%     set(gca,'fontsize',16, 'fontWeight','bold');
%     xlabel('Time (year)');
%     ylabel ('Cover (%)');
%     title('fire_maxseed')
%     filename=strcat(['fire_maxseed',num2str(maxseedSeed(k)),'_',num2str(irun),'.tif' ]); %saves the files
%     saveas(gcf,filename,'tif')% saves one figure for each output
    %end
end % maxseeds
end
    
