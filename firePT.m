%Model PT fire 
%CASCADE project
%Asynchronous CA
%Vasques et al. 2014
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

%%% SEEDER
AgeMS=1;                  % Age of maturity seeder % field obs Calluna% [year]
                          %Cistus 3 years ref
SeedFS=1000;              % Seed production per plant/occupied cell approx value
                          %check in lit
LSS=30;                   % Life span calluna % in woodland education centre [year]

%%% OAK
AgeMO=50;                 % Age of maturity seeder % Kew % [year] 
%!!! here we assume that the oaks are always dominated in the understory
%so seed production needs correction
SeedFQ=10;                % Seed production oak per occupied cell - value that is not fundamented by literature assumes that the tree is dominated and small if under cover of pinus
                          % 120 acorns per tree refered in Martin?k et al. 2014% [n/m2/year] 
BirdSeedN=50;             % Annual seed input by birds - based on average values Q. suber Pons and Pausas 2007 - this value depends on surrounding populations
LSO= 1000;                % Life span quercus robur % in forestar

%GENERAL
minG= [0 0 0.3];          % minumum germination first pine second seeder third oak
maxG= [0.9 0.9 0.9];      % maximum germination first pine second seeder third oak
SeedLoss= [0.50 0.05 1];  % rate seed loss first pine second seeder third oak

lrate=0.42;               % rate of litter deposition [cm/year] Fernandes et al 2004
ProbG=[0 0 0];            % probability of germination first pine second seeder third oak
amp=[0 0.3 0];            % amplitude of curve interaction with litter

SBP1=0;                   % !SBP1 and SBP2 are only ways of initializing the seed bank every year
SBP2=0;
SBPC=0;                    % Inicialization seed bank pine canopy
SB= [0 0 0];              % seed bank first pine second seeder third oak
ProbS= [0 0 0];           % to calculate probability based on seed prod
maxsedl= [7 400 1];       % max number of seedlings per cell CCD field

mort=1./[LSP,LSS,LS0];    % MORTALITY of pine, seeder, oak = 1/lifespan


AR= [0,0,1];              % ability to resprout: pine=0, seeder=0, oak=1;

% Control constants
StartTime= 0;                   % [year]
EndTime= 150;                   % [year]

% NOTE: put dt smaller than one year in a way that the probabilities are <1 but not too small otherwise the model runs slowly
dt= 1;                          % [year]
m= 100;                         % for size of lattice [meter]

%Control variables
Time = StartTime;

D= 0;                           % initialization only the regime is defined underneath with interval of 10 years
Pine=0;                         % will count the number of cells with pine
Seeder=0;                       % will count the number of cells with pine
Oak=0;                          % will count the number of cells with pine

%%%Initialization of the matrices 
%TC - Type of Cover: 1- pine; 2- seeder; 3- oak; Age- plant age; SB- Seed Bank); Lit- Litter
% -------------------------------------------------------------------------
TC= zeros(m,m);  % Creates a matrix of size m*m filled with zeros       
Age= zeros(m,m); % Creates a matrix of size m*m filled with zeros    
Lit= zeros(m,m); % Creates a matrix of size m*m filled with zeros
z= 8;            % Number of neighbours

% Fills the matrix TC with planted pines
%--------------------------------------------------------------------------
TC(4:4:96,4:4:96)= 1; % plants 1 pine every 4 meters - dense prodution stand excluding the borders

%%%%%%%%%%%%%%%%%%%%%DYNAMIC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
while Time < EndTime
    for k=1:round(1/dt)
        % Creates litter in the neighborhod of pine (8 neighbors)+ the pine
        % site itself
        [x,y]=find(TC==1);
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
                    
                    ProbG=ProbG*dt; %this is the trick to get probability smal
                  
                    ProbS=1.*SB>0; % this term includes the seeds (ProbS) into the probability of establishment
                    
                    %%% Term for the relation between available seeds and
                    %%% Prob S. careful that the probability is lower than
                    %%% 1!
                    
                    %ProbS(1)=SB(1)*1/maxSedl(1);
                    %ProbS(2)=SB(2)*1/maxSedl(2); 
                    %ProbS(3)=SB(3)*1/maxSedl(3);
                    
                    ProbG=ProbG.*ProbS;
                    
                    if sum(ProbG)*dt>1
                    %this step has the reference  of Alains' MSc thesis
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
    SB(2)=SB(2)+SeedFS*(sum(sum(TC==2)))-SeedLoss(2)*SB(2); % LONG SEED LIFE
    SB(3)=SeedFQ*(sum(sum(TC(Age>AgeMO)==3)))+randi(BirdSeedN,1); % NO MEMORY
    SB(3)=SB(3)-SeedLoss(3)*SB(3);
    
    SBP2=SBP1; % PINE SEED BANK OF TWO YEARS BEFORE
    SBP1=SB(1);% PINE SEED BANK OF 1 YEAR BEFORE
   
    % accumulation of seeds in the canopy
    SBPC=SeedFP*canopyBank*sum(TC(Age>AgeMP)==1); %canoyBank=% of the seeds that stay in the
    % canopy and accumulate over time
    
    %%% DISTURBANCE
    D=randi(10,1);
    if D == 1 % if Time/10 is an integer there is disturbance % interval of 10 years
        Lit(:,:)=0;
        for i=1:m
            for j=1:m
                TC(i,j)= TC(i,j)*AR(TC(i,j)); % check - the idea is to have 0 for pine and seeder and 1 for oak
                Age(i,j)= Age(i,j)*AR(TC(i,j)); % check - the idea is to have 0 for pine and seeder and 1 for oak
                %SBP1=SBC
            end
        end
    end
end
    
    %%% update abundance of different species in the lattice
%   Pine=sum(sum(TC==1));
%   Seeder=sum(sum(TC==2));
%   Oak= sum(sum(TC==3));
    
    % %Plot the relative abundance of plant species at each time step
%     % This is done now with TC but could be done with Pine, Seeder and Oak
    % white=[1 1 1];
    % lightgreen=[0.5 1 0.5];
    % green=[0 1 0];
    % darkgreen=[0 0.5 0.5];
    % TCmap=[white; green; lightgreen; darkgreen];
    % %Plot image
    % colormap(TCmap);
    % color
    % showimage(TC);
    % drawnow;
    
    % %Plot the relative abundance of plant species over time
    % Creates colormap
        figure
        white=[1 1 1];
        green=[0 1 0];
        red=[1 0 0];
        blue=[0 0 1];
        VegetationColormap=[white; green; red; blue];
        % Plot image
        colormap(VegetationColormap);
    %     %color 
    %     N=4;
        showimage(TC);
    %     L = line(ones(N),ones(N), 'LineWidth',2);
    %     set(L,{'color'},mat2cell(VegetationColormap,ones(1,N),3));
    %     legend('empty','Pine','Seeder', 'Resprouter');
        drawnow;

    % Creates movie
    % showimagesc(TC);
    % movie(Frame)=getframe;
    % Frame=Frame+1;

                
                