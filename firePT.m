%Model PT fire CASCADE
%Asynchronous CA
%Vasques et al. 2014
% ------------------------------------------------------------------------

close all
clear all

%Parameter values
% ------------------------------------------------------------------------

%Pine
%-------------------------------------------------------------------------
% %AgeMP=10;     % Age of maturity pine %start at 6 and regularly 10-15 % in Cronk and Fuller, 1995 [year] 
% %SeedPP=0.5;   % Seed production pine per occupied cell % should prob increase with age [n/m2/year]
% SeedLongP=2;  % Seed longevity pine one year % in "Dean et al. (1986)" [year]
% LSP= 100;     % Life span of pine % in "practices centro pinus" [year]
% %Seeder
% %-------------------------------------------------------------------------
% AgeMS=1;       % Age of maturity seeder % field obs % [year] 
% SeedPS=500;    % Seed production seeder per occupied cell %invented value% increases with age [n/m2/year]
% %SeedLongS=200;% Seed longevity seeder %invented value [year]
% %mSS=XX         % Define mortality of seeds of seeder species
% LSS=30;        % Life span calluna % in woodland education centre [year]
% %Oak
% %-------------------------------------------------------------------------
% AgeMO=50;       % Age of maturity seeder % Kew % [year] 
% SeedPO=120;     % Seed production oak per occupied cell isolated tree in Martin?k et al. 2014% [n/m2/year]
% SeedLongO=1;    %Does not subsist from one year to another %based on literature [year]
% LSO= 1000;      % Life span quercus robur % in forestar

%%% !!!!!!!!!
%%%%% INITIALIZE ALL THE PARAMETERS

mort=1./[LSP,LSS,LS0]; % MORTALITY of Pine, seeder, Oak = 1/lifespan
AR= [0,0,1]; % ability to resprout: pine=0, seeder=0, oak=1;

%------- Probability of establishment - max and min
maxG=[0.5 0.5 0.5];
minG=[0.1 0.1 0.1];

%Control constants
StartTime= 0;                   % [year]
EndTime= 150;                   % [year] % TO SOLVE: time frame is problematic because the quercus will also become trees and start producing litter
dt= 1;                          % [year]
m= 100;                         % [meter]

%Control variables
%--------------------------------------------------------------------------
Time = StartTime;
D=Time/10; % disturbance regime defined with ifs
Pine=0; % will count the number of cells with pine
Seeder=0; % will count the number of cells with pine
Oak=0; % will count the number of cells with pine

%%%Initialisation of the matrices 
%TC - type cover; Age- plant age; SB- Seed Bank); Lit- Litter
% -------------------------------------------------------------------------
TC= zeros(m,m);  % Creates a matrix of size m*m filled with zeros       
Age= zeros(m,m); % Creates a matrix of size m*m filled with zeros    
SBS= zeros(m,m);  % Creates a matrix of size m*m filled with zeros  
Lit= zeros(m,m); % Creates a matrix of size m*m filled with zeros
z= 8;            % Number of neighbours

% Fills the matrix TC with planted pines
%--------------------------------------------------------------------------
TC(4:4:96,4:4:96)= 1; % plants 1 pine every 4 meters - prodution stand
Pine=sum(sum(TC));
      
%Dynamic
%--------------------------------------------------------------------------
while Time < EndTime
    % GERMINATION AND MORTALITY EACH TIME STEP
    for k=1:round(1/dt) 
        % Creates litter in the neighborhod of pine (8 neighbors)+ the pine site itself
        % ---------------------------------------------------------------------------
        [x,y]=find(TC==1);
        for i=1:length(x)
            Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)=Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)+lrate*dt;
        end
        
        
        % Colonize an Empty cell & mortality for vegetated cells
        %-----------------------------------------------------------------------------------------------------
        % Test who is going to colonize an empty cell in the lattice based on the
        % seed bank
        
        %Include here a term for the equation of the interaction of each seed type with
        %Litter in the cell and decide who wins;
        
        for i = 1 : m
            for j=1:m
                test=rand;
                if TC(i,j)==0 % colonization/germination
                    % LITTER DEPENDENCE:
                    ProbG(1)=(maxG(1)+minG(1))/2;%(maxG(1)+minG(1))/2+(maxG(1)-minG(1))/2*tanh((Lit(i,j)-2)/amp(1)); % PINE
                    ProbG(2)=(maxG(2)+minG(2))/2+(maxG(2)-minG(2))/2*tanh((2-Lit(i,j))/ampS); % SEEDER ampS=0.3 max=.9 min=0.
                    ProbG(3)=maxG(3);%(maxG(3)+minG(3))/2+(maxG(3)-minG(3))/2*tanh((Lit(i,j)-2)/amp(3)); % QUERCUS
                    
                    ProbG=ProbG*dt;
                    
                    % SEED DEPENDENCE
                    % DEFINE THIS RELATION FOR EACH SPECIES TAKING INTO
                    % ACCOUNT THE SPATIAL OCUPATION AND ALSO THE RESULT OF
                    % THE PROBABILITY IN THE TERM THAT MULTIPLIES WITH
                    % LITTER PROBABILITY (GETS SMALLER) - THINK
                    
                    %ProbG(2)=ProbG*%MULTIPLY BY PROBABILITY DUE TO SEED NUMBER
                    
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
    %include the values for these parameters
    
    
    SB(1)=SBP1+SBP2+SeedFP*sum(sum(TC(Age>AgeMaturityPine)==1));
    SB(2)=SB(2)+SeedFS*(sum(sum(TC==2)))-SeedLoss*SB(2);
    SB(3)=SeedFQ*(sum(sum(TC==3)))+randi(BirdSeedN,1);
    
    SBP2=SBP1;
    SBP1=SB(1);
    
  
    %%% DISTURBANCE
    %if D == (int); % if Time/10 is an integer there is disturbance % interval of 10 years
    %Lit(i,j)=0;
    %TC(i,j)= TC(i,j)*AR; % check - the idea is to have 0 for pine and seeder and 1 for oak
    %Age(i,j)= Age(i,j)*AR; % check - the idea is to have 0 for pine and seeder and 1 for oak
    %else
    %    Lit(i,j)= Lit(i-1,j-1)+ 1; % this value should be the factor of accumulation of litter according to real values
    %   TC(i,j)= TC(i,j);
    %  Age(i,j)= Age(i-1, j-1)+ dt;
    %end
    

   
    % % Plot the CA at each time step
    % %--------------------------------------------------------------------------------------
    % white=[1 1 1];
    % lightgreen=[0.5 1 0.5];
    % green=[0 1 0];
    % darkgreen=[0 0.5 0.5];
    % TCmap=[white; green; lightgreen; darkgreen];
    % % Plot image
    % colormap(TCmap);
    % %color
    % showimage(TC);
    % drawnow;
    
    
    
end

                
                