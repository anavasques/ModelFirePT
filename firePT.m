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
AgeMP=10;     % Age of maturity pine %start at 6 and regularly 10-15 % in Cronk and Fuller, 1995 [year] 
SeedPP=0.5;   % Seed production pine per occupied cell % should prob increase with age [n/m2/year]
SeedLongP=2;  % Seed longevity pine one year % in "Dean et al. (1986)" [year]
LSP= 100;     % Life span of pine % in "practices centro pinus" [year]
%Seeder
%-------------------------------------------------------------------------
AgeMS=1;       % Age of maturity seeder % field obs % [year] 
SeedPS=500;    % Seed production seeder per occupied cell %invented value% increases with age [n/m2/year]
%SeedLongS=200;% Seed longevity seeder %invented value [year]
%mSS=XX         % Define mortality of seeds of seeder species
LSS=30;        % Life span calluna % in woodland education centre [year]
%Oak
%-------------------------------------------------------------------------
AgeMO=50;       % Age of maturity seeder % Kew % [year] 
SeedPO=120;     % Seed production oak per occupied cell isolated tree in Martin?k et al. 2014% [n/m2/year]
SeedLongO=1;    %Does not subsist from one year to another %based on literature [year]
LSO= 1000;      % Life span quercus robur % in forestar

%mP=0.5;        % Random probability factor mortality pine (yearly)
%mS=0.5;        % Random probability factor mortality pine (yearly)
%mO=0.2;        % Probability factor seed mortality pine   (yearly)

mort=1./[LSP,LSS,LS0]; % MORTALITY of Pine, seeder, Oak = 1/lifespan
AR= [0,0,1]; % ability to resprout: pine=0, seeder=0, oak=1;

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
TC(4:4:100,4:4:100)= 1; % plants 1 pine every 4 meters - prodution stand
Pine=sum(sum(TC));          
           
%Dynamic
%--------------------------------------------------------------------------

while Time < EndTime
% Creates litter in the neighborhod of pine (8 neighbors)  
% ---------------------------------------------------------------------------
[x,y]=find(TC==1);
for i=1:length(x)
    Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)=Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)+lrate*dt;    
end
% for i = 1 : m
%     for j = 1 : m
%         z = mod(j-1+m,m) ;
%         if z==0
%             z=m;
%         end
%         if TC(i,z) == 1
%             Lit(i,j) = Lit(i,j)+1 ;
%         end
%         z = mod(i-1+m,m) ;
%         if z==0
%             z=m;
%         end
%         if TC(z,j) == 1
%             Lit(i,j) = Lit(i,j)+1 ;
%         end
%         z = mod(j+1,m) ;
%         if z==0
%             z=m;
%         end
%         if TC(i,z) == 1
%             Lit(i,j) = Lit(i,j)+1 ;
%         end
%         z = mod(i+1,m);
%         if z==0
%             z=m;
%         end
%         if TC(z,j) == 1
%             Lit(i,j) = Lit(i,j)+1 ;
%         end
%     end
% end  

            
            
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
            % ProbG(1)=; % PINE
            ProbG(2)=(maxG(2)+minG(2))/2+(maxG(2)-minG(2))/2*tanh(2-Lit(i,j)/amp(2)); % SEEDER amp(2)=0.3
            % ProbG(3)=;% QUERCUS
            
            % SEED DEPENDENCE
            %ProbG(2)=ProbG*%PROBABILITY DUT O SEED NUMBER
            
            if sum(ProbG)*dt>1
                'sum of probability higher than 1! Please decrease dt'
                break
            else
                ProbG=ProbG/sum(ProbG);
                if test< litterfun
                %%%%
                end
            elseif test< mort(TC(i,j))*dt
               TC(i,j)=0; 
            end
        end
    end
end
        
%if TC==0;
% SONIA'S PART OF CODE: TESTS TO COLONIZE AN EMPTY CELL   
% identification of a cell randomly picked in the lattice
  
%end
%end
    
%test=rand; % how to do this? there is a seedbank from the seeder but not
%from the other species
%SBS % seed bank from the seeder
%SPP % seed production from year before
%SPO % seed production from year before

%Change(Row,Col)=2;
            
%Change(Row,Col)=1;

%Change(Row,Col)=3;
           
%Update number of cells occupied per each species
 %if Change (Row,Col)==1;
 %   TC(Row-1,Col-1)=1;
 %  Pine=Pine+1;
 %elseif Change (Row,Col)==2;
 %    TC(Row-1,Col-1)=2;
 %   Seeder=Seeder+1;
 %elseif Change (Row,Col)==3;
 %    TC(Row-1,Col-1)=3;
 %   Oak=Oak+1;              
 %else TC=TC;
 %end

%Update age
 %if Change > 0 % Starts Age in the occupied cells
 %   Age=1
 %else Age=0;
 %end
       
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
 
 
% Seed production - based on local input or overall input? Dependent on the
% local population??
%---------------------------------------------------------------------------------
%for i=1:m
%for j= 1:m                  
%if TC(i,j)==1 and Age>AgeMP; % CHECK AND!! produces seeds in the cells occupied by
%pine that are mature
%SeedPP=SeedPP;
%else SeedPP=0
%if TC(i,j)==2 and Age>AgeMS; % CHECK AND!! produces seeds in the cells occupied by
%seeder that are mature
%else SeedPS=0
%if TC(i,j)==2 and Age>AgeMO; % CHECK AND!! produces seeds in the cells occupied by
%oak that are mature
%else SeedPS=0
%end

% Seeds go to soil seed bank
%--------------------------------------------------------------------------------------
%SBS= SBS+SPS-mSS; % should this value be included in the random test for colonization or is
%more correct to d the test independently from the number of seeds
%available


% Plot the CA at each time step
%--------------------------------------------------------------------------------------
white=[1 1 1];
lightgreen=[0.5 1 0.5];
green=[0 1 0];
darkgreen=[0 0.5 0.5];
TCmap=[white; green; lightgreen; darkgreen];
% Plot image
colormap(TCmap);
%color 
showimage(TC);
drawnow;


Time= Time+dt
end

                
                