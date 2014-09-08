%% How do the interactions between fire regime and pre-fire structure
%% determine potential shifts in plant community composition?
%% Model PT fire CASCADE developed by Vasques et al. 2014


%%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%Constants

%Control constants
StartTime= 0;                   % [year]
EndTime= 100;                   % [year]
dt= 1;                          % [year]
m= 100;                         % [meter]

%%%Control variables
Time = StartTime;
Dist= 
Pine=0 % to be able to count the number of cells with each type and introduce the formulas of AgeMat, SeedProd, LifeSpan
Seeder=0
Oak=0

%%%System parameters
AgeMatP=
SeedProdP=
%MaxSeedLongP
%SeedSizeP
%SeedAgeP
LifeSpanP=

AgeMatS=
SeedProdS=
%MaxSeedLongS
%SeedSizeS
%SeedAgeS
LifeSpanS=

AgeMatO=
SeedProdO=
%MaxSeedLongO
%SeedSizeO
%SeedAgeO
LifeSpanO=
ResType= [0,1];

% Probability factor of seed mortality (yearly)
MortP=0.5;
MortS=0.1;
MortO=0.8;

%%%Initialisation of the lattice
PlantCover= zeros(m,m);         %1=pine; 2=oak; 3=seeder (first only seeder and then use it as resprouter modofy mortality)
PlantAge= zeros(m,m);
SeedBank= zeros(m,m);
LitterCover= zeros(m,m);

% CHECK HOW TO DO THIS FOR THE LITTER CASE
%Create convolution matrix
%ConvMatrix=[0 0 f 0 0; 0 0 e 0 0;f e 0 e f; 0 0 e 0 0; 0 0 f 0 0];
%local environment parameters
%e= %closer neigh.
%f= %distant neigh.



% Create colormap
white=[1 1 1];
lightgreen=[0.5 1 0.5];
green=[0 1 0];
darkgreen=[0 0.5 0.5];
PlantCovermap=[white; lightgreen; darkgreen; darkgreen];
% Plot image
colormap(PlantCovermap);
%color 
showimage(PlantCover);
drawnow;

%%%%%%%%%%% DYNAMIC %%%%%%%%%%%%%%%%%%%%%%

while Time < EndTime
    
for i=1:100
    for j= 1:100
        for Row=3:(m-2)
        for Col=3:(m-2)
            
if Dist= 1; 
  LitterCover (i,j)=0;
  PlantCover (i,j)= PlantCover (i,j)* ResType
  PlantAge (i,j)= PlantAge (i,j)* ResType
else
    LitterCover(i,j)= LitterCover (i-1,j-1)+ % include factor of accumulation a * PineCover
    PlantCover (i,j)= PlantCover (i,j)
    PlantAge (i,j)= PlantAge (i-1, j-1)+ dt
end
    
    
    %colonization
    
  for i=1:100
    for j= 1:100
        
        %PlantCover=1 % in a pecentage of the cells
            
        if PlantCover=0; 
                
                %check the seedbank and calculate a probability of
                %establishment (equation of the interaction of each seed type with
                %Littercover);
                
                %update plant cover as a result of colonization
                
            else PlantCover
                if PlantCover(Row-1,Col-1)==1
                    Pine=Pine+1;
                elseif PlantCover(Row-1,Col-1)==2
                    Seeder=Seeder+1;
                elseif PlantCover(Row-1,Col-1)==3
                    Oak=Oak+1;
                end
    
                % calculate litter for the local environment of pine (start
                % accumulating after pine is 4 years old? until a maximum depth)
                %local environment
for i=1:100
    for j= 1:100
            %this is only for pine
            %%%SIMPLIFY%% LitterCover=[min(1,PlantCover(Row-2,Col-2)) min(1,PlantCover(Row-2,Col-1)) min(1,PlantCover(Row-2,Col)) min(1,PlantCover(Row-2,Col+1)) min(1,PlantCover(Row-2,Col+2)); min(1,PlantCover(Row-1,Col-2)) min(1,PlantCover(Row-1,Col-1)) min(1,PlantCover(Row-1,Col)) min(1,PlantCover(Row-1,Col+1)) min(1,PlantCover(Row-1,Col+2)); min(1,PlantCover(Row,Col-2)) min(1,PlantCover(Row,Col-1)) min(1,PlantCover(Row,Col)) min(1,PlantCover(Row,Col+1)) min(1,PlantCover(Row,Col+2)); min(1,PlantCover(Row+1,Col-2)) min(1,PlantCover(Row+1,Col-1));min(1,PlantCover(Row+1,Col)) min(1,PlantCover(Row+1,Col+1)) min(1,PlantCover(Row+1,Col+2)); min(1,PlantCover(Row+2,Col-2)) min(1,PlantCover(Row+2,Col-1)) min(1,PlantCover(Row+2,Col)) min(1,PlantCover(Row+2,Col+1)) min(1,PlantCover(Row+2,Col+2))]; 
            %ConvLocEnv=LitterCover.*ConvMatrix;
            %Change(Row,Col)=round(max(-1,min(1,(c*sum(sum(ConvLocEnv))))));
                
                % Produce seeds - local input? use local environment
                if PlantAgeP (i,j)>= AgeMatP; SeedProdP;
                else SeedProdP=0
                    if PlantAgeS (i,j)>= AgeMatS; SeedProdS;
                else SeedProdS=0
                    if PlantAgeO (i,j)>= AgeMatO; SeedProdO;
                else SeedProdO=0
                    end
                % put seeds in the soil seed bank % number of seeds of each
                % species per cell need different lattices?
                SeedBank(i,j)= SeedBank(ij)= SeedBank(i-1,j-1)+SeedProd-SeedBank(i,j)*MortP)
                
                % update plant cover
                
                % not for now	
                %Check SeedAge; if > MaxSeedLong then increase SeedMortality; else SeedBank
                %Check PlantAge to get resprouting ability


                
                    end
                    
                    
  Time= Time+dt
                
                