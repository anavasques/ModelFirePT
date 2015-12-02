clear all

nruns=100;
m=100;
eeping all other values fixed
BirdSeedNv=[1 5 20 100];
for FR=[7 15 30 2000]; % fire recurrence - needs to be taken from the names of the files
for k=1:length(BirdSeedNv)
       for irun=1:nruns
           filename=strcat(['./firePT_FIRE30LS_BIRS_FS',num2str(BirdSeedNv(k)),'_',num2str(irun),'.mat' ]);%loads the file in the current directory %%name of directory should correspond; can load many directories
           load (filename)% command to load only certain variables: load(filename,variables) e.g. %'StorePine','StoreSeeder','StoreOak','StoreLitter','VectorTime')
           
           %saves a matrix for each variable and value of k with all the repeated runs
           eval(['matrixStoreTime',int2str(k),'(:,irun)=VectorTime;'])
           eval(['matrixStorePine',int2str(k),'(:,irun)=StorePine;']) %this command (eval) is used to store the variables that can be plotted later
           eval(['matrixStoreSeeder',int2str(k),'(:,irun)=StoreSeeder;']) 
           eval(['matrixStoreOak',int2str(k),'(:,irun)=StoreOak;'])
           
           %%%%%%calculates the time when pine cover is equal to zero
           pineTime(irun,k)=(StoreTime(StorePine==0));
           
           %%%%%%calculates the time when oak cover is equal or higher than
           %%%%%%50%
           oakTime(irun,k)=(StoreTime(StoreOak>=50));
       end

end
       %%% stores the matrices outside the loop - to later calculate
       %%% averages
       eval(['pinezero=pineTime',int2str(k),';'])
       eval(['oak50=oakTime',int2str(k),';'])
   
filename=strcat(['pine_and_oak',num2str(BirdSeedNv(k)),'.mat' ]);
avPineTime=mean(pinezero); %makes the mean per row the command simple (without 2) makes mean per column
stdPineTime=std(pinezero,0); %different command for std (than that one of mean) but does the same
avOakTime=mean(oak50); %makes the mean per row the command simple (without 2) makes mean per column
stdOakTime=std(oak50,0); %different command for std (than that one of mean) but does the same

save(filename,'avPineTime','stdPineTime','avOakTime','stdOakTime')
    
      
       end
%gets the files and makes the plotting

for k=1:length(BirdSeedNv)
  
filename=strcat(['pine_and_oak',num2str(BirdSeedNv(k)),'.mat']);
load (filename,'avPineTime','stdPineTime','avOakTime','stdOakTime')


figure(100) %here we can put any number, but a number should be in the parenthesis for it to open always the same figure and plot all in the same figure
plot (FR,'avPineTime')
hold on % hold on should be after the plot
legend('1','5','20','100')%, 'Average age')
       set(gca,'fontsize',14, 'fontWeight','bold');
       set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
       set(gca,'fontsize',16, 'fontWeight','bold');
       xlabel('Fire recurrence (years)');
       ylabel ('Cover (%)');
       axis([0 2000 0 1000]);
errorbar(FR,avPineTime,stdPineTime,'k') %same plot than before but with errorbar (stdev)
title('Average time when pine is locally estinguished')

pause
end
hold off


