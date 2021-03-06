clear all

nruns=10;
m=100;
%maxseedSeed=100:100000:1000000; % makes runs changing the parameter of SB (2) between the three values determined and keeping all other values fixed
BirdSeedNv=1:5:50
%pineTimeM=zeros(nruns,length(maxseedSeed));pineTimestd=pineTimeM;

<<<<<<< Updated upstream:RUNS May 16 (now working) Maras computer/PLOT_IND_runSept.m
for k=1:length(BirdSeedNv)
       for irun=1:nruns
           filename=strcat(['./firePT_NOF_BIRS_OS',num2str(BirdSeedNv(k)),'_',num2str(irun),'.mat' ]);%loads the file in the current directory %%name of directory should correspond; can load many directories
           load (filename)% command to load only certain variables: load(filename,variables) e.g. %'StorePine','StoreSeeder','StoreOak','StoreLitter','VectorTime')
           
           %saves a matrix for each variable and value of k with all the repeated runs
           eval(['matrixStoreTime',int2str(k),'(:,irun)=VectorTime;'])
           eval(['matrixStorePine',int2str(k),'(:,irun)=StorePine;']) %this command (eval) is used to store the variables that can be plotted later
           eval(['matrixStoreSeeder',int2str(k),'(:,irun)=StoreSeeder;']) 
           eval(['matrixStoreOak',int2str(k),'(:,irun)=StoreOak;'])
           
%            %%%%%%calculates a mean over a certain time span (e.g. year 3 and 9)
%            pineTimeM(irun,k)=mean(StorePine(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
%            pineTimestd(irun,k)=std(StorePine(3:9));
%            seederTimeM(irun,k)=mean(StoreSeeder(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
%            seederTimestd(irun,k)=std(StoreSeeder(3:9));
%            oakTimeM(irun,k)=mean(StoreOak(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
%            oakTimestd(irun,k)=std(StoreOak(3:9));
       end

       
       %%% stores the matrices outside the loop - to later calculate
       %%% averages
       eval(['time=matrixStoreTime',int2str(k),';'])
       eval(['matrix=matrixStorePine',int2str(k),';'])
       eval(['matrix2=matrixStoreSeeder',int2str(k),';'])
       eval(['matrix3=matrixStoreOak',int2str(k),';'])
       
   
filename=strcat(['average_birds',num2str(BirdSeedNv(k)),'.mat' ]);
avPine=mean(matrix,2); %makes the mean per row the command simple (without 2) makes mean per column
stdPine=std(matrix,0,2); %different command for std (than that one of mean) but does the same
avSeeder=mean(matrix2,2); %makes the mean per row the command simple (without 2) makes mean per column
stdSeeder=std(matrix2,0); %different command for std (than that one of mean) but does the same
avOak=mean(matrix3,2); %makes the mean per row the command simple (without 2) makes mean per column
stdOak=std(matrix3,0,2); %different command for std (than that one of mean) but does the same

save(filename,'avPine','avSeeder','avOak','stdPine','stdSeeder','stdOak','VectorTime')
    
      
       end
%gets the files and makes the plotting

for k=1:length(BirdSeedNv)
  
filename=strcat(['average_birds',num2str(BirdSeedNv(k)),'.mat']);
load (filename,'VectorTime','avPine','stdPine','avSeeder','stdSeeder','avOak','stdOak')



figure(100) %here we can put any number, but a number should be in the parenthesis for it to open always the same figure and plot all in the same figure
%plot (VectorTime,avPine/m/m*100,'b', VectorTime,avSeeder/m/m*100,'r--.', VectorTime,avOak/m/m*100,'g*')
%plot (VectorTime, avPine/m/m*100,'k')
errorbar(VectorTime,avPine/m/m*100,stdPine/m/m*100,'k') %same plot than before but with errorbar (stdev)
% plot(VectorTime, avSeeder/m/m*100,'k')
errorbar(VectorTime,avSeeder/m/m*100,stdPine/m/m*100,'k') %same plot than before but with errorbar (stdev)
% plot(VectorTime, avOak/m/m*100,'k')
errorbar(VectorTime,avOak/m/m*100,stdOak/m/m*100,'k') %same plot than
hold on % hold on should be after the plot
%legend('Pine')
%legend('Seeder')
%legend('Oak')
legend('Pine','Seeder','Oak')
       set(gca,'fontsize',14, 'fontWeight','bold');
       set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
       set(gca,'fontsize',16, 'fontWeight','bold');
       xlabel('Time (year)');
       ylabel ('Cover (%)');
       axis([0 2000 0 100]);
title('NO DISTURBANCE - individual plotting')

pause
=======
for k=1:length(maxseedSeed)
    for irun=1:nruns
        filename=strcat(['./firePT_NOF_MAXSS_OS',num2str(maxseedSeed(k)),'_',num2str(irun),'.mat' ]);%loads the file in the current directory %%name of directory should correspond; can load many directories
        load (filename)% command to load only certain variables: load(filename,variables) e.g. %'StorePine','StoreSeeder','StoreOak','StoreLitter','VectorTime')
        
        %saves a matrix for each variable and value of k with all the repeated runs
        eval(['matrixStoreTime',int2str(k),'(:,irun)=VectorTime;'])
        eval(['matrixStorePine',int2str(k),'(:,irun)=StorePine;']) %this command (eval) is used to store the variables that can be plotted later
        eval(['matrixStoreSeeder',int2str(k),'(:,irun)=StoreSeeder;'])
        eval(['matrixStoreOak',int2str(k),'(:,irun)=StoreOak;'])
        
        %            %%%%%%calculates a mean over a certain time span (e.g. year 3 and 9)
        %            pineTimeM(irun,k)=mean(StorePine(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
        %            pineTimestd(irun,k)=std(StorePine(3:9));
        %            seederTimeM(irun,k)=mean(StoreSeeder(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
        %            seederTimestd(irun,k)=std(StoreSeeder(3:9));
        %            oakTimeM(irun,k)=mean(StoreOak(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
        %            oakTimestd(irun,k)=std(StoreOak(3:9));
    end
    
    
    %%% stores the matrices outside the loop - to later calculate
    %%% averages
    eval(['time=matrixStoreTime',int2str(k),';'])
    eval(['matrix=matrixStorePine',int2str(k),';'])
    eval(['matrix2=matrixStoreSeeder',int2str(k),';'])
    eval(['matrix3=matrixStoreOak',int2str(k),';'])
    
    
    filename=strcat(['average_maxseed',num2str(maxseedSeed(k)),'.mat' ]);
    avPine=mean(matrix,2); %makes the mean per row the command simple (without 2) makes mean per column
    stdPine=std(matrix,0,2); %different command for std (than that one of mean) but does the same
    avSeeder=mean(matrix2,2); %makes the mean per row the command simple (without 2) makes mean per column
    stdSeeder=std(matrix2,0); %different command for std (than that one of mean) but does the same
    avOak=mean(matrix3,2); %makes the mean per row the command simple (without 2) makes mean per column
    stdOak=std(matrix3,0,2); %different command for std (than that one of mean) but does the same
    
    save(filename,'avPine','avSeeder','avOak','stdPine','stdSeeder','stdOak','VectorTime')
    
    
>>>>>>> Stashed changes:Save_av_and_plotting_(allin1).m
end
% %gets the files and makes the plotting
% 
% for k=1:length(maxseedSeed)
%   
%     filename=strcat(['average_maxseed',num2str(maxseedSeed(k)),'.mat']);
%     load (filename)
%     
%     
%     hold on
%     figure()
%     plot (VectorTime,avPine/m/m*100,'b', VectorTime,avSeeder/m/m*100,'r--.', VectorTime,avOak/m/m*100,'g*')
%     legend('Pine','Seeder','Oak')%, 'Average age')
%     set(gca,'fontsize',14, 'fontWeight','bold');
%     set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
%     set(gca,'fontsize',16, 'fontWeight','bold');
%     xlabel('Time (year)');
%     ylabel ('Cover (%)');
%     axis([0 500 0 100]);
%     title(filename)
%     
%     pause
% end
% hold off
% % plot (VectorTime, avPine/m/m*100,'k')
% % %errorbar(VectorTime,avPine/m/m*100,stdPine/m/m*100,'k') %same plot than before but with errorbar (stdev)
% % plot(VectorTime, avSeeder/m/m*100,'k')
% % %errorbar(VectorTime,avSeeder/m/m*100,stdPine/m/m*100,'k') %same plot than before but with errorbar (stdev)
% % plot(VectorTime, avOak/m/m*100,'k')
% %errorbar(VectorTime,avOak/m/m*100,stdPine/m/m*100,'k') %same plot than
% %before but with errorbar (stdev)
% 
% 
% %   figure
% % plot(maxseedSeed,pineTimeM/m/m*100)
% %
% % %%% THE SAME BUT WITH STD
% % figure
% % for i=1:nruns
% %     errorbar(maxseedSeed,pineTimeM(i,:)/m/m*100,pineTimestd(i,:)/m/m*100)
% %     errorbar(maxseedSeed,seederTimeM(i,:)/m/m*100,seederTimestd(i,:)/m/m*100)
% %     errorbar(maxseedSeed,oakTimeM(i,:)/m/m*100,oakTimestd(i,:)/m/m*100)
% %     hold on
% % end


