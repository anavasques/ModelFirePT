%this script should be inside the folder with the files to be plotted - one
%figure per folder (each condition of fire frequency, resprouting ability
%of shrub and fire severity is plotted in each of the figures)

close all
clear all

nruns=20;
m=100;
%maxseedSeed=100:100000:1000000; % makes runs changing the parameter of SB (2) between the three values determined and keeping all other values fixed
b=[142:1:151]
%pineTimeM=zeros(nruns,length(maxseedSeed));pineTimestd=pineTimeM;

    for k=1:length(b)
     figure(10)
        for irun=1:nruns
       filename=strcat(['firePT_F7HS_BIRS_OS5','_',num2str(irun),'_',num2str(b(k)),'.mat' ]);%loads the file in the current directory %%name of directory should correspond; can load many directories
%        if (exist(filename,'file')==0)% this command is for the cases
%        where various files are missing
%            continue
%        else
           load (filename)

      plot (VectorTime,StorePine/m/m*100,'b', VectorTime,StoreSeeder/m/m*100,'r--.', VectorTime, StoreOak/m/m*100,'g*')
        hold on % hold on should be after the plot
        %legend('Pine','Seeder','Oak')%, 'Average age')
        set(gca,'fontsize',14, 'fontWeight','bold');
        set(gcf,'Position',[2 2 750 500],'PaperPositionMode','auto');
        set(gca,'fontsize',16, 'fontWeight','bold');
        xlabel('Time (year)');
        ylabel ('Cover (%)');
        axis([0 2000 0 100]);
        %title('XXX')
        %pause
    end
end
%            %saves a matrix for each variable and value of k with all the repeated runs
%            eval(['matrixStoreTime',int2str(k),'(:,irun)=VectorTime;'])
%            eval(['matrixStorePine',int2str(k),'(:,irun)=StorePine;']) %this command (eval) is used to store the variables that can be plotted later
%            eval(['matrixStoreSeeder',int2str(k),'(:,irun)=StoreSeeder;']) 
%            eval(['matrixStoreOak',int2str(k),'(:,irun)=StoreOak;'])
%            
% %            %%%%%%calculates a mean over a certain time span (e.g. year 3 and 9)
% %            pineTimeM(irun,k)=mean(StorePine(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
% %            pineTimestd(irun,k)=std(StorePine(3:9));
% %            seederTimeM(irun,k)=mean(StoreSeeder(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
% %            seederTimestd(irun,k)=std(StoreSeeder(3:9));
% %            oakTimeM(irun,k)=mean(StoreOak(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
% %            oakTimestd(irun,k)=std(StoreOak(3:9));
%        end

%        
%        %%% stores the matrices outside the loop - to later calculate
%        %%% averages
%        eval(['time=matrixStoreTime',int2str(k),';'])
%        eval(['matrix=matrixStorePine',int2str(k),';'])
%        eval(['matrix2=matrixStoreSeeder',int2str(k),';'])
%        eval(['matrix3=matrixStoreOak',int2str(k),';'])
%        
%    
% filename=strcat(['average_birds',num2str(BirdSeedNv(k)),'.mat' ]);
% avPine=mean(matrix,2); %makes the mean per row the command simple (without 2) makes mean per column
% stdPine=std(matrix,0,2); %different command for std (than that one of mean) but does the same
% avSeeder=mean(matrix2,2); %makes the mean per row the command simple (without 2) makes mean per column
% stdSeeder=std(matrix2,0); %different command for std (than that one of mean) but does the same
% avOak=mean(matrix3,2); %makes the mean per row the command simple (without 2) makes mean per column
% stdOak=std(matrix3,0,2); %different command for std (than that one of mean) but does the same
% 
% save(filename,'avPine','avSeeder','avOak','stdPine','stdSeeder','stdOak','VectorTime')
%     
%       
%        end
% %gets the files and makes the plotting
% 
% for k=1:length(BirdSeedNv)
%   
% filename=strcat(['average_birds',num2str(BirdSeedNv(k)),'.mat']);
% load (filename,'VectorTime','avPine','stdPine','avSeeder','stdSeeder','avOak','stdOak')
% 
% 
% 
% figure(100) %here we can put any number, but a number should be in the parenthesis for it to open always the same figure and plot all in the same figure
% plot (VectorTime,avPine/m/m*100,'b', VectorTime,avSeeder/m/m*100,'r--.', VectorTime,avOak/m/m*100,'g*')
% hold on % hold on should be after the plot
% legend('Pine','Seeder','Oak')%, 'Average age')
%        set(gca,'fontsize',14, 'fontWeight','bold');
%        set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
%        set(gca,'fontsize',16, 'fontWeight','bold');
%        xlabel('Time (year)');
%        ylabel ('Cover (%)');
%        axis([0 2000 0 100]);
% title('number seeds oak')
% 
% pause
% end
% hold off
% plot (VectorTime, avPine/m/m*100,'k')
% %errorbar(VectorTime,avPine/m/m*100,stdPine/m/m*100,'k') %same plot than before but with errorbar (stdev)
% plot(VectorTime, avSeeder/m/m*100,'k')
% %errorbar(VectorTime,avSeeder/m/m*100,stdPine/m/m*100,'k') %same plot than before but with errorbar (stdev)
% plot(VectorTime, avOak/m/m*100,'k')
%errorbar(VectorTime,avOak/m/m*100,stdPine/m/m*100,'k') %same plot than
%before but with errorbar (stdev)


%   figure
% plot(maxseedSeed,pineTimeM/m/m*100)
%
% %%% THE SAME BUT WITH STD
% figure
% for i=1:nruns
%     errorbar(maxseedSeed,pineTimeM(i,:)/m/m*100,pineTimestd(i,:)/m/m*100)
%     errorbar(maxseedSeed,seederTimeM(i,:)/m/m*100,seederTimestd(i,:)/m/m*100)
%     errorbar(maxseedSeed,oakTimeM(i,:)/m/m*100,oakTimestd(i,:)/m/m*100)
%     hold on
% end


