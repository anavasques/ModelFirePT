clear all

nruns=100;
m=100;
%keeping all other values fixed
BirdSeedNv=[1 5 20 100];
FR=[7 15 30 2000];

% IMAGINE YOU SAVE THIS .M FILE IN THE SAME FOLDER WHERE THERE ARE ALL THE
% FOLDERS OF FIRE RECURRENCE THEN AFTER THE FIRST FOR LOOP YOU CAN MOVE IN
% THE FIRE RECURRENCE FOLDER OF YOUR CHOICE

for i=1:length(FR) % fire recurrence 
    foldername=strcat% I don't know your foldername but you can here use the FR(i) to create your folder name
    cd(foldername) % to be tested for windows..
    for k=1:length(BirdSeedNv)
         
        for irun=1:nruns
            filename=strcat(['firePT_FIRE',num2str(FR(i),'LS_BIRS_FS',num2str(BirdSeedNv(k)),'_',num2str(irun),'.mat' ]);%loads the file in the current directory %%name of directory should correspond; can load many directories
            load (filename)% command to load only certain variables: load(filename,variables) e.g. %'StorePine','StoreSeeder','StoreOak','StoreLitter','VectorTime')
            
            % this is not needed (to me, for the purpose of the plotting)
%             %saves a matrix for each variable and value of k with all the repeated runs
%             eval(['matrixStoreTime',int2str(k),'(:,irun)=VectorTime;'])
%             eval(['matrixStorePine',int2str(k),'(:,irun)=StorePine;']) %this command (eval) is used to store the variables that can be plotted later
%             eval(['matrixStoreSeeder',int2str(k),'(:,irun)=StoreSeeder;'])
%             eval(['matrixStoreOak',int2str(k),'(:,irun)=StoreOak;'])
            
            %%%%%%calculates the time when pine cover is equal to zero
            pineTime(irun)=(StoreTime(StorePine==0));
            
            %%%%%%calculates the time when oak cover is equal or higher than
            %%%%%%50%
            oakTime(irun)=(StoreTime(StoreOak>=50));
        end
        avPineTime=mean(pineTime); %makes the mean per row the command simple (without 2) makes mean per column
        stdPineTime=std(pineTime); %different command for std (than that one of mean) but does the same
        avOakTime=mean(oakTime); %makes the mean per row the command simple (without 2) makes mean per column
        stdOakTime=std(oakTime); %different command for std (than that one of mean) but does the same
        
        figure(k)
        bar(i,avPineTime), hold on
        errorbar(i,avPineTime,stdPineTime)
        xlabel('Fire recurrence (years)');
        ylabel ('Cover (%)');
        title(['oak seeds= ',num2str(k)])
        figure(100+k)
        bar(FR(i),avOakTime), hold on
        errorbar(FR(i),stdOakTime)
        xlabel('Fire recurrence (years)');
        ylabel ('Cover (%)');
    end
end

for k=1:length(BirdSeedNv)
    figure(k)
    set(gca,'xtick',1:4,'xticklabel',['  7 '; ' 15 '; ' 30 '; '2000'])
    figure(100+k)
    set(gca,'xtick',1:4,'xticklabel',['  7 '; ' 15 '; ' 30 '; '2000'])
end
%     %%% stores the matrices outside the loop - to later calculate
%     %%% averages
%       this needs to be corrected to make it work
%     eval(['pinezero=pineTime',int2str(k),';'])
%     eval(['oak50=oakTime',int2str(k),';'])
    
%     filename=strcat(['pine_and_oak',num2str(BirdSeedNv(k)),'.mat' ]);
%     avPineTime=mean(pinezero); %makes the mean per row the command simple (without 2) makes mean per column
%     stdPineTime=std(pinezero,0); %different command for std (than that one of mean) but does the same
%     avOakTime=mean(oak50); %makes the mean per row the command simple (without 2) makes mean per column
%     stdOakTime=std(oak50,0); %different command for std (than that one of mean) but does the same
%     
%     save(filename,'avPineTime','stdPineTime','avOakTime','stdOakTime')
%     
%     
% end
%gets the files and makes the plotting

% for k=1:length(BirdSeedNv)
%     
%     filename=strcat(['pine_and_oak',num2str(BirdSeedNv(k)),'.mat']);
%     load (filename,'avPineTime','stdPineTime','avOakTime','stdOakTime')
%     
%     
%     figure(100) %here we can put any number, but a number should be in the parenthesis for it to open always the same figure and plot all in the same figure
%     plot (FR,'avPineTime')
%     hold on % hold on should be after the plot
%     legend('1','5','20','100')%, 'Average age')
%     set(gca,'fontsize',14, 'fontWeight','bold');
%     set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
%     set(gca,'fontsize',16, 'fontWeight','bold');
%     xlabel('Fire recurrence (years)');
%     ylabel ('Cover (%)');
%     axis([0 2000 0 1000]);
%     errorbar(FR,avPineTime,stdPineTime,'k') %same plot than before but with errorbar (stdev)
%     title('Average time when pine is locally estinguished')
%     
%     pause
% end
% hold off


