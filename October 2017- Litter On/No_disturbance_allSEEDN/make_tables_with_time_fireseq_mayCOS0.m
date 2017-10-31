%%%% File to make tables and figures with different fire sequences

clear all
close all

nruns=20;
m=100;
%keeping all other values fixed
BirdSeedNv=[1 2 5 50];
F = 0;
b =[122:141]; %b is the fire sequence

storeMatrOS0=zeros(length(F)*length(b)*nruns,6);
storeMatrAvOS0=zeros(length(F)*length(b),6);
iii=0;
ii=0;

% SAVE THIS .M FILE IN THE SAME FOLDER WHERE THERE ARE ALL THE
% FOLDERS OF FIRE RECURRENCE THEN AFTER THE FIRST FOR LOOP YOU CAN MOVE IN
% THE FIRE RECURRENCE FOLDER OF YOUR CHOICE

for i=1:length(BirdSeedNv) % fire recurrence 
    foldername= strcat('firePT_', num2str(BirdSeedNv(i)),'F0_BIRS_OS'); %if file exists
    % I don't know your foldername but you can here use the FR(i) to create
    % your folder name
% figure(i)
    cd(foldername) % tested for windows..
    for k=1:length(b)
        ii=ii+1;
        for irun=1:nruns
         iii=iii+1;
        
       filename=strcat(['firePT_F0_BIRS_OS',num2str(BirdSeedNv(i)),'_',num2str(irun),'_',num2str(b(k)),'.mat' ]);%loads the file in the current directory %%name of directory should correspond; can load many directories
%        if (exist(filename,'file')==0)% this command is for the cases
%        where various files are missing
%            continue
%        else
       load (filename)% command to load only certain variables: load(filename,variables) e.g. %'StorePine','StoreSeeder','StoreOak','StoreLitter','VectorTime')
           % this is not needed (to me, for the purpose of the plotting)
           %             %saves a matrix for each variable and value of k with all the repeated runs
           %             eval(['matrixStoreTime',int2str(k),'(:,irun)=VectorTime;'])
           %             eval(['matrixStorePine',int2str(k),'(:,irun)=StorePine;']) %this command (eval) is used to store the variables that can be plotted later
           %             eval(['matrixStoreSeeder',int2str(k),'(:,irun)=StoreSeeder;'])
           %             eval(['matrixStoreOak',int2str(k),'(:,irun)=StoreOak;'])
%         plot (VectorTime,StorePine/m/m*100,'b', VectorTime,StoreSeeder/m/m*100,'r--.', VectorTime, StoreOak/m/m*100,'g*')
%         hold on % hold on should be after the plot
%         %legend('Pine','Seeder','Oak')%, 'Average age')
%         set(gca,'fontsize',14, 'fontWeight','bold');
%         set(gcf,'Position',[2 2 750 500],'PaperPositionMode','auto');
%         set(gca,'fontsize',16, 'fontWeight','bold');
%         xlabel('Time (year)');
%         ylabel ('Cover (%)');
%         axis([0 1000 0 100]);
%         title('firePT_F7LS_BIRS_OS')
           %pause

        %%%%%%calculates the time when pine cover is equal to zero
            vec=VectorTime(StorePine==0);
            if isempty(vec)
                pineTime(irun)=NaN; % if you don' like NaN you could also put -9999 (any negative number I mean)
            else
                pineTime(irun)=vec(1);
            end
            
            %%%%%%calculates the time when oak cover is equal or higher than
            %%%%%%50%
            vec=VectorTime(StoreOak>=50);
            if isempty(vec)
                oakTime50(irun)=NaN;
            else
                oakTime50(irun)=vec(1);
            end
            
             vec=VectorTime(StoreOak>=80);
            if isempty(vec)
                oakTime80(irun)=NaN;
            else
                oakTime80(irun)=vec(1);
            end
            
 	storeMatrFS0(iii,1)= BirdSeedNv(i);
        storeMatrFS0(iii,2)=b(k);
        storeMatrFS0(iii,3)=irun;
        storeMatrFS0(iii,4)=pineTime(irun);
        storeMatrFS0(iii,5)=oakTime50(irun);
        storeMatrFS0(iii,6)=oakTime80(irun);
        

       %end %End of if
       end % end of for

        avPineTime=mean(pineTime); %makes the mean per row the command simple (without 2) makes mean per column
        stdPineTime=std(pineTime); %different command for std (than that one of mean) but does the same
        avOakTime=mean(oakTime50); %makes the mean per row the command simple (without 2) makes mean per column
        stdOakTime=std(oakTime50); %different command for std (than that one of mean) but does the same
        
       storeMatrAvFS0(ii,1)=BirdSeedNv(i);
       storeMatrAvFS0(ii,2)=b(k);
       storeMatrAvFS0(ii,3)=avPineTime;
       storeMatrAvFS0(ii,4)=stdPineTime;
       storeMatrAvFS0(ii,5)=avOakTime;
       storeMatrAvFS0(ii,6)=stdOakTime;
        
        
%         figure(k)
%         
%         bar(i,avPineTime), hold on
%         errorbar(i,avPineTime,stdPineTime)
%         xlabel('Fire recurrence (years)');
%         ylabel ('Time (y)');
%         
%         figure(100+k)
%         bar(i,avOakTime), hold on
%         errorbar(i,avOakTime,stdOakTime)
%         xlabel('Fire recurrence (years)');
%         ylabel ('Time (y)');
    end
    cd ..
end

% for k=1:length(b)
%     figure(k)
%     title(['Pine time =0; oak seeds= ',num2str(b(k))])
%     set(gca,'xtick',1:4,'xticklabel',['  7 '; ' 15 '; ' 30 '; '2000'],'ylim',[0 500])
%     figure(100+k)
%     title(['Oak time <50; oak seeds= ',num2str(b(k))])
%     set(gca,'xtick',1:4,'xticklabel',['  7 '; ' 15 '; ' 30 '; '2000'],'ylim',[0 500])
% end


save storeMatrOS0.txt storeMatrOS0 -ASCII
save storeMatrAvOS0.txt storeMatrAvOS0 -ASCII

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


