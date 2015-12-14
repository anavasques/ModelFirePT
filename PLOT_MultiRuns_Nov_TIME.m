clear all

nruns=100;
m=100;
%keeping all other values fixed
BirdSeedNv=[1 5 20 100];
F=[0 7 15 30];

storeMatr=zeros(length(F)*length(BirdSeedNv)*nruns,5);
storeMatrAv=zeros(length(F)*length(BirdSeedNv),6);
iii=0;
ii=0;

% IMAGINE YOU SAVE THIS .M FILE IN THE SAME FOLDER WHERE THERE ARE ALL THE
% FOLDERS OF FIRE RECURRENCE THEN AFTER THE FIRST FOR LOOP YOU CAN MOVE IN
% THE FIRE RECURRENCE FOLDER OF YOUR CHOICE

for i=1:length(F) % fire recurrence 
    foldername=strcat('firePT_F',num2str(F(i)),'HS_BIRS_OS');
    % I don't know your foldername but you can here use the FR(i) to create your folder name
    cd(foldername) % to be tested for windows..
    for k=1:length(BirdSeedNv)
         ii=ii+1;
        for irun=1:nruns
         iii=iii+1;
            filename=strcat(['firePT_F',num2str(F(i)),'HS_BIRS_OS',num2str(BirdSeedNv(k)),'_',num2str(irun),'.mat' ]);%loads the file in the current directory %%name of directory should correspond; can load many directories
            load (filename)% command to load only certain variables: load(filename,variables) e.g. %'StorePine','StoreSeeder','StoreOak','StoreLitter','VectorTime')
            
            % this is not needed (to me, for the purpose of the plotting)
%             %saves a matrix for each variable and value of k with all the repeated runs
%             eval(['matrixStoreTime',int2str(k),'(:,irun)=VectorTime;'])
%             eval(['matrixStorePine',int2str(k),'(:,irun)=StorePine;']) %this command (eval) is used to store the variables that can be plotted later
%             eval(['matrixStoreSeeder',int2str(k),'(:,irun)=StoreSeeder;'])
%             eval(['matrixStoreOak',int2str(k),'(:,irun)=StoreOak;'])
            
            %%%%%%calculates the time when pine cover is equal to zero
            
% ANA, SORRY. THIS IS MAYBE CORRECT, BUT FOR A STRANGE THING WITH GITUHUB I
% DIDN' FIND THIS FILE AT FIRST AND MODIFIED IT MYSELF IN A DIFFERENT WAY..
% SORRY (SEE BELOW). IF YOU THINK YOUR NOTATION WORKS, YOU CAN KEEP IT...!
% targetPine = 0;
% XPine = StorePine;
% ind = find((StorePine==0));
% XtargetPine= XPine(ind);
% StoretimeTargetPine=Storetime(XtargetPine(1,1));
% pineTime(irun)= StoretimeTargetPine(irun);
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
                oakTime(irun)=NaN;
            else
                oakTime(irun)=vec(1);
            end
 	storeMatr(iii,1)=F(i);
        storeMatr(iii,2)=BirdSeedNv(k);
        storeMatr(iii,3)=irun;
        storeMatr(iii,4)=pineTime(irun);
        storeMatr(iii,5)=oakTime(irun);
        

        end

	

        avPineTime=mean(pineTime); %makes the mean per row the command simple (without 2) makes mean per column
        stdPineTime=std(pineTime); %different command for std (than that one of mean) but does the same
        avOakTime=mean(oakTime); %makes the mean per row the command simple (without 2) makes mean per column
        stdOakTime=std(oakTime); %different command for std (than that one of mean) but does the same
       
       storeMatrAv(ii,1)=F(i);
       storeMatrAv(ii,2)=BirdSeedNv(k);
       storeMatrAv(ii,3)=avPineTime;
       storeMatrAv(ii,4)=stdPineTime;
       storeMatrAv(ii,5)=avOakTime;
       storeMatrAv(ii,6)=stdOakTime;
        
        
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

for k=1:length(BirdSeedNv)
    figure(k)
    title(['Pine time =0; oak seeds= ',num2str(BirdSeedNv(k))])
    set(gca,'xtick',1:4,'xticklabel',['  7 '; ' 15 '; ' 30 '; '2000'],'ylim',[0 500])
    figure(100+k)
    title(['Oak time <50; oak seeds= ',num2str(BirdSeedNv(k))])
    set(gca,'xtick',1:4,'xticklabel',['  7 '; ' 15 '; ' 30 '; '2000'],'ylim',[0 500])
end


save storeMatr.txt storeMatr -ASCII
save storeMatrAv.txt storeMatrAv -ASCII

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


