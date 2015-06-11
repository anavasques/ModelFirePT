% PLOTTING THE FILES OF THE MULTIRUNS - AVERAGING

clear all

nruns=10;
m=100;
maxseedSeed=100:100000:1000000; % makes runs changing the parameter of SB (2) between the three values determined and keeping all other values fixed
pineTimeM=zeros(nruns,length(maxseedSeed));pineTimestd=pineTimeM;

for k=1:2%length(maxseedSeed)
    figure(k)
       for irun=1:nruns
           filename=strcat(['firePT_NOF_MAXSS_OS',num2str(maxseedSeed(k)),'_',num2str(irun),'.mat' ]);%loads the file in the current directory %%name of directory should correspond; can load many directories
           load (filename)% command to load only certain variables: load(filename,variables) e.g. %'StorePine','StoreSeeder','StoreOak','StoreLitter','VectorTime')
%            figure
           plot(VectorTime,StorePine/m/m*100,'b', VectorTime,StoreSeeder/m/m*100, 'r--.', VectorTime,StoreOak/m/m*100, 'g*','LineWidth',5)%, VectorTime,StoreAge,'gr')
           hold on % the info between hold on and hold off is added to the plot
           
           legend('Pine','Seeder','Oak')%, 'Average age')
           set(gca,'fontsize',20, 'fontWeight','bold');
           set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
           set(gca,'fontsize',16, 'fontWeight','bold');
           xlabel('Time (year)');
           ylabel ('Cover (%)');
           pause
%in this case eval stores the matrix for each run

%             eval(['matrixStorePine',int2str(k),'(:,irun)=StorePine;']) %this command (eval) is used to store the variables that can be plotted later
%             eval(['matrixStoreSeeder',int2str(k),'(:,irun)=StoreSeeder;']) 
%             eval(['matrixStoreOak',int2str(k),'(:,irun)=StoreOak;']) 
%             
%             % CALCULATE THE MEAN OVER A CERTAIN TIME SPAN
%             %%% this is an example many combinations can be made
%             pineTimeM(irun,k)=mean(StorePine(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
%             pineTimestd(irun,k)=std(StorePine(3:9));
%             seederTimeM(irun,k)=mean(StoreSeeder(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
%             seederTimestd(irun,k)=std(StoreSeeder(3:9));
%             oakTimeM(irun,k)=mean(StoreOak(3:9)); % the average and stdev are calculated between timei:timef for each value of irun and k
%             oakTimestd(irun,k)=std(StoreOak(3:9));
       end
%        legend('Pine','Seeder','Oak')%, 'Average age')
%        set(gca,'fontsize',14, 'fontWeight','bold');
%        set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
%        set(gca,'fontsize',16, 'fontWeight','bold');
%        xlabel('Time (year)');
%        ylabel ('Cover (%)');
%        title(filename)
       %title(['seed number=',int2str(maxseedSeed(k))])
     
%        eval(['matrix=matrixStorePine',int2str(k),';'])% stores the matrices for each k outside the loop - t calculare the average
%        eval(['matrix2=matrixStoreSeeder',int2str(k),';'])
%        eval(['matrix3=matrixStoreOak',int2str(k),';'])
%        % example to do average over a certain time only between different
%        % runs - not so useful!!
%        eval(['matrix(30:90,:)=matrixStorePine',int2str(k),';'])
% avPine=mean(matrix,2); %makes the mean per row the command simple (without 2) makes mean per column
% stdPine=std(matrix,0,2); %different command for std (than that one of mean) but does the same
% avSeeder=mean(matrix2,2); %makes the mean per row the command simple (without 2) makes mean per column
% stdSeeder=std(matrix2,0,2); %different command for std (than that one of mean) but does the same
% avOak=mean(matrix3,2); %makes the mean per row the command simple (without 2) makes mean per column
% stdOak=std(matrix3,0,2); %different command for std (than that one of mean) but does the same
% 
% plot (VectorTime, avPine/m/m*100,'b', VectorTime, avSeeder/m/m*100,'r--.', VectorTime, avOak/m/m*100,'g*', 'LineWidth',5)
%            
%            legend('Pine','Seeder','Oak')%, 'Average age')
%            set(gca,'fontsize',14, 'fontWeight','bold');
%            set(gcf,'Position',[374 407 981 410],'PaperPositionMode','auto');
%            set(gca,'fontsize',20, 'fontWeight','bold');
%            xlabel('Time (year)');
%            ylabel ('Average cover (%)');
%            pause
% 
% plot(VectorTime, avPine/m/m*100,'k')
% %errorbar(VectorTime,avPine/m/m*100,stdPine/m/m*100,'k') %same plot than before but with errorbar (stdev)
% plot(VectorTime, avSeeder/m/m*100,'k')
% %errorbar(VectorTime,avSeeder/m/m*100,stdSeeder/m/m*100,'k') %same plot than before but with errorbar (stdev)
% plot(VectorTime, avOak/m/m*100,'k')
% %errorbar(VectorTime,avOak/m/m*100,stdOak/m/m*100,'k') %same plot than before but with errorbar (stdev)

hold off
pause     
end

%figure
% plot(maxseedSeed,pineTimeM/m/m*100)

% THE SAME BUT WITH STD
figure
for i=1:nruns
    errorbar(maxseedSeed,pineTimeM(i,:)/m/m*100,pineTimestd(i,:)/m/m*100)
    errorbar(maxseedSeed,seederTimeM(i,:)/m/m*100,seederTimestd(i,:)/m/m*100)
    errorbar(maxseedSeed,oakTimeM(i,:)/m/m*100,oakTimestd(i,:)/m/m*100)
    hold on
end
        