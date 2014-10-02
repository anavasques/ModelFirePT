% FILE TO TRY TO IMPROVE firePT.m - OCT 2ND 2014
m=100;
TC= zeros(m,m);  % Creates a matrix of size m*m filled with zeros

% Fills the matrix TC with planted pines
%--------------------------------------------------------------------------
TC(4:4:96,4:4:96)= 1; % plants 1 pine every 4 meters - prodution stand
Pine=sum(sum(TC));

[x,y]=find(TC==1);

for i=1:length(x)
    Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)=Lit(x(i)-1:x(i)+1,y(i)-1:y(i)+1)+lrate*dt;
%     Lit(x(i),y(i))=Lit(x(i),y(i))-lrate*dt;
end

%  Si=randi(m,1);
%  Sj=randi(m,1);