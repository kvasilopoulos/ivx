%17 April 2018
%
%
% Alexandros Kostakis, Tassos Magdalinos and Michalis P. Stamatogiannis, 
% "Robust econometric inference for stock return predictability", 
% The Review of Financial Studies, Vol. 28, No. 5, 2014, pp. 1506-1553.

%Replication of Tables 6-9, 11-13
% This program requires the data files monthly.xlsx, quarterly.xlsx and 
% annual.xlsx to run, as well as the function ivxlh.m to be placed in the
% same folder.
%Ordering of variables in the excel files:
%
%monthly.xlsx: 
% 1         2       3       4       5       6       7       8       9    10	11      12      13               
% yyyy      D/E     LTY     D/Y     D/P     TBL     E/P     B/M     INF	DFY	NTIS	TMS     LOG EXCESS VW 
%
%quarterly.xlsx
% 1         2       3       4       5       6       7       8       9    10	11      12       13             14         
% yyyy      D/E     LTY     D/Y     D/P     TBL     E/P     B/M     INF	DFY	NTIS	TMS     LOG EXCESS VW   CAY   

clc;
clear; 
datam = xlsread('monthly.xlsx',1);
dataq = xlsread('quarterly.xlsx',1);


%Table 6: Panel A, Monthly Data 1927-2012
vnam6a={'DE', 	'LTY', 	'DY', 	'DP', 	'TBL', 	'EP', 	'BM', 'DFY', 	'NTIS', 	'TMS', 	'INF'};
vnam6b={'Aols', 	'Aivx', 	'IVXWald', 	'delta'};
data=datam;
y=data(:,13);

xv=[2,3,4,5,6,7,8,10,11,12,9];
resmatrix=zeros(11,4);
for i=1:(size(xv,2));
x=data(:,xv(i));
[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,1,0);
resmatrix(i,1)=Aols(2,1);
resmatrix(i,2)=Aivx(1,1);
resmatrix(i,3)=Wivx(1,1);
resmatrix(i,4)=corr2(2,1);

end;

disp('Table 6; Panel A: Monthly Data January 1927-December 2012');
sTable6A = array2table(resmatrix,'VariableNames',vnam6b);
disp([vnam6a' sTable6A]);


%Table 6: Panel B, Monthly Data 1952-2012
subsample2=301:1033;
y=data(subsample2,13);
xv=[2,3,4,5,6,7,8,10,11,12,9];
resmatrix=zeros(11,4);
for i=1:(size(xv,2));
x=data(subsample2,xv(i));


[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,1,0);

resmatrix(i,1)=Aols(2,1);
resmatrix(i,2)=Aivx(1,1);
resmatrix(i,3)=Wivx(1,1);
resmatrix(i,4)=corr2(2,1);

end;

disp('Table 6; Panel B: Monthly Data January 1952-December 2012');
sTable6A = array2table(resmatrix,'VariableNames',vnam6b);
disp([vnam6a' sTable6A]);


%Table 7: Panel A, Panel B, Quarterly Data 1927-2012
   
vnam6a={'DE', 	'LTY', 	'DY', 	'DP', 	'TBL', 	'EP', 	'BM', 	'DFY', 	'NTIS', 	'TMS', 	'INF'};
vnam6b={'Aols', 	'Aivx', 	'IVXWald', 	'delta'};
data=dataq;
y=data(:,13);
xv=[2,3,4,5,6,7,8,10,11,12,9];

resmatrix=zeros(11,4);
for i=1:(size(xv,2));
x=data(:,xv(i));


[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,1,0);

resmatrix(i,1)=Aols(2,1);
resmatrix(i,2)=Aivx(1,1);
resmatrix(i,3)=Wivx(1,1);
resmatrix(i,4)=corr2(2,1);

end;

disp('Table 7; Panel A: Quarterly Data Q1 1927-Q4 2012');
sTable6A = array2table(resmatrix,'VariableNames',vnam6b);
disp([vnam6a' sTable6A]);



%Table 7: Panel B, Quarterly Data 1952-2012
vnam8={'DE', 	'LTY', 	'DY', 	'DP', 	'TBL', 	'EP', 	'BM', 	'DFY', 	'NTIS', 	'TMS', 	'INF', 	'CAY'};
data=dataq;
xv=[2,3,4,5,6,7,8,10,11,12,9,14];
resmatrix=zeros(12,4);
for i=1:(size(xv,2));

    if xv(i)==14;        
    subsample2=102:345;
    else;
       subsample2=101:345;
    end;

    y=data(subsample2,13);
    x=data(subsample2,xv(i));
    

[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,1,0);

resmatrix(i,1)=Aols(2,1);
resmatrix(i,2)=Aivx(1,1);
resmatrix(i,3)=Wivx(1,1);
resmatrix(i,4)=corr2(2,1);

end;

disp('Table 7; Panel B: Quarterly Data Q1 1952-Q4 2012');
sTable7A = array2table(resmatrix,'VariableNames',vnam6b);
disp([vnam8' sTable7A]);




%Table 8: Panel A, Monthly Data 1927-2012

vnam4={'yyyymm', 	'DE', 	'LTY', 	'DY', 	'DP', 	'TBL', 	'EP', 	'BM', 	'INF', 	'DFY', 	'NTIS', 	'TMS', 	'JointWald'};
data=datam;
y=data(:,13);

xdata=[5 6 0 0 0 0 0; 5 6 10 12 0 0 0; 5 8 0 0 0 0 0; 5 2 0 0 0 0 0; 7 8 12 0 0 0 0; 7 6 0 0 0 0 0];
xdata1=[0 0 0 0 1 1 0 0 0 0 0 0;
        0 0 0 0 1 1 0 0 0 1 0 1; 
        0 0 0 0 1 0 0 1 0 0 0 0;
        0 1 0 0 1 0 0 0 0 0 0 0; 
        0 0 0 0 0 0 1 1 0 0 0 1; 
        0 0 0 0 0 1 1 0 0 0 0 0];

resmatrix=zeros(6,13);

for j=1:size(xdata1,1);
    predictors=find(xdata1(j,:));
    x=data(:,predictors );    
    
   [Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,1,0);
    resmatrix(j,predictors)=Aivx;
    resmatrix(j,13)=Wivx(1,1);
end;


disp('Table 8; Panel A: Monthly Data January 1927-December 2012');

sTable = array2table(resmatrix,'VariableNames',vnam4);
%exclude unnecessary columns from the Table;
sTable8A =horzcat(sTable(:,2),sTable(:,5:8),sTable(:,10),sTable(:,12:end));

disp(sTable8A);

%Table 8: Panel B, Monthly Data 1952-2012

vnam4={'yyyymm', 	'DE', 	'LTY', 	'DY', 	'DP', 	'TBL', 	'EP', 	'BM', 	'INF', 	'DFY', 	'NTIS', 	'TMS', 	'JointWald'};
data=datam;
subsample2=301:1033;
y=data(subsample2,13);
xdata=[5 6 0 0 0 0 0; 5 6 10 12 0 0 0; 5 8 0 0 0 0 0; 5 2 0 0 0 0 0; 7 8 12 0 0 0 0; 7 6 0 0 0 0 0];

resmatrix=zeros(6,13);
for j=1:size(xdata1,1);
    predictors=find(xdata1(j,:));
    x=data(subsample2,predictors );    
    
   [Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,1,0);
    resmatrix(j,predictors)=Aivx;
    resmatrix(j,13)=Wivx(1,1);
end;
disp('Table 8; Panel B: Monthly Data January 1952-December 2012');
sTable = array2table(resmatrix,'VariableNames',vnam4);
%exclude unnecessary columns from the Table;
sTable8B =horzcat(sTable(:,2),sTable(:,5:8),sTable(:,10),sTable(:,12:end));
disp(sTable8B);


%Table 9: Panel A, Quarterly Data 1927-2012
data=dataq;
y=data(:,13);
xdata=[5 6 0 0 0 0 0; 5 6 10 12 0 0 0; 5 8 0 0 0 0 0; 5 2 0 0 0 0 0; 7 8 12 0 0 0 0; 7 6 11 0 0 0 0];
xdata1=[0 0 0 0 1 1 0 0 0 0 0 0;
        0 0 0 0 1 1 0 0 0 1 0 1; 
        0 0 0 0 1 0 0 1 0 0 0 0;
        0 1 0 0 1 0 0 0 0 0 0 0; 
        0 0 0 0 0 0 1 1 0 0 0 1; 
        0 0 0 0 0 1 1 0 0 0 1 0];


resmatrix=zeros(6,13);
for j=1:size(xdata1,1);
    predictors=find(xdata1(j,:));
    x=data(:,predictors );    
    
   [Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,1,0);
    resmatrix(j,predictors)=Aivx;
    resmatrix(j,13)=Wivx(1,1);
end;
disp('Table 9; Panel A: Quarterly Data Q1 1952-Q4 2012');

sTable = array2table(resmatrix,'VariableNames',vnam4);
%exclude unnecessary columns from the Table;
sTable9A =horzcat(sTable(:,2),sTable(:,5:8),sTable(:,10:end));
disp(sTable9A);



%Table 9: Panel B, Quarterly Data 1952-2012
vnam5={'yyyymm', 	'DE', 	'LTY', 	'DY', 	'DP', 	'TBL', 	'EP', 	'BM', 	'INF', 	'DFY', 	'NTIS', 	'TMS', 	'JointWald', 	'CAY'};
xdata=[5 6 0 0 0 0 0 0; 5 6 10 12 0 0 0 0; 5 8 0 0 0 0 0 0; 5 2 0 0 0 0 0 0; 7 8 12 0 0 0 0 0; 5 2 14 0 0 0 0 0; 7 6 10 14 0 0 0 0];
xdata1=[0 0 0 0 1 1 0 0 0 0 0 0 0 0;
        0 0 0 0 1 1 0 0 0 1 0 1 0 0; 
        0 0 0 0 1 0 0 1 0 0 0 0 0 0;
        0 1 0 0 1 0 0 0 0 0 0 0 0 0; 
        0 0 0 0 0 0 1 1 0 0 0 1 0 0; 
        0 1 0 0 1 0 0 0 0 0 0 0 0 1; 
        0 0 0 0 0 1 1 0 0 1 0 0 0 1];

resmatrix=zeros(6,14);
for j=1:size(xdata1,1);
    predictors=find(xdata1(j,:));
    
     if xdata1(j,14)==1;        
    subsample2=102:345;
    else;
       subsample2=101:345;
    end;
    
    y=data(subsample2,13);
    x=data(subsample2,predictors );    
    
   [Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,1,0);
    resmatrix(j,predictors)=Aivx;
    resmatrix(j,13)=Wivx(1,1);
end;
disp('Table 9; Panel B: Quarterly Data Q1 1952-Q4 2012');
sTable = array2table(resmatrix,'VariableNames',vnam5);
%exclude unnecessary columns from the Table;
sTable9B =horzcat(sTable(:,2),sTable(:,5:8),sTable(:,10),sTable(:,12),sTable(:,14),sTable(:,13));
disp(sTable9B)

%Table 11: Panel A, Monthly Data 1927-2012
vnam6={'Kmths', 	'DE', 	'LTY', 	'DY', 	'DP', 	'TBL', 	'EP', 	'BM', 'DFY', 	'NTIS', 	'TMS', 	'INF'};
data=datam;
Kv=[4,12,24,36,48,60];
y=data(:,13);

xv=[2,3,4,5,6,7,8,10,11,12,9];
resmatrix=zeros(size(Kv,2),11);
for j=1:(size(xv,2));
x=data(:,xv(j));

for i=1:(size(Kv,2));

% K: horizon, l: number of regressors 
K=Kv(1,i);
[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,K,0);
resmatrix(i,j)=Wivx(1,1);

end;

end;
resmatrix=[Kv' resmatrix];
disp('Table 11; Panel A: Monthly Data January 1927-December 2012');
sTable11A = array2table(resmatrix,'VariableNames',vnam6);
disp(sTable11A);


%Table 11: Panel B, Monthly Data 1952-2012
data=datam;
subsample2=301:1033;

y=data(subsample2,13);

xv=[2,3,4,5,6,7,8,10,11,12,9];
resmatrix=zeros(size(Kv,2),11);
for j=1:(size(xv,2));
x=data(subsample2,xv(j));

for i=1:(size(Kv,2));
% K: horizon, l: number of regressors 
K=Kv(1,i);
[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,K,0);
resmatrix(i,j)=Wivx(1,1);

end;

end;
disp('Table 11; Panel B: Monthly Data January 1952-December 2012');
resmatrix=[Kv' resmatrix];
sTable11B = array2table(resmatrix,'VariableNames',vnam6);
disp(sTable11B);





%Table 12: Panel A, Quarterly Data 1927-2012
data=dataq;
Kv=[4,8,12,16,20];
y=data(:,13);
vnam7={'Kqtrs', 'DE', 	'LTY', 	'DY', 	'DP', 	'TBL', 	'EP', 	'BM', 	'DFY', 	'NTIS', 	'TMS', 	'INF'};
xv=[2,3,4,5,6,7,8,10,11,12,9];
resmatrix=zeros(size(Kv,2),11);
for j=1:(size(xv,2));
x=data(:,xv(j));

for i=1:(size(Kv,2));

% K: horizon, l: number of regressors 
K=Kv(1,i);
[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,K,0);
resmatrix(i,j)=Wivx(1,1);

end;

end;
disp('Table 12; Panel A: Quarterly Data Q1 1927-Q4 2012');
resmatrix=[Kv' resmatrix];
sTable12A = array2table(resmatrix,'VariableNames',vnam7);
disp(sTable12A);



%Table 12: Panel B, Quarterly Data 1927-2012
vnam8={'Kqtrs', 'DE', 	'LTY', 	'DY', 	'DP', 	'TBL', 	'EP', 	'BM', 	'DFY', 	'NTIS', 	'TMS', 	'INF', 	'CAY'};
data=dataq;
Kv=[4,8,12,16,20];
xv=[2,3,4,5,6,7,8,10,11,12,9,14];
resmatrix=zeros(size(Kv,2),11);
for j=1:(size(xv,2));
    if xv(j)==14;        
    subsample2=102:345;
    else;
       subsample2=101:345;
    end;
    y=data(subsample2,13);
    x=data(subsample2,xv(j));

for i=1:(size(Kv,2));

% K: horizon, l: number of regressors 
K=Kv(1,i);

[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,K,0);

resmatrix(i,j)=Wivx(1,1);

end;

end;
disp('Table 12; Panel B: Quarterly Data Q1 1952-Q4 2012');
resmatrix=[Kv' resmatrix];
sTable12B = array2table(resmatrix,'VariableNames',vnam8);
disp(sTable12B);


%Table 13: 1927-2012
%Monthly Data
vnam9={'Kqtrs', 'EP', 'TBL', 'JointWald'};

data=datam;
Kv=[4,12,24,36,48,60];

x=[data(:,7) data(:,6) ];
y=data(:,13);
resmatrix=zeros(size(Kv,2),size(x,2)+1);

for i=1:(size(Kv,2));

% K: horizon, l: number of regressors 
K=Kv(1,i);

[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,K,0);
resmatrix(i,:)=[WivxInd(1,:).^2 Wivx(1)];

end;
disp('Table 13; Panel A: Monthly Data');
disp('January 1927-December 2012');
resmatrix=[Kv' resmatrix];
sTable13A = array2table(resmatrix,'VariableNames',vnam9);
disp(sTable13A);




%Table 13: 1952-2012
subsample2=301:1033;
x=[data(subsample2,7) data(subsample2,6) ];

y=data(subsample2,13);

resmatrix=zeros(size(Kv,2),size(x,2)+1);

for i=1:(size(Kv,2));

% K: horizon, l: number of regressors 
K=Kv(1,i);
[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,K,0);
resmatrix(i,:)=[WivxInd(1,:).^2 Wivx(1)];
end;
disp('Table 13; Panel A: Monthly Data');
disp('January 1952-December 2012');
resmatrix=[Kv' resmatrix];
sTable13B = array2table(resmatrix,'VariableNames',vnam9);
disp(sTable13B);



%Table 13: 1927-2012
%Panel B
%Quarterly Data

vnam10={'Kqtrs', 'EP', 'TBL', 'NTIS','JointWald'};
data=dataq;
Kv=[4,8,12,16,20];
%subsample2=301:1033;
x=[data(:,7) data(:,6) data(:,11)];

y=data(:,13);

resmatrix=zeros(size(Kv,2),size(x,2)+1);

for i=1:(size(Kv,2));

% K: horizon, l: number of regressors 
K=Kv(1,i);
[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,K,0);
resmatrix(i,:)=[WivxInd(1,:).^2 Wivx(1)];
end;
disp('Table 13; Panel B: Quarterly Data');
disp('Q1 1927-Q4 2012');
resmatrix=[Kv' resmatrix];
sTable13C = array2table(resmatrix,'VariableNames',vnam10);
disp(sTable13C);



%Table 13: 1952-2012
%Panel B
vnam10={'Kqtrs', 'EP', 'TBL', 'NTIS','CAY','JointWald'};
subsample2=102:345;
x=[data(subsample2,7) data(subsample2,6) data(subsample2,10) data(subsample2,14)];

y=data(subsample2,13);

resmatrix=zeros(size(Kv,2),size(x,2)+1);

for i=1:(size(Kv,2));

% K: horizon, l: number of regressors 
K=Kv(1,i);

[Aols,Aivx,Wivx,WivxInd,Q,corr2]=ivxlh(y,x,K,0);

resmatrix(i,:)=[WivxInd(1,:).^2 Wivx(1)];
%resmatrix(i,:)=[WivxInd(1,:).^2 Wivx(1)];

end;
disp('Table 13; Panel B: Quarterly Data');
disp('Q1 1952-Q4 2012');
resmatrix=[Kv' resmatrix];
sTable13D = array2table(resmatrix,'VariableNames',vnam10);
disp(sTable13D);



