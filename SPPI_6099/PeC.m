clc; %Clear command window
clear all;



%Import data
PIN='filterPPI1.txt';
DATAfile_r=PIN;
fid_r = fopen(DATAfile_r,'r'); %open file
Df_r = textscan(fid_r,'%s%s'); %read file
fclose(fid_r); %close file success 0,fail -1
Pro = union(Df_r{1},Df_r{2}); 
number = length(Pro); 
IN='yeastgenedata.txt';
fin = fopen(IN,'r'); %open file
Din = textscan(fin,'%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'); %read file
fclose(fin); %close file success 0,fail -1
express=cat(2,Din{3},Din{4},Din{5},Din{6},Din{7},Din{8},Din{9},Din{10},Din{11},Din{12},Din{13},Din{14},Din{15},Din{16},Din{17},Din{18},Din{19},Din{20},Din{21},Din{22},Din{23},Din{24},Din{25},Din{26},Din{27},Din{28},Din{29},Din{30},Din{31},Din{32},Din{33},Din{34},Din{35},Din{36},Din{37},Din{38});


% to build genedata matrix
[bool1,edge1(:,1)] = ismember(Pro,Din{2});
Matrix=zeros(number,36);
for i=1:number
    if bool1(i,1)==1
       for j =1:36
            Matrix(i,j)=express(edge1(i,1),j);  % genedata matrix
       end
    end
end

standarddeviation=zeros(number,1);
standarddeviation=std(Matrix,0,2);
meanvalue=zeros(number,1);
meanvalue=mean(Matrix,2);

%Create the adjacency matrix
[bool,edge(:,1)] = ismember(Df_r{1},Pro); %edge以蛋白质编号的形式存放网络中的所有边
[bool,edge(:,2)] = ismember(Df_r{2},Pro);
AdjMatrix = zeros(number); %AdjMatrix为邻接矩阵
for i = 1:length(Df_r{1}) %length(Df_r{1})是网络的边数
    AdjMatrix(edge(i,1),edge(i,2)) = 1;
    AdjMatrix(edge(i,2),edge(i,1)) = 1;
end


%PCC matrix
PCCMatrix=zeros(number);
ppc=zeros(length(Df_r{1}),1);
for i = 1:length(Df_r{1}) 
    for j=1:36
        if standarddeviation(edge(i,1),1) ~=0 & standarddeviation(edge(i,2),1) ~=0
            ppc(i,1) =  ppc(i,1)+(( Matrix(edge(i,1),j)-meanvalue(edge(i,1),1))/standarddeviation(edge(i,1),1) * ( Matrix(edge(i,2),j)-meanvalue(edge(i,2),1))/standarddeviation(edge(i,2),1)) ;
        end
    end 
    PCCMatrix(edge(i,1),edge(i,2)) = ppc(i,1)/35; 
    PCCMatrix(edge(i,2),edge(i,1)) = ppc(i,1)/35;
    
end

sumpcc=sum(PCCMatrix,2);
Sumpcc=atan(sumpcc)*2/pi;


%Establishment of ecc matrix
ECC= zeros(number,number);  %ECC为边聚集系数矩阵

for i = 1:number
    for j = 1:number
    neighbornumber=0;
    if ( AdjMatrix(i,j) ~=0)
      for k=1:number
       if(k~=j)&&(k~=i)&&(AdjMatrix(i,k)~=0)&&(AdjMatrix(j,k)~=0)
       neighbornumber=neighbornumber+1;
       end
      end
     
      ECC(i,j)=neighbornumber/min(sum(AdjMatrix(i,:)) ,sum(AdjMatrix(j,:)));
    
     end 
    end
 end
SOECC= sum(ECC,2); 

ecc=SOECC/max(SOECC);


%Compare the results
ESpro = importdata('Essential_yeast.txt');
%filter =[51,255,510,764,1019,1274];
filter =[100,200,300,400,500,600]; 



walkPF=zeros(number,11);

for i=0:10
    for j=1:number
      
       walkPF(j,i+1)=Sumpcc(j)*ecc(j);    
        
    end
end


value2= zeros(number,11);
inx2=zeros(number,11);
for i=1:11
    [value2(:,i),inx2(:,i)]=sort(walkPF(:,i),'descend');
end   


walkESnum = zeros(6,11);

ESnumten = zeros(6,11);

for j=1:11
  for i = 1:6
    
     walkESnum(i,j) = sum(ismember(Pro(inx2(1:filter(i),j)),ESpro)); 
     
  end
end


%The output file ;The file path needs to be changed based on the storage location
Top100=Pro(inx2(1:100,1));
Top100_TP=intersect(Top100,ESpro);
writecell(Top100_TP, 'Model analysis-PeC_W/SPPI6099_gene.csv');


Top100_ne=setdiff(Top100,Top100_TP);
writecell(Top100_ne, 'Model analysis-PeC_W/SPPI6099_no_essentialgene.csv');

for i=1:number
    rank(i,1)=i;
end

Topall=Pro(inx2(1:number,1));
Topall_rank=table(Topall,rank);
writetable(Topall_rank, 'Model analysis-PeC_W/SPPI6099_rank.csv');




%计算RWPOC的Precision-Recall  
  RWPOCPR=zeros(number,11);
 for n=1:11
    for i1 = 1:number
     
    RWPOCPR(i1,n) = sum(ismember(Pro(inx2(1:i1,n)),ESpro))/i1; 
     
    end
end
  
    
  %计算RWPOC的Jackknife 
  jkRWPOC=zeros(600,11);
  for m=1:11
   for w = 1:600
      jkRWPOC(w,m) = sum(ismember(Pro(inx2(1:w,m)),ESpro)); 
   end
  end
  
 
  
  
  
%计算程序运行时间
%elapsed_time = etime(clock,t0);
%save dataDC/elapsed_time.mat elapsed_time;

%输出到文件		
save dataWDC/RWPOCPR.mat RWPOCPR;
save dataWDC/jkRWPOC.matjkRWPOC;
save dataWDC/walkESnum.mat walkESnum;

%OUT=strcat('locs/',locs,'_ION.mat'); 
%save (OUT,'results')
%save Data/OUT results '-ascii' '-tab';

