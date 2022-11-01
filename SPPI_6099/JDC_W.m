clc; %Clear command window
clear all;



%Import data
PIN='krogan2006_extended.txt';
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


meanvalue=zeros(number,1);
meanvalue=mean(Matrix,2);
standarddeviation=zeros(number,1);  
standarddeviation=std(Matrix,1,2); 
variancedeviation=zeros(number,1);  
variancedeviation=var(Matrix,1,2);

for i=1:number
  
  BT(i,1)=1/(1+variancedeviation(i,1));
  % threshold parameter
  %G=zeros(number,1);
  G(i,1)=meanvalue(i,1)+2*standarddeviation(i,1)*BT(i,1);
end
Matrix1=zeros(number,36);
for i=1:number
    for j=1:36
        if  Matrix(i,j)>G(i,1)
            Matrix1(i,j)=1;
        end
    
    end
end


%Create the adjacency matrix
[bool,edge(:,1)] = ismember(Df_r{1},Pro); 
[bool,edge(:,2)] = ismember(Df_r{2},Pro);
AdjMatrix = zeros(number); 
for i = 1:length(Df_r{1}) 
    AdjMatrix(edge(i,1),edge(i,2)) = 1;
    AdjMatrix(edge(i,2),edge(i,1)) = 1;
end
SparseAdjMatrix = sparse(AdjMatrix); 
degree = sum(AdjMatrix,2); 



%Establishment of ecc matrix
ECC= zeros(number,number);  

for i = 1:number
    for j = 1:number
        
        
        D=[Matrix1(i,:);Matrix1(j,:)];
        Jaccard(i,j)=1-pdist(D,'Jaccard');
       
       neighbornumber=0;
       if ( AdjMatrix(i,j) ~=0)
          for k=1:number
            if(k~=j)&&(k~=i)&&(AdjMatrix(i,k)~=0)&&(AdjMatrix(j,k)~=0)
             neighbornumber=neighbornumber+1;
            end
          end
      
       if((degree(i)>1) && (degree(j)>1))
       
        ECC(i,j)=neighbornumber/min(sum(AdjMatrix(i,:))-1 ,sum(AdjMatrix(j,:))-1);
        end
    
 end
    
      
    end
end
 
ECC(isnan(ECC)==1)=0;
ecc=sum(ECC,2);  
Jaccard(isnan(Jaccard)==1)=0;
JD=sum(Jaccard,2);
 s=2;
walkPF=zeros(number,11);

for k=0:10
    for i=1:number
       for j=1+i:number

         if JD(i)>JD(j)-s && ecc(i)>ecc(j)-s
        
         walkPF(i,k+1)=walkPF(i,k+1)+(JD(i)-JD(j))*k/10+(ecc(i)-ecc(j))*(1-k/10);
     
        
         
           end
           if  JD(j)>JD(i)-s && ecc(j)>ecc(i)-s
           walkPF(j,k+1)=walkPF(j,k+1)+(JD(j)-JD(i))*k/10+(ecc(j)-ecc(i))*(1-k/10);
          
          
           end
          
        
       end   
    end
end


%results
ESpro = importdata('Essential_yeast.txt');
filter =[100,200,300,400,500,600]
value2= zeros(number,1);
inx2=zeros(number,1);
for i=1:11
    [value2(:,i),inx2(:,i)]=sort(walkPF(:,i),'descend');
end   


walkESnum = zeros(6,1);

ESnumten = zeros(6,1);

for j=1:11
  for i = 1:6
    
     walkESnum(i,j) = sum(ismember(Pro(inx2(1:filter(i),j)),ESpro)); 
      
  end
end


%Model analysis ;The file path needs to be changed based on the storage location
Top100_TP=importdata('Model analysis-JDC_W/SPPI6099_gene.csv');
weak_Topall=Pro(inx2(1:number,1));
for i=1:number
    JDC_W_rank(i,1)=i;
end
list=table(weak_Topall,JDC_W_rank);
weak_Top100=Pro(inx2(1:100,1));
weak_Top100_TP=intersect(weak_Top100,ESpro);
intersection=intersect(weak_Top100_TP,Top100_TP);
weak_Top100_TP_unique=setdiff(weak_Top100_TP,intersection);
Top100_ne=importdata('Model analysis-JDC_W/SPPI6099_no_essentialgene.csv');
weak_Top100_ne=setdiff(weak_Top100,weak_Top100_TP);
ne_intersection=intersect(Top100_ne,weak_Top100_ne);
Top100_ne_unique=setdiff(Top100_ne,ne_intersection);
Union=union(Top100_ne_unique,weak_Top100_TP_unique);
list_name=table2cell(list(:,1));
[bool1,edge5(:,1)] = ismember(Union,list_name);  
list2=list(edge5(:,1),:);
JDC_rank=cell(height(list2),1);
EPPI_rank=readtable('Model analysis-JDC_W/SPPI6099_rank.csv');
EPPI_rank_name=table2cell(EPPI_rank(:,1))
[bool4,edge4(:,1)] = ismember(Union,EPPI_rank_name);
for i=1:length(edge4)
    JDC_rank{i,1}=EPPI_rank{edge4(i,1),2};
end
list2.JDC_rank=JDC_rank;
genename=table2cell(list2(:,1))
label=cell(height(list2),1)
[bool2,edge2(:,1)] = ismember(Top100_ne_unique,genename);
for i=1:length(edge2)
    label{edge2(i,1),1}='JDC-FP';
end
[bool3,edge3(:,1)] = ismember(weak_Top100_TP_unique,genename);
for i=1:length(edge3)
    label{edge3(i,1),1}='JDC_W-TP';
end
list2.label=label
list2=sortrows(list2,2)
writetable(list2, 'Model analysis-JDC_W/SPPI6099_result.csv');
