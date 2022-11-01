clc; %Clear command window
clear all;



%Import data
PIN='Ecoli_DIP_ECC1.txt';
DATAfile_r=PIN;
fid_r = fopen(DATAfile_r,'r'); %open file
Df_r = textscan(fid_r,'%s%s%s'); %read file
fclose(fid_r); %close file success 0,fail -1
Pro = union(Df_r{1},Df_r{2});
number = length(Pro); 
IN='Ecoli_gene expression data1.txt';
fin = fopen(IN,'r'); %open file
Din = textscan(fin,'%s%f%f%f%f%f%f%f%f%f'); %read file
fclose(fin); %close file success 0,fail -1
express=cat(2,Din{2},Din{3},Din{4},Din{5},Din{6},Din{7},Din{8},Din{9},Din{10});

% to build genedata matrix
[bool1,edge1(:,1)] = ismember(Pro,Din{1});
Matrix=zeros(number,9);
for i=1:number
    if bool1(i,1)==1
       for j =1:9
            Matrix(i,j)=express(edge1(i,1),j);  % genedata matrix
       end
    end
end

standarddeviation=zeros(number,1);
standarddeviation=std(Matrix,0,2);
meanvalue=zeros(number,1);
meanvalue=mean(Matrix,2);

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
%PCC matrix
PCCMatrix=zeros(number);
ppc=zeros(length(Df_r{1}),1);
for i = 1:length(Df_r{1}) 
    for j=1:9
        if standarddeviation(edge(i,1),1) ~=0 & standarddeviation(edge(i,2),1) ~=0
            ppc(i,1) =  ppc(i,1)+(( Matrix(edge(i,1),j)-meanvalue(edge(i,1),1))/standarddeviation(edge(i,1),1) * ( Matrix(edge(i,2),j)-meanvalue(edge(i,2),1))/standarddeviation(edge(i,2),1)) ;
        end
    end 
    PCCMatrix(edge(i,1),edge(i,2)) = ppc(i,1)/8; 
    PCCMatrix(edge(i,2),edge(i,1)) = ppc(i,1)/8;
end

sumpcc=sum(PCCMatrix,2);
Sumpcc=atan(sumpcc)*2/pi;


%Establishment of ecc matrix
ECC= zeros(number,number);  

for i = 1:number
    for j = 1:number
    neighbornumber=0;
     if ( AdjMatrix(i,j) ~=0)
      for k=1:number
       if(k~=j)&&(k~=i)&&(AdjMatrix(i,k)~=0)&&(AdjMatrix(j,k)~=0)
       neighbornumber=neighbornumber+1;
       end
      end
      if((degree(i)>1) && (degree(j)>1))
      ECC(i,j)=neighbornumber/min(sum(AdjMatrix(i,:)) ,sum(AdjMatrix(j,:)));
      end
     end 
    end
 end
SOECC= sum(ECC,2); 
ecc=SOECC/max(SOECC);



%Compare the results
ESpro = importdata('Essential_Ecoli1.txt');
filter =[100,200,300,400,500,600]; 

count=zeros(number,11)
 s=3;
walkPF=zeros(number,11);

for k=0:10
    for i=1:number
       for j=1+i:number
           
           if Sumpcc(i)>Sumpcc(j)-s/10 && ecc(i)>ecc(j)-s/10
         
        walkPF(i,k+1)=walkPF(i,k+1)+(Sumpcc(i)-Sumpcc(j))*k/10+(ecc(i)-ecc(j))*(1-k/10);
        count(i,k+1)=count(i,k+1)+1;
         
          
           end
           if  Sumpcc(j)>Sumpcc(i)-s/10 && ecc(j)>ecc(i)-s/10
          
       walkPF(j,k+1)=walkPF(j,k+1)+(Sumpcc(j)-Sumpcc(i))*k/10+(ecc(j)-ecc(i))*(1-k/10);
           count(j,k+1)=count(j,k+1)+1;

           end
          
        
       end   
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

Top100_TP=importdata('Model analysis-PeC_W/EPPI_gene.csv');
weak_Topall=Pro(inx2(1:number,2));
for i=1:number
    PeC_W_rank(i,1)=i;
end
list=table(weak_Topall,PeC_W_rank);
weak_Top100=Pro(inx2(1:100,2));
weak_Top100_TP=intersect(weak_Top100,ESpro);
intersection=intersect(weak_Top100_TP,Top100_TP);
weak_Top100_TP_unique=setdiff(weak_Top100_TP,intersection);
Top100_ne=importdata('Model analysis-PeC_W/EPPI_no_essentialgene.csv');
weak_Top100_ne=setdiff(weak_Top100,weak_Top100_TP);
ne_intersection=intersect(Top100_ne,weak_Top100_ne);
Top100_ne_unique=setdiff(Top100_ne,ne_intersection);
Union=union(Top100_ne_unique,weak_Top100_TP_unique);
list_name=table2cell(list(:,1));
[bool1,edge5(:,1)] = ismember(Union,list_name);  
list2=list(edge5(:,1),:);
PeC_rank=cell(height(list2),1);
EPPI_rank=readtable('Model analysis-PeC_W/EPPI_rank.csv');
EPPI_rank_name=table2cell(EPPI_rank(:,1))
[bool4,edge4(:,1)] = ismember(Union,EPPI_rank_name);
for i=1:length(edge4)
    PeC_rank{i,1}=EPPI_rank{edge4(i,1),2};
end
list2.PeC_rank=PeC_rank;
genename=table2cell(list2(:,1))
label=cell(height(list2),1)
[bool2,edge2(:,1)] = ismember(Top100_ne_unique,genename);
for i=1:length(edge2)
    label{edge2(i,1),1}='PeC-FP';
end
[bool3,edge3(:,1)] = ismember(weak_Top100_TP_unique,genename);
for i=1:length(edge3)
    label{edge3(i,1),1}='PeC_W-TP';
end
list2.label=label
list2=sortrows(list2,2)
writetable(list2, 'Model analysis-PeC_W/EPPI_result.csv');
