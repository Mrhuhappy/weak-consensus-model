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


%results
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
writecell(Top100_TP, 'Model analysis-PeC_W/EPPI_gene.csv');
Top100_ne=setdiff(Top100,Top100_TP);
writecell(Top100_ne, 'Model analysis-PeC_W/EPPI_no_essentialgene.csv');
for i=1:number
    rank(i,1)=i;
end

Topall=Pro(inx2(1:number,1));
Topall_rank=table(Topall,rank);
writetable(Topall_rank, 'Model analysis-PeC_W/EPPI_rank.csv');

