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


meanvalue=zeros(number,1);
meanvalue=mean(Matrix,2);
standarddeviation=zeros(number,1);  
standarddeviation=std(Matrix,1,2);  
variancedeviation=zeros(number,1);  
variancedeviation=var(Matrix,1,2);
for i=1:number

  BT(i,1)=1/(1+variancedeviation(i,1));
  
end 

% threshold parameter
for i=1:number
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
        
        if  AdjMatrix(i,j)==1 && i~=j
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
end
 
ECC(isnan(ECC)==1)=0;
ecc= sum(ECC,2);
Jaccard(isnan(Jaccard)==1)=0;
JD=sum(Jaccard,2);
jdc=zeros(number,1);
for i=1:number

       jdc(i,1)=JD(i,1)*ecc(i,1)
      
end








%results
ESpro = importdata('Essential_yeast.txt');

filter =[100,200,300,400,500,600]
value2= zeros(number,1);
inx2=zeros(number,1);
[value2(:,1),inx2(:,1)]=sort(jdc(:,1),'descend');
 


walkESnum = zeros(6,1);

ESnumten = zeros(6,1);

  for i = 1:6
    
     walkESnum(i,1) = sum(ismember(Pro(inx2(1:filter(i),1)),ESpro)); 
     
  end
  
%Model analysis ;The file path needs to be changed based on the storage location
Top100=Pro(inx2(1:100,1));
Top100_TP=intersect(Top100,ESpro);
writecell(Top100_TP, 'Model analysis-JDC_W/SPPI6099_gene.csv');

Top100_ne=setdiff(Top100,Top100_TP);
writecell(Top100_ne, 'Model analysis-JDC_W/SPPI6099_no_essentialgene.csv');

for i=1:number
    rank(i,1)=i;
end

Topall=Pro(inx2(1:number,1));
Topall_rank=table(Topall,rank);
writetable(Topall_rank, 'Model analysis-JDC_W/SPPI6099_rank.csv');

  
 