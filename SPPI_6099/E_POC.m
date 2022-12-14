clc; %Clear command window
clear all;


%Import data
allorthology='Nucleus_protein_allorthology.txt';
DATAfile =allorthology;
fid = fopen(DATAfile,'r'); %open the file 
Df = textscan(fid,'%s%s%d%s'); %read file
fclose(fid); %close file, success 0, fail -1 
protein=cat(2,Df{1},Df{2},Df{4});
essentialproteinnumbers = sum(Df{3});  
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
[bool,edge(:,1)] = ismember(Df_r{1},Pro); 
[bool,edge(:,2)] = ismember(Df_r{2},Pro);
AdjMatrix = zeros(number); 
for i = 1:length(Df_r{1})
    AdjMatrix(edge(i,1),edge(i,2)) = 1;
    AdjMatrix(edge(i,2),edge(i,1)) = 1;
end
SparseAdjMatrix = sparse(AdjMatrix); 

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
degree = sum(AdjMatrix,2); 



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
      ECC(i,j)=neighbornumber/min(sum(AdjMatrix(i,:))-1 ,sum(AdjMatrix(j,:))-1);
      end
     end 
    end
   end
SOECC = sum(ECC,2); 


%calculate homology scores were calculated
orthology= zeros(number,1); 
lengthprotein=length(protein);
for j=1:lengthprotein
    for i=1:number   
        if strcmp(protein{j,2},Pro{i})
            
            orthology(i,1)=str2num(protein{j,3})/99; %??????????????
            continue;
        end
    end    
end



%Compare the results
ESpro = importdata('Essential_yeast.txt');
filter =[100,200,300,400,500,600];

%Combined with the model
RWPOCPF=zeros(number,11);
for k=0:10
   for i=1:number
      for j=i+1:number
           if tp1(i,1)>tp1(j,1) && orthology(i,1)>orthology(j,1)
              RWPOCPF(i,k+1)=RWPOCPF(i,k+1)+(tp1(i,1)-tp1(j,1))*k/10+(orthology(i,1)-orthology(j,1))*(1-k/10);
           end
           if tp1(i,1)<tp1(j,1) && orthology(i,1)<orthology(j,1)
              RWPOCPF(j,k+1)=RWPOCPF(j,k+1)+(tp1(j,1)-tp1(i,1))*k/10+(orthology(j,1)-orthology(i,1))*(1-k/10);
           end
      end     
  end
end
value1 = zeros(number,11); 
inx1 = zeros(number,11);
RWPOCESnum = zeros(6,11);

for i=1:11
  [value1(:,i),inx1(:,i)] = sort(RWPOCPF(:,i),'descend'); 
end

for j=1:11
  for i = 1:6
     
     RWPOCESnum(i,j) = sum(ismember(Pro(inx1(1:filter(i),j)),ESpro)); 
     
  end
end

%The output file ;The file path needs to be changed based on the storage location
Top100=Pro(inx1(1:100,7));
Top100_TP=intersect(Top100,ESpro);
writecell(Top100_TP, 'Model analysis-E_POC_W/SPPI6099_gene.csv');
Top100_ne=setdiff(Top100,Top100_TP);
writecell(Top100_ne, 'Model analysis-E_POC_W/SPPI6099_no_essentialgene.csv');
for i=1:number
    rank(i,1)=i;
end
Topall=Pro(inx1(1:number,7));
Topall_rank=table(Topall,rank);
writetable(Topall_rank, 'Model analysis-E_POC_W/SPPI6099_rank.csv');



