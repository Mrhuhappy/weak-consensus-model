
clc; %Clear command window
clear all;
t0 = clock; %Current date and time as date vector


%Import data
allorthology='protein_allorthology1.txt';
DATAfile =allorthology;
fid = fopen(DATAfile,'r'); %open the file 
Df = textscan(fid,'%s%s%d%s'); %read file
fclose(fid); %close file, success 0, fail -1 
protein=cat(2,Df{1},Df{2},Df{4});
essentialproteinnumbers = sum(Df{3});  
PIN='Ecoli_DIP_ECC1.txt';
DATAfile_r=PIN;
fid_r = fopen(DATAfile_r,'r'); %open file
Df_r = textscan(fid_r,'%s%s%s'); %read file
fclose(fid_r); %close file success 0,fail -1
Pro = union(Df_r{1},Df_r{2}); 
Pro=union(Pro,Df{2});
number = length(Pro);
PF=zeros(number,1);

%Create the adjacency matrix
[bool,edge(:,1)] = ismember(Df_r{1},Pro); 
[bool,edge(:,2)] = ismember(Df_r{2},Pro);
AdjMatrix = zeros(number); 
for i = 1:length(Df_r{1}) 
    AdjMatrix(edge(i,1),edge(i,2)) = 1;
    AdjMatrix(edge(i,2),edge(i,1)) = 1;
end
degree = sum(AdjMatrix,2); 
SparseAdjMatrix = sparse(AdjMatrix);
ratio=99; 
%Random walk parameter
alpha=1;
%calculate homology scores were calculated
PR_0= zeros(number,1); 
lengthprotein=length(protein);
for j=1:lengthprotein
    for i=j:number   
        if strcmp(protein{j,2},Pro{i})
            PR_0(i,1)=str2num(protein{j,3})/99; 
            continue;
        end
    end    
end

ECC= zeros(number,number); 
%Establishment of ecc matrix
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
HECC=zeros(number,1);
for i=1:number
    for j=1:number
       if SOECC(i)>0
            HECC(i,j)=ECC(i,j)/SOECC(i);
       end
    end
end
tp1=zeros(number,1);
tp1=sum(HECC,2);
H_0= tp1;  %H_0 is transfer matrix
tNodes= zeros(number,1);
for i=1:number
tNodes(i)=(1-alpha)*PR_0(i);
 
end
H_0=alpha*H_0;

%iteration until convergence 
%if error is less than the value of epsilon£¬the values of pagerank are got
k=0;
epsilon=0.00000001;
error=1;
while(error >=epsilon)  
      tmpPR=PR_0; 
      PR_0=tNodes + H_0.* PR_0 ;
      error=0;
     for i = 1:number
         error=error+abs(PR_0(i)-tmpPR(i));
     end
     k=k+1;
end   

%results
ESpro = importdata('Essential_Ecoli1.txt');

filter =[100,200,300,400,500,600]; 
ESnum = zeros(6,1);
value1 = zeros(number,1); 
inx1 = zeros(number,1);
[value1(:,1),inx1(:,1)] = sort(PR_0(:,1),'descend'); 
  for i = 1:6
      a1=inx1(1:filter(i),1);
      a2=Pro(a1);
     ESnum(i,1) = sum(ismember(a2,ESpro));     
  end
%The output file;;The file path needs to be changed based on the storage location
Top100=Pro(inx1(1:100,1));
Top100_TP=intersect(Top100,ESpro);
writecell(Top100_TP, 'Model analysis-ION_W/EPPI_gene.csv');
Top100_ne=setdiff(Top100,Top100_TP);
writecell(Top100_ne, 'Model analysis-ION_W/EPPI_no_essentialgene.csv');
for i=1:number
    rank(i,1)=i;
end
Topall=Pro(inx1(1:number,1));
Topall_rank=table(Topall,rank);
writetable(Topall_rank, 'Model analysis-ION_W/EPPI_rank.csv');
