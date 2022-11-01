
clc; %Clear command window
clear all;

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


%calculate homology scores were calculated
OS= zeros(number,1); 
lengthprotein=length(protein);
for j=1:lengthprotein
    for i=1:number   
        if strcmp(protein{j,2},Pro{i})
            OS(i,1)=str2num(protein{j,3})/99; 
            continue;
        end
    end    
end


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
tp1=SOECC/max(SOECC);


s=4;
 walkPF=zeros(number,11);
for k=0:10
 for i=1:number
     for j=i+1:number  
         if tp1(i,1)>tp1(j,1)-s/10&&OS(i,1)>OS(j,1)-s/10
         walkPF(i,k+1)=walkPF(i,k+1)+(tp1(i,1)-tp1(j,1))*k/10+(OS(i,1)-OS(j,1))*(1-k/10);
         end
         
         if tp1(i,1)-s/10<tp1(j,1) && OS(i,1)-s/10<OS(j,1)
         walkPF(j,k+1)=walkPF(j,k+1)+(tp1(j,1)-tp1(i,1))*k/10+(OS(j,1)-OS(i,1))*(1-k/10);
         end
         
        end
    
   end

end

%results
ESpro = importdata('Essential_Ecoli1.txt');
filter =[100,200,300,400,500,600]; 
ESnum = zeros(6,11);
value1 = zeros(number,11); 
inx1 = zeros(number,11);
for i=1:11
  [value1(:,i),inx1(:,i)] = sort(walkPF(:,i),'descend'); 
end
for j=1:11
  for i = 1:6
      
     ESnum(i,j) = sum(ismember(Pro(inx1(1:filter(i),j)),ESpro)); 
     
  end
end

%Model analysis ;The file path needs to be changed based on the storage location
Top100_TP=importdata('Model analysis-ION_W/EPPI_gene.csv');
weak_Topall=Pro(inx1(1:number,7));
for i=1:number
    ION_W_rank(i,1)=i;
end
list=table(weak_Topall,ION_W_rank);

weak_Top100=Pro(inx1(1:100,7));
weak_Top100_TP=intersect(weak_Top100,ESpro);
intersection=intersect(weak_Top100_TP,Top100_TP);
weak_Top100_TP_unique=setdiff(weak_Top100_TP,intersection);

Top100_ne=importdata('Model analysis-ION_W/EPPI_no_essentialgene.csv');
weak_Top100_ne=setdiff(weak_Top100,weak_Top100_TP);
ne_intersection=intersect(Top100_ne,weak_Top100_ne);
Top100_ne_unique=setdiff(Top100_ne,ne_intersection);

Union=union(Top100_ne_unique,weak_Top100_TP_unique);
list_name=table2cell(list(:,1));
[bool1,edge5(:,1)] = ismember(Union,list_name);  
list2=list(edge5(:,1),:);
ION_rank=cell(height(list2),1);
EPPI_rank=readtable('Model analysis-ION_W/EPPI_rank.csv');
EPPI_rank_name=table2cell(EPPI_rank(:,1))
[bool4,edge4(:,1)] = ismember(Union,EPPI_rank_name);
for i=1:length(edge4)
    ION_rank{i,1}=EPPI_rank{edge4(i,1),2};
end
list2.ION_rank=ION_rank;
genename=table2cell(list2(:,1))
label=cell(height(list2),1)
[bool2,edge2(:,1)] = ismember(Top100_ne_unique,genename);
for i=1:length(edge2)
    label{edge2(i,1),1}='ION-FP';
end

[bool3,edge3(:,1)] = ismember(weak_Top100_TP_unique,genename);
for i=1:length(edge3)
    label{edge3(i,1),1}='ION_W-TP';
end
list2.label=label
list2=sortrows(list2,2)
writetable(list2, 'Model analysis-ION_W/EPPI_result.csv');
