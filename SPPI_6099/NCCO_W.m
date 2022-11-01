
clc; %Clear command window
clear all;
t0 = clock; %Current date and time as date vector


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




%cc
Dist = graphallshortestpaths(SparseAdjMatrix); 
Dist(Dist==Inf) = number;
cc = (number-1)./sum(Dist,2); 
degree = sum(AdjMatrix,2); 


%Establishment of ecc matrix
ECC= zeros(number,number);  

for i = 1:number
    for j =1:number
     trianglescore=0;
     if ( AdjMatrix(i,j) ~=0)
      for k=1:number
       if(k~=j)&&(k~=i)&&(AdjMatrix(i,k)~=0)&&(AdjMatrix(j,k)~=0) %是否形成三角形
           
           trianglescore=trianglescore+cc(i,1)+cc(j,1)+cc(k,1);
       
         
       end
      end
      if((degree(i)>1) && (degree(j)>1))
         
          ECC(i,j)=trianglescore/min(sum(AdjMatrix(i,:)) ,sum(AdjMatrix(j,:)));
      end
     end 
    end
 end
NCC= sum(ECC,2); 




orthology= zeros(number,1); 
OS=zeros(number,1);
lengthprotein=length(protein);
for j=1:lengthprotein
    for i=j:number   
        if strcmp(protein{j,2},Pro{i})
            OS(i,1)=str2num(protein{j,3});
            orthology(i,1)=str2num(protein{j,3}); %字符转换成数字
            continue;
        end
    end    
end


s=1.6;
PF=zeros(number,11);

for k=0:10
   for i=1:number
      for j=i+1:number
           if NCC(i,1)>NCC(j,1)-s && OS(i,1)>OS(j,1)-s
               PF(i,k+1)=PF(i,k+1)+(NCC(i,1)-NCC(j,1))*k/10+(OS(i,1)-OS(j,1))*(1-k/10);
               
           end
           if NCC(i,1)-s<NCC(j,1) && OS(i,1)-s<OS(j,1)
               PF(j,k+1)=PF(j,k+1)+(NCC(j,1)-NCC(i,1))*k/10+(OS(j,1)-OS(i,1))*(1-k/10);
               
           end
      end     
  end
end

% %results
ESpro = importdata('Essential_yeast.txt');

filter =[100,200,300,400,500,600];
ESnum = zeros(6,11);
value1 = zeros(number,11); 
inx1 = zeros(number,11);


for i=1:11
  [value1(:,i),inx1(:,i)] = sort(PF(:,i),'descend'); %按分值序按降排列
end


for j=1:11
  for i = 1:6
     
     ESnum(i,j) = sum(ismember(Pro(inx1(1:filter(i),j)),ESpro)); 
     
  end
  
end

%Model analysis ;The file path needs to be changed based on the storage location
Top100_TP=importdata('Model analysis-NCCO_W/SPPI6099_gene.csv');
weak_Topall=Pro(inx1(1:number,2));
for i=1:number
    NCCO_W_rank(i,1)=i;
end
list=table(weak_Topall,NCCO_W_rank);

weak_Top100=Pro(inx1(1:100,2));
weak_Top100_TP=intersect(weak_Top100,ESpro);
intersection=intersect(weak_Top100_TP,Top100_TP);
weak_Top100_TP_unique=setdiff(weak_Top100_TP,intersection);

Top100_ne=importdata('Model analysis-NCCO_W/SPPI6099_no_essentialgene.csv');
weak_Top100_ne=setdiff(weak_Top100,weak_Top100_TP);
ne_intersection=intersect(Top100_ne,weak_Top100_ne);
Top100_ne_unique=setdiff(Top100_ne,ne_intersection);

Union=union(Top100_ne_unique,weak_Top100_TP_unique);
list_name=table2cell(list(:,1));
[bool1,edge5(:,1)] = ismember(Union,list_name);  
list2=list(edge5(:,1),:);

NCCO_rank=cell(height(list2),1);
EPPI_rank=readtable('Model analysis-NCCO_W/SPPI6099_rank.csv');
EPPI_rank_name=table2cell(EPPI_rank(:,1))
[bool4,edge4(:,1)] = ismember(Union,EPPI_rank_name);
for i=1:length(edge4)
    NCCO_rank{i,1}=EPPI_rank{edge4(i,1),2};
end
list2.NCCO_rank=NCCO_rank;



genename=table2cell(list2(:,1))
label=cell(height(list2),1)

[bool2,edge2(:,1)] = ismember(Top100_ne_unique,genename);
for i=1:length(edge2)
    label{edge2(i,1),1}='NCCO-FP';
end

[bool3,edge3(:,1)] = ismember(weak_Top100_TP_unique,genename);
for i=1:length(edge3)
    label{edge3(i,1),1}='NCCO_W-TP';
end
list2.label=label
list2=sortrows(list2,2)
writetable(list2, 'Model analysis-NCCO_W/SPPI6099_result.csv');



    

