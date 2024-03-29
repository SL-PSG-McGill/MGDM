%{ 
 Author: Shanglin Li, Photonic DataCom Team, McGill University 
 Date: April 2023
 Note:
    The script uses The Genetic Algorithm Toolbox for MATLAB developed by
	The University of Sheffield (http://uos-codem.github.io/GA-Toolbox/).
%}
clear variables
addpath(genpath(pwd)) 
min_channel_power=0.5; % the power coupling coefficient of the signal in each channel is required to be not less than 0.5
max_channel_XT=0.2;
Dataset_No = input("Which case is chosen? Please input 1 or 2.Input 0 to exit the program.\n");
while ~ismember(Dataset_No, [1,2,0])
    Dataset_No =input("Your input is WRONG! Please input 1 or 2.Input 0 to exit the program.\n");
end
if Dataset_No==1
    Dataset_name='Case_I_dataset';
elseif Dataset_No==2
    Dataset_name='Case_II_dataset';
elseif Dataset_No==0
    fprintf(2,'Exit the program!\n');
    return;
end
load(Dataset_name)
power_beam_PMN=power_beam_PMN';
[PMN_num,beam_info_num]=size(power_beam_PMN);
tx_num = input("How many fiber channels are excited? Please input 5 or 6.Input 0 to exit the program.\n");
while ~ismember(tx_num, [5,6,0])
    tx_num =input("Your input is WRONG! Please input 5 or 6. Input 0 to exit the program. \n");
end
split_PMN=[];
loop1=1;
if tx_num==5
	 for sec1=1:1:PMN_num-tx_num+1
	  for  sec2=1:1:PMN_num-sec1-tx_num+2
		  for sec3=1:1:PMN_num-sec1-sec2-tx_num+3
		   for sec4=1:1:PMN_num-sec1-sec2-sec3-tx_num+4
			   sec5=PMN_num-sec1-sec2-sec3-sec4;
			   split_PMN(loop1,1:tx_num)=[sec1,sec2,sec3,sec4,sec5];
			   loop1=loop1+1;
		   end
		  end
	  end
	end
elseif tx_num==6
    for sec1=1:1:PMN_num-tx_num+1
        for  sec2=1:1:PMN_num-sec1-tx_num+2
            for sec3=1:1:PMN_num-sec1-sec2-tx_num+3
                for sec4=1:1:PMN_num-sec1-sec2-sec3-tx_num+4
                    for sec5=1:1:PMN_num-sec1-sec2-sec3-sec4-tx_num+5
                        sec6=PMN_num-sec1-sec2-sec3-sec4-sec5;
                        split_PMN(loop1,1:tx_num)=[sec1,sec2,sec3,sec4,sec5,sec6];
                        loop1=loop1+1;
                    end
                end
            end
        end
    end
elseif tx_num==0
    fprintf(2,'Exit the program!\n');
    return;
end
[split_num,~]=size(split_PMN);

%% Multi-group Genetic Algorithm
f = waitbar(0,'The calculation of launch conditions starts! Please Wait...');
NIND=40000;
MAXGEN=200; % Max numbers of the evolution generation
NVAR=tx_num+1; % The number of variables, where '1' denotes how to split the mode groups
GGAP=0.8; % Generation gap
XOVR=1; % Crossover rate
INSR=0.9; % Insertion rate
SUBPOP=NIND/400; % the number of sub-populations
MIGR=0.4; % Migration rate
MIGGEN=2; % Migration occurs for every MIGGEN generations
MUTR=0.2; % Mutation rate
BaseV=crtbase([tx_num,1],[beam_info_num,split_num]);
[Chrom,Lind]=crtbp(NIND,BaseV);% Initial population
Chrom_addone=Chrom+ones(NIND,tx_num+1);
gen_count=0;
ObjV=XT_all(power_beam_PMN,split_PMN,Chrom_addone,tx_num,min_channel_power);
[Best_ObjV,Index]=min(ObjV);
Optimal_Chrom=Chrom(Index,:);
while gen_count<MAXGEN
    FitnV=ranking(ObjV,[2,1],SUBPOP); % Calculate fitness scores
    Selch1=select('sus',Chrom,FitnV,GGAP,SUBPOP); % Selection
    Selch2=recombin('xovsh',Selch1,XOVR,SUBPOP); % Recombination
    Selch3=mutate('mut',Selch2,BaseV,MUTR,SUBPOP);% Mutation
    ObjV_Sel=XT_all(power_beam_PMN,split_PMN,Selch3+1,tx_num,min_channel_power); % Calculate the overall inter-channel crosstalk
    [Chrom, ObjV]=reins(Chrom,Selch3,SUBPOP,[1,INSR],ObjV,ObjV_Sel); % Re-insertion
    gen_count=gen_count+1;
    if (rem(gen_count,MIGGEN)==0)
        [Chrom, ObjV]=migrate(Chrom,SUBPOP,[MIGR,1,1],ObjV); % Migration
    end
     [Best_ObjV_new,Index]=min(ObjV);
    if Best_ObjV_new<=Best_ObjV
        Best_ObjV=Best_ObjV_new;
        Optimal_Chrom=Chrom(Index,:);
    else
        Chrom(1,:)=Optimal_Chrom;
        ObjV(1,:)=Best_ObjV;
    end       
    waitbar(gen_count/MAXGEN,f,['Running ',num2str(gen_count/MAXGEN*100,'%.1f'),'%'])
end
waitbar(1,f,'The calculation of launch conditions is done!')
Best_Ind_GA=Optimal_Chrom+1;
% Output mode-group distribution in each fiber channel
Best_PMN_Dist_GA=split_PMN(Best_Ind_GA(end),:); 
MG_start=zeros(tx_num,1); MG_end=zeros(tx_num,1); Best_MG_Channel=strings(1,tx_num);
for loop2=1:tx_num
    if loop2==1
        MG_start(loop2)=1;
    else
        MG_start(loop2)=1+MG_end(loop2-1);
    end
    MG_end(loop2)=MG_start(loop2)+Best_PMN_Dist_GA(loop2)-1;
    Best_MG_Channel(loop2)=[num2str(MG_start(loop2)),'-',num2str(MG_end(loop2))];
end
disp('The mode group distributions are:'); disp(Best_MG_Channel);
% Output selected launch conditions
Best_beam_info_GA=beam_info(Best_Ind_GA(1:tx_num)',:);
Best_beam_info_result(:,1:2)=Best_beam_info_GA(:,1:2)*1e6;
Best_beam_info_result(:,3:4)=Best_beam_info_GA(:,3:4)*180/pi;
fprintf("The launch conditions to excite fiber channels are\nspot size(micron)/radial offset(micron)/incidence tilt(degree)/azimuthal tilt(degree):\n");
disp(Best_beam_info_result);
% Output the power coupling coefficient matrix
Best_powerPMN_tx_GA=power_beam_PMN(:,Best_Ind_GA(1:tx_num));
Best_powerPMN_tx_cell_GA=mat2cell(Best_powerPMN_tx_GA,Best_PMN_Dist_GA);
Best_block_sum_GA=cell2mat(cellfun(@(x)sum(x,1),Best_powerPMN_tx_cell_GA,'UniformOutput',false));
Best_block_sum_db_GA=10*log10(Best_block_sum_GA);
disp('The power coupling coefficient matrix is (in dB):'); disp(Best_block_sum_db_GA);

%% sub-function
function Eval=XT(powerPMN_tx,split_PMN_ind,min_channel_power) % Target Function
powerPMN_tx_cell=mat2cell(powerPMN_tx,split_PMN_ind);
block_sum=cell2mat(cellfun(@(x)sum(x,1),powerPMN_tx_cell,'UniformOutput',false));
target_power=diag(block_sum);
if all(target_power>=min_channel_power)  % the requirements of the power coupling coefficient matrices
    Eval=sum(-10*log10(target_power+1e-60));
else
    Eval=sum(-10*log10(target_power+1e-60)+1e10); % 1e10 is the penalty when the requirements are not met.
end
end

function Eval_all=XT_all(power_beam_PMN,split_PMN,chrom,tx_num,min_channel_power)
[nind,~]=size(chrom);
Eval_all=ones(nind,1);
for i=1:nind
    chrom_ind=chrom(i,:);
    tx_No=chrom_ind(:,1:tx_num);
    SMG=chrom_ind(:,tx_num+1);
    powerPMN_tx=power_beam_PMN(:,tx_No);
    split_PMN_ind=split_PMN(SMG,1:tx_num);
    Eval_all(i)=XT(powerPMN_tx,split_PMN_ind,min_channel_power);
end
end