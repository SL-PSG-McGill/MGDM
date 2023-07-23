%{ 
 Author: Shanglin Li, Photonic DataCom Team, McGill University 
 Date: 04/22/2022
 Note:
 Please install The Genetic Algorithm Toolbox for MATLAB developed by
 The University of Sheffield (http://uos-codem.github.io/GA-Toolbox/),
 before running this script.
%}
clear variables
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
    error('Exit the program!')
end
load(Dataset_name)
power_beam_PMN=power_beam_PMN';
[PMN_num,beam_info_num]=size(power_beam_PMN);
tx_num = input("How many fiber channels are excited? Please input 3 or 4.Input 0 to exit the program.\n");
while ~ismember(tx_num, [3,4,0])
    tx_num =input("Your input is WRONG! Please input 3 or 4. Input 0 to exit the program. \n");
end
split_PMN=[];
loop1=1;
if tx_num==3
    for sec1=1:1:PMN_num-tx_num+1
        for  sec2=1:1:PMN_num-sec1-tx_num+2
            sec3=PMN_num-sec1-sec2;
            split_PMN(loop1,1:tx_num)=[sec1,sec2,sec3];
            loop1=loop1+1;
        end
    end
elseif tx_num==4
    for sec1=1:1:PMN_num-tx_num+1
        for  sec2=1:1:PMN_num-sec1-tx_num+2
            for sec3=1:1:PMN_num-sec1-sec2-tx_num+3
                sec4=PMN_num-sec1-sec2-sec3;
                split_PMN(loop1,1:tx_num)=[sec1,sec2,sec3,sec4];
                loop1=loop1+1;
            end
        end
    end
elseif tx_num==0
    error('Exit the program!')
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
MIGGEN=1; % Migration occurs for every MIGGEN generations
MUTR=0.2; % Mutation rate
BaseV=crtbase([tx_num,1],[beam_info_num,split_num]);
[Chrom,Lind]=crtbp(NIND,BaseV);% Initial population
Chrom_addone=Chrom+ones(NIND,tx_num+1);
gen_count=0;
ObjV=XT_all(power_beam_PMN,split_PMN,Chrom_addone,tx_num,min_channel_power,max_channel_XT);
while gen_count<MAXGEN
    FitnV=ranking(ObjV,[2,1],SUBPOP); % Calculate fitness scores
    Selch1=select('sus',Chrom,FitnV,GGAP,SUBPOP); % Selection
    Selch2=recombin('xovsh',Selch1,XOVR,SUBPOP); % Recombination
    Selch3=mutate('mut',Selch2,BaseV,MUTR,SUBPOP);% Mutation
    ObjV_Sel=XT_all(power_beam_PMN,split_PMN,Selch3+1,tx_num,min_channel_power,max_channel_XT); % Calculate the overall inter-channel crosstalk
    [Chrom, ObjV]=reins(Chrom,Selch3,SUBPOP,[1,INSR],ObjV,ObjV_Sel); % Re-insertion
    gen_count=gen_count+1;
    if (rem(gen_count,MIGGEN)==0)
        [Chrom, ObjV]=migrate(Chrom,SUBPOP,[MIGR,1,1],ObjV); % Migration
    end
    waitbar(gen_count/MAXGEN,f,['Running ',num2str(gen_count/MAXGEN*100,'%.1f'),'%'])
end
waitbar(1,f,'The calculation of launch conditions is done!')
[Best_XT_result,I]=min(ObjV);
Best_Ind_GA=Chrom(I,:)+1;
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
function Eval=XT(powerPMN_tx,split_PMN_ind,min_channel_power,max_channel_XT) % Target Function
powerPMN_tx_cell=mat2cell(powerPMN_tx,split_PMN_ind);
block_sum=cell2mat(cellfun(@(x)sum(x,1),powerPMN_tx_cell,'UniformOutput',false));
block_sum_XT=tril(block_sum,-1)+triu(block_sum,1);
target_power=diag(block_sum);
Rel_XT=(sum(block_sum,2)-target_power)./target_power; % Relative crostalk
if all(target_power>=min_channel_power) && all(block_sum_XT(:)<=max_channel_XT)  % the requirements of the power coupling coefficient matrices
    Eval=10*log10(sum(Rel_XT)+1e-60);
else
    Eval= 10*log10(sum(Rel_XT)+1e-60)+1e6; % 1e6 is the penalty when the requirements are not met.
end
end

function Eval_all=XT_all(power_beam_PMN,split_PMN,chrom,tx_num,min_channel_power,max_channel_XT)
[nind,~]=size(chrom);
Eval_all=ones(nind,1);
for i=1:nind
    chrom_ind=chrom(i,:);
    tx_No=chrom_ind(:,1:tx_num);
    SMG=chrom_ind(:,tx_num+1);
    powerPMN_tx=power_beam_PMN(:,tx_No);
    split_PMN_ind=split_PMN(SMG,1:tx_num);
    Eval_all(i)=XT(powerPMN_tx,split_PMN_ind,min_channel_power,max_channel_XT);
end
end