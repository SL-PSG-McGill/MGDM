% runtime=datestr(now)
addpath(genpath(pwd)) 
% clear variables
if ~exist('Optimal_Chrom','var')
    Optimal_Chrom=[]; % initialize Optimal Chrome
end
% Optimal_Chrom=[37240,23271,67021,78250,13938,89367,913]; % 强行引入最佳的种群，帮助收敛

newFilename = '0dataset_1310nm_w0_2-0_5-20_aoff_0-0_5-30_psi_-15-0_5-15_theta_0-5-0_radius_31_25';
load(newFilename)

% 直接到min_channel_power_final有时无法找到最优点，所以以0.1为步进值，从channel_power_start始一步步迭代
channel_power_start=0.3; % tx的功率迭代的起始点
channel_power_final=0.5; %要求目标tx功率要大于min_channel_power_final
max_channel_XT=0.2;

power_beam_PMN=power_beam_PMN';
[PMN_num,beam_info_num]=size(power_beam_PMN);
[split_num,tx_num]=size(split_PMN);
max_bin_beam_str=dec2bin(beam_info_num);
max_bin_beam=max_bin_beam_str-'0';
max_bin_beam_length=length(max_bin_beam);
max_bin_MG_str=dec2bin(split_num);
max_bin_MG=max_bin_MG_str-'0';
max_bin_MG_length=length(max_bin_MG);
%% Multi-group Genetic Algorithm  MATLAB遗传算法工具箱及应用_雷英杰 P117
tic
% max_dec_beam=bin2dec(max_bin_beam_str);
% max_dec_MG=bin2dec(max_bin_MG_str);
% NIND=floor(beam_info_num*split_num/1e4)*10; % the number of individual; comment:数字越大越好,可以避免局部收敛, 但是运算速度会变慢 (200 for tx_num=5)
NIND=40000;
MAXGEN=200; % max numbers of generations； comment: 收敛到一定地步, 增大MAXGEN已无意义
NVAR=tx_num+1; % 变量数目，tx_num+1 (1为PMN的分割方法)
GGAP=0.8; %代沟(Generation gap)
XOVR=1; %交叉率
INSR=0.9; %插入率
SUBPOP=NIND/400; %子种群数
MIGR=0.4; %迁移率
MIGGEN=2; %每MIGGEN代迁移一次
MUTR=0.2; %变异率 还有用1/NVAR;默认0.7/NVAR, comment:感觉越大越不准
BaseV=crtbase([tx_num,1],[beam_info_num,split_num]);

channel_power=channel_power_start;
%----------------------------------------
min_channel_power=channel_power_final;

% for loop=1:10 %循环多次，避免受初始值Chrom影响过大，陷入局部最优
trace=zeros(MAXGEN,2);%遗传算法性能跟踪初始值
% FieldDR=[repmat([1;beam_info_num],1,tx_num),[1;split_num]];
% Chrom=round(crtrp(NIND,FieldDR));

[Chrom,Lind]=crtbp(NIND,BaseV);%初始种群
if ~isempty(Optimal_Chrom)
    Chrom(1,:)=Optimal_Chrom; %把头n行替换为最优化的值
end
Chrom_addone=Chrom+ones(NIND,tx_num+1); %对应matlab向量坐标从1开始

gen=0; %代数
ObjV=XT_all(power_beam_PMN,split_PMN,Chrom_addone,tx_num,min_channel_power,max_channel_XT);
[Best_ObjV,I]=min(ObjV);
Optimal_Chrom=Chrom(I,:);
while gen<MAXGEN
    FitnV=ranking(ObjV,[2,1],SUBPOP); %分配适应度值
    Selch1=select('sus',Chrom,FitnV,GGAP,SUBPOP); %选择
    Selch2=recombin('xovsh',Selch1,XOVR,SUBPOP); %重组
    Selch3=mutate('mut',Selch2,BaseV,MUTR,SUBPOP);%变异
    Selch3_addone= Selch3+1;
    ObjV_Sel=XT_all(power_beam_PMN,split_PMN,Selch3_addone,tx_num,min_channel_power,max_channel_XT); %计算子代目标函数值
    [Chrom, ObjV]=reins(Chrom,Selch3,SUBPOP,[1,INSR],ObjV,ObjV_Sel); %重插入
    if (rem(gen,MIGGEN)==0)
        [Chrom, ObjV]=migrate(Chrom,SUBPOP,[MIGR,1,1],ObjV);
    end
    [Best_ObjV_new,I]=min(ObjV);
    if Best_ObjV_new<=Best_ObjV
        Best_ObjV=Best_ObjV_new;
        Optimal_Chrom=Chrom(I,:);
    else
        Chrom(1,:)=Optimal_Chrom;
        ObjV(1,:)=Best_ObjV;
    end       
    gen=gen+1;
    trace(gen,1)=Best_ObjV;  %遗传算法性能跟踪
    trace(gen,2)=sum(ObjV)/length(ObjV);
end
%      end
Best_Ind_GA=Optimal_Chrom+1; %note:已经加1了
Best_beam_info_GA=beam_info(Best_Ind_GA(1:tx_num)',:);
Best_PMN_Dist_GA=split_PMN(Best_Ind_GA(end),:)
Best_powerPMN_tx_GA=power_beam_PMN(:,Best_Ind_GA(1:tx_num));
Best_powerPMN_tx_cell_GA=mat2cell(Best_powerPMN_tx_GA,Best_PMN_Dist_GA); % 按照super-mode-group 概念分成块
Best_block_sum_GA=cell2mat(cellfun(@(x)sum(x,1),Best_powerPMN_tx_cell_GA,'UniformOutput',false)) % 按照super-mode-group 概念求和,(super-mode-group*tx);最佳信道功率分布
Best_XT_GA=(sum(Best_block_sum_GA,2)-diag(Best_block_sum_GA))./diag(Best_block_sum_GA);
Best_block_sum_db_GA=10*log10(Best_block_sum_GA)
% beam info translation
Best_beam_info_result=zeros(size(Best_beam_info_GA));
Best_beam_info_result(:,1:2)=Best_beam_info_GA(:,1:2)*1e6;
Best_beam_info_result(:,3:4)=Best_beam_info_GA(:,3:4)*180/pi;
Best_beam_info_result
toc

% figure
% plot(trace(:,1),'-.');hold on
% plot(trace(:,2));
% xlabel('Generation');
% ylabel('Crosstalk (dB)');
% legend('Minimum','Mean');

%% Simulated Annealing Algorithm
% tic
% chain_length=200; %马尔科夫链长度
% Step=0.05; % 步长因子
% T0=1000; %初始温度
% Tend=1e-3; %终止温度
% At=0.9; %衰减参数
% YZ=1e-8; % 容差
% 
% Pre_beam_info=Best_beam_info_GA;
% Pre_PMN=Best_PMN_Dist_GA;
% Pre_powerPMN=Best_powerPMN_tx_GA;
% Best_beam_info=Pre_beam_info;
% Best_PMN=Pre_PMN;
% Best_powerPMN=Pre_powerPMN;
% Pre_Best_beam_info=Best_beam_info;
% Pre_Best_powerPMN=Best_powerPMN;
% Pre_Best_PMN=Best_PMN;
% 
% delta=1;
% while  (delta>YZ && T0>Tend)
%     for Markov_loop=1:chain_length
%         Next_beam_info=zeros(size(Pre_beam_info));
%         Next_powerPMN=zeros(size(Pre_powerPMN));
%         %  Next_PMN=Pre_PMN;
%         [Next_PMN,Next_PMN_Split]=rand_PMN(Pre_PMN,PMN_num,tx_num);
%         
%         for loop1=1:tx_num
%             p=0;
%             while p==0
%                 Pre_w0=Pre_beam_info(loop1,1);
%                 Pre_a_off=Pre_beam_info(loop1,2);
%                 Pre_psi=Pre_beam_info(loop1,3);
%                 Pre_theta=Pre_beam_info(loop1,4);
%                 Next_w0=Pre_w0+Step*(rand-0.5)*(w0_end-w0_start)*1e-6;
%                 Next_a_off=Pre_a_off+Step*(rand-0.5)*(a_off_end-a_off_start)*1e-6;
%                 Next_psi=Pre_psi+Step*(rand-0.5)*(psi_end-psi_start)*pi/180;
%                 Next_theta=Pre_theta+Step*(rand-0.5)*(theta_end-theta_start)*pi/180;
%                 if (Next_w0>=w0_start*1e-6 && Next_w0<=w0_end*1e-6 && Next_a_off>=a_off_start*1e-6&& Next_a_off<=a_off_end*1e-6...
%                         && Next_psi>=psi_start*pi/180 && Next_psi<=psi_end*pi/180 && Next_theta>=theta_start*pi/180 && Next_theta<=theta_end*pi/180)
%                     p=1;
%                     Next_beam_info(loop1,1)=Next_w0;
%                     Next_beam_info(loop1,2)=Next_a_off;
%                     Next_beam_info(loop1,3)=Next_psi;
%                     Next_beam_info(loop1,4)=Next_theta;
%                 end
%             end
%             Next_powerPMN(:,loop1)=MPD_func(para,Next_w0,Next_a_off,Next_psi,Next_theta);
%         end
%         
%         if (XT(Next_powerPMN,Next_PMN,min_channel_power,max_channel_XT)...
%                 <XT(Best_powerPMN,Best_PMN,min_channel_power,max_channel_XT) )
%             Pre_Best_beam_info=Best_beam_info;
%             Pre_Best_powerPMN=Best_powerPMN;
%             Pre_Best_PMN=Best_PMN;
%             Best_beam_info=Next_beam_info;
%             Best_powerPMN=Next_powerPMN;
%             Best_PMN=Next_PMN;
%         end
%         % Metropolis Process
%         if (XT(Next_powerPMN,Next_PMN,min_channel_power,max_channel_XT)...
%                 <XT(Pre_powerPMN,Pre_PMN,min_channel_power,max_channel_XT) )
%             Pre_beam_info= Next_beam_info;
%             Pre_powerPMN=Next_powerPMN;
%             Pre_PMN=Next_PMN;
%         else
%             changer=-1*(XT(Next_powerPMN,Next_PMN,min_channel_power,max_channel_XT)...
%                 -XT(Pre_powerPMN,Pre_PMN,min_channel_power,max_channel_XT))/T0;
%             p1=exp(changer);
%             %%%%%%接受较差的解%%%%%%%%%%
%             if p1>rand
%                 Pre_beam_info= Next_beam_info;
%                 Pre_powerPMN=Next_powerPMN;
%                 Pre_PMN=Next_PMN;
%             end
%         end
%     end
%     
%     delta=abs(XT(Best_powerPMN,Best_PMN,min_channel_power,max_channel_XT) ...
%         -XT(Pre_Best_powerPMN,Pre_Best_PMN,min_channel_power,max_channel_XT));
%     T0=T0*At;
% end
% toc
% Best_powerPMN_tx_cell_SA=mat2cell(Best_powerPMN,Best_PMN); % 按照super-mode-group 概念分成块
% Best_block_sum_SA=cell2mat(cellfun(@(x)sum(x,1),Best_powerPMN_tx_cell_SA,'UniformOutput',false)) % 按照super-mode-group 概念求和,(super-mode-group*tx);最佳信道功率分布
% Best_XT_SA=(sum(Best_block_sum_SA,2)-diag(Best_block_sum_SA))./diag(Best_block_sum_SA);
% Best_block_sum_db_SA=10*log10(Best_block_sum_SA);
%% sub-function
function Eval=XT(powerPMN_tx,split_PMN_ind,min_channel_power,max_channel_XT) %目标函数
powerPMN_tx_cell=mat2cell(powerPMN_tx,split_PMN_ind); % 按照super-mode-group 概念分成块
block_sum=cell2mat(cellfun(@(x)sum(x,1),powerPMN_tx_cell,'UniformOutput',false));% 按照super-mode-group 概念求和,(super-mode-group*tx)
block_sum_XT=tril(block_sum,-1)+triu(block_sum,1);
target_power=diag(block_sum);
% block_sum_norm=block_sum./target_power;，，，，，，，
% block_sum_norm(logical(eye(size(block_sum_norm))))=-1;%对角线元素变为-1
XT=(sum(block_sum,2)-target_power)./target_power; % default, relative XT
% XT=(sum(block_sum,2)-target_power); % absolute XT
%%
% if all(target_power>=min_channel_power) && all(block_sum_XT(:)<=max_channel_XT)  %要求目标tx功率要大于0.1 即大于-10 dB
 if all(target_power>=min_channel_power)  %要求目标tx功率要大于0.1 即大于-10 dB
%     Eval=10*log10(sum(XT)+1e-60); % XT; the smaller, the better
    Eval=sum(-10*log10(target_power+1e-60));
else
%     Eval= 10*log10(sum(XT)+1e-60)+1e6; % 1e6为不满足条件的惩罚
    Eval=sum(-10*log10(target_power+1e-60)+1e10);
end
%%
end

function Eval_all=XT_all(power_beam_PMN,split_PMN,chrom,tx_num,min_channel_power,max_channel_XT)
[nind,~]=size(chrom);
Eval_all=ones(nind,1);
for i=1:nind
    chrom_ind=chrom(i,:);
    tx_No=chrom_ind(:,1:tx_num);
    SMG=chrom_ind(:,tx_num+1);   % Super Mode Group
    powerPMN_tx=power_beam_PMN(:,tx_No);
    split_PMN_ind=split_PMN(SMG,1:tx_num);
    Eval_all(i)=XT(powerPMN_tx,split_PMN_ind,min_channel_power,max_channel_XT);
end
end

function  y=istxsame(bin_variables,max_bin_beam_length,tx_num)  %判断前5个max_bin_beam_length不同
tx_matrix=reshape(bin_variables(:,1:tx_num*max_bin_beam_length),max_bin_beam_length,tx_num);
C=unique(tx_matrix','rows');
if size(C,1)==tx_num
    y=true(1);
else
    y=false(1);
end
end

function [New_PMN_Dist,New_PMN_Split]=rand_PMN(PMN_Dist,PMN_num,tx_num)
PMN_Split=zeros(1,tx_num);
PMN_Split(1)=PMN_Dist(1);
for ite1=2:tx_num
    PMN_Split(ite1)=PMN_Split(ite1-1)+PMN_Dist(ite1);
end

New_PMN_Split=zeros(1,tx_num);
New_PMN_Dist=zeros(1,tx_num);
PMN_Split_value=PMN_Dist(1)+round(rand*2-1);
if PMN_Split_value>0 && PMN_Split_value<=PMN_num
    New_PMN_Split(1)=PMN_Split_value;
    New_PMN_Dist(1)=PMN_Split_value;
else
    New_PMN_Split(1)=PMN_Split(1);
    New_PMN_Dist(1)=PMN_Dist(1);
end
for ite1=2:tx_num-1
    PMN_Split_value=PMN_Split(ite1)+round(rand*2-1);
    if PMN_Split_value>New_PMN_Split(ite1-1)&& PMN_Split_value<PMN_num
        New_PMN_Split(ite1)=PMN_Split_value;
        New_PMN_Dist(ite1)=PMN_Split_value-New_PMN_Split(ite1-1);
    else
        New_PMN_Split(ite1)=PMN_Split(ite1);
        New_PMN_Dist(ite1)=PMN_Split(ite1)-New_PMN_Split(ite1-1);
    end
end
New_PMN_Split(tx_num)=PMN_num;
New_PMN_Dist(tx_num)=PMN_num-New_PMN_Split(tx_num-1);
end