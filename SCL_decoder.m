function [polar_info_esti,u_all,active_length] = SCL_decoder(llr, L, K, ...
    frozen_bits,  lambda_offset, llr_layer_vec, bit_layer_vec,c,M)
%LLR-based SCL deocoder, 一个函数，没有调用子函数


%常量
N = length(llr);
m = log2(N);

%memory declared
%Lazy Copy 记录该译码器所需数据来自哪个译码器 
% Here,data refer to LLRs and partial sums.
%If you do not understand such operation, you can directly copy data.
lazy_copy = zeros(m, L);
P = zeros(N - 1, L); %不包含信道llr初始值因此N-1行足够
 C = zeros(2*N-1,2*L);%不估计 (x1, x2, ... , xN), so N - 1 is enough.
%C = ones(N - 1, 2 * L)*6;
u_all = zeros(K-M, L);%记录信息比特译码结果的数组，可有可无
activepath = zeros(L, 1);%标注每个SC译码器是否被激活. '1'被激活; '0' otherwise.
cnt_u = 1;%信息位计数
u_p=1;%校验位计数
active_length=zeros(1,M);
p_2=zeros(M,L);
% error=0;
%initialize，SC1译码器被激活
PM=zeros(L,1);%路径度量值
activepath(1) = 1;
lazy_copy(:, 1) = 1;%lazy_copy中的值代表新生路径其所需要的数据从哪里读取
Lp=4;

%%
%decoding starts
%default: in the case of path clone, the origianl path...
%always corresponds to bit 0,while the new path bit 1.
for phi = 0 : N - 1
    layer = llr_layer_vec(phi + 1);
    phi_mod_2 = mod(phi, 2);%可判断u下标的奇偶性
    
    for l_index = 1 : L%SC译码器序号
        if activepath(l_index) == 0%若该序号的SC译码器未被激活，则结束此次循环
            continue;
        end
        
        %Decoding bits u_0 and u_N/2 needs channel LLR,...
        %so the decoding of them is separated from other bits. 
        switch phi
            case 0
                index_1 = lambda_offset(m);
                %用llr的f运算
                for beta = 0 : index_1 - 1
                    P(beta + index_1, l_index) = sign(llr(beta + 1)) *...
                        sign(llr(beta + index_1 + 1)) *...
                        min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                end
                for i_layer = m - 2 : -1 : 0%判断剩下f运算的次数
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
                        %用P向量的f运算
                        P(beta + index_1, l_index) = sign(P(beta + index_2, l_index)) *...
                            sign(P(beta + index_1 + index_2, l_index)) *...
                            min(abs(P(beta + index_2, l_index)), abs(P(beta + index_1 + index_2, l_index)));
                    end
                end
                
            case N/2
                index_1 = lambda_offset(m);
                 %用llr的g运算
                for beta = 0 : index_1 - 1
                    x_tmp = C(beta + index_1, 2 * l_index - 1);
                    P(beta + index_1, l_index) = (1 - 2 * x_tmp) *...
                        llr(beta + 1) + llr(beta + 1 + index_1);
                end
                for i_layer = m - 2 : -1 : 0
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    %用P的f运算
                    for beta = 0 : index_1 - 1
                        P(beta + index_1, l_index) = ...
                            sign(P(beta + index_2, l_index)) *...
                            sign(P(beta + index_1 + index_2, l_index)) * ...
                            min(abs(P(beta + index_2, l_index)), abs(P(beta + index_1 + index_2, l_index)));
                    end
                end
                
            otherwise
                index_1 = lambda_offset(layer + 1);
                index_2 = lambda_offset(layer + 2);
                %用P向量的g运算
                for beta = 0 : index_1 - 1
                    P(beta + index_1, l_index) = (1 - 2 *...
                        C(beta + index_1, 2 * l_index - 1)) *...
                        P(beta + index_2, lazy_copy(layer + 2, l_index)) +...
                        P(beta + index_1 + index_2, lazy_copy(layer + 2, l_index));
                end
                % 用P向量的f运算                
                for i_layer = layer - 1 : -1 : 0
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
                        P(beta + index_1, l_index) = sign(P(beta + index_2, l_index)) *...
                            sign(P(beta + index_1 + index_2, l_index)) * min(abs(P(beta + index_2, l_index)),...
                            abs(P(beta + index_1 + index_2, l_index)));
                    end
                end
        end
    end
    
    %%
    %当该u是不可靠信息位时
    if frozen_bits(phi + 1) == 3%if now we decode an info bit
        %信息比特需要进行克隆，激活一个SC译码器
        PM_pair = realmax * ones(2, L);%realmax返回指定浮点数类型所能表示的正的最大值。
        
        for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            %计算复制的PM值，判断有无加惩罚值
            if P(1, l_index) >= 0%LNk>=0时
                %uk=0,满足硬判决，原PM值
                PM_pair(1, l_index) = PM(l_index);
                %uk=1,不满足硬判决，+|LNk|
                PM_pair(2, l_index) = PM(l_index) + P(1, l_index);
            else%LNk<0时
                %uk=0,不满足硬判决，+|LNk|
                PM_pair(1, l_index) = PM(l_index) - P(1, l_index);
                %uk=1,满足硬判决，原PM值
                PM_pair(2, l_index) = PM(l_index);
            end
        end
        
        %%
        %不能超过L个PM_pair值，若超过就要删选取较小的L个PM_pair值
%         activepath
        middle = min(2 * sum(activepath), L);
        PM_sort = sort(PM_pair(:));%将矩阵PM_pair转换成列向量，再由小到大排列
        PM_cv = PM_sort(middle);%可保留路径的最大值
        compare = PM_pair <= PM_cv; 
        
        %记录被删除的路径
        kill_index = zeros(L, 1);%to record the index of the path that is killed
        %被删除路径的总个数
        kill_cnt = 0;%the total number of killed path
        %上述两个变量由一个堆栈组成
        
        for i = 1 : L
            %判断是否有对应2条路径都要被删除的SC译码器
            if (compare(1, i) == 0)&&(compare(2, i) == 0)
                activepath(i) = 0;%对应SC译码器死亡
                kill_cnt = kill_cnt + 1;%堆栈，死亡译码器的个数
                kill_index(kill_cnt) = i;%将死亡的SC译码器序号i放入kill_index中
            end
        end
        
        %%
        %判断保留下来的SC译码器的两条路径是否分别保留，若保留对应的compare值为1
        for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            %一个*2是为了区分两条路径,说明2是保留第一条路径，1是保留第二条路径
            path_state = compare(1, l_index) * 2 + compare(2, l_index);
            
            %不考虑path_state=0，因为在上面已经考虑过了
            switch path_state
                case 1%只保留第2条路径
                    u_all(cnt_u, l_index) = 1;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 1;
                    PM(l_index) = PM_pair(2, l_index);
                case 2%只保留第1条路径
                    u_all(cnt_u, l_index) = 0;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 0;
                    PM(l_index) = PM_pair(1, l_index);
                case 3%保留2条路径，需弹出一个译码器存储uk=1
                    index = kill_index(kill_cnt);%将堆栈里最后一个译码器弹出来，index为该序号
                    kill_cnt = kill_cnt - 1;%剩余死亡SC译码器个数
                    activepath(index) = 1;%将弹出来的译码器激活
                    
                    %lazy copy更新
                    %将原译码器LC的一列复制给激活译码器对应的列
                    lazy_copy(:, index) = lazy_copy(:, l_index);
                    u_all(:, index) = u_all(:, l_index);%u复制到被激活的译码器对应u中
                    u_all(cnt_u, l_index) = 0;%原译码器，令u=0
                    u_all(cnt_u, index) = 1;%被激活译码器，令u=1
                    C(1, 2 * l_index - 1 + phi_mod_2) = 0;%对应C，u=0
                    C(1, 2 * index - 1 + phi_mod_2) = 1;%对应C，u=1
                    PM(l_index) = PM_pair(1, l_index);%将PM_pair值存入对应PM中
                    PM(index) = PM_pair(2, l_index);
            end
        end
        cnt_u = cnt_u + 1;
     
       %%
      % 当是冻结位时
    else if frozen_bits(phi + 1) == 1||frozen_bits(phi + 1) == 11%frozen bit operation
        for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            if frozen_bits(phi + 1) == 11
            if P(1, l_index) < 0%LNk<0与uk=0,不满足硬判决需加惩罚
                PM(l_index) = PM(l_index) - P(1, l_index);%PM+|LNk|
            end
            end
            if phi_mod_2 == 0%判断u下标是否为奇数
                C(1, 2 * l_index - 1) = 0;%是奇数，将uk=0放在C[1]的左边一位
            else
                C(1, 2 * l_index) = 0;%是偶数，将uk=0放在C[1]的右边一位
            end 
        end
     else if frozen_bits(phi + 1) == 2%校验比特 
             active_index=[];%记录能通过奇偶校验比特的路径数
             for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            p=0;
               for k=1:length(c{u_p})              
                   p=p+u_all(c{u_p}(k),l_index);
               end
               p_2(u_p,l_index)=mod(p,2);
                if P(1,l_index)>0
                   p_3=0;
               else
                   p_3=1;
                end
%               p_1(u_p)
%                if(p_1(u_p))~=p_3
%                    error=error+1
%                end
%                if p_2(u_p,l_index)~=p_1(u_p)%将不通过奇偶校验的路径删除，仅保留通过奇偶校验的路径
               if p_2(u_p,l_index)~=p_3%不通过奇偶校验的路径加惩罚值
                    PM(l_index) = PM(l_index) +abs(P(1, l_index));
%                      activepath(l_index) = 0;%对应SC译码器死亡
                  continue
               end
               active_index=[active_index,l_index];%记录通过奇偶校验的路径序号 
             end
             active_length(u_p)=length(active_index);
             
           %保留Lp条通过奇偶校验的路径     
           if length(active_index)>=Lp %有多条路径通过奇偶校验
               for l_index = 1:L
                   if activepath(l_index) == 0
                continue;
                   end
               end
                     [~,index]=sort(PM(active_index));
                    PM_min_index = active_index(index(1:Lp));
%                   PM_min=PM(active_index(1));
%                   PM_min_index=active_index(1);
%                   for j=2:length(active_index)
%                       if PM(active_index(j))<PM_min;
%                           PM_min=PM(active_index(j));
%                           PM_min_index=active_index(j);
%                       end
%                   end
                  for l_index=1:L
                  if ismember(l_index,PM_min_index)==0   %保留路径度量最小的Lp条通过路径
                      activepath(l_index) = 0;%对应SC译码器死亡
%                       kill_cnt = kill_cnt + 1;%堆栈，死亡译码器的个数
%                       kill_index(kill_cnt) = l_index;%将死亡的SC译码器序号l_index放入kill_index中
                      continue
%                         lazy_copy(:, l_index) = lazy_copy(:, PM_min_index);
%                         u_all(:,l_index)=u_all(:,PM_min_index);
                  end
                  if phi_mod_2 == 0%判断u下标是否为奇数
                     C(1, 2 * l_index - 1) = p_2(u_p,l_index);%是奇数，将uk=0放在C[1]的左边一位
                    else
                    C(1, 2 * l_index) = p_2(u_p,l_index);%是偶数，将uk=0放在C[1]的右边一位
                  end  
                     
                 
                  end
                  
              
           elseif length(active_index)==0%奇偶校验均未通过,保留PM最小的Lp条
%                index=3
%                activepath
           
                 row=find(activepath==1);%在存活的路径中将PM排序
%                  PM
%                   PM(row)
                 
     
               [~, min_index_11]=sort(PM(row));%min_index_1为最佳路径的序号
               min_index_1=row(min_index_11(1:Lp));
%                min_index_1
               for l_index = 1 : L
                    if activepath(l_index) == 0
                continue;
                    end
                     if ismember(l_index,min_index_1)==0
                     activepath(l_index) = 0;%对应SC译码器死亡
%                       kill_cnt = kill_cnt + 1;%堆栈，死亡译码器的个数
%                       kill_index(kill_cnt) = l_index;%将死亡的SC译码器序号l_index放入kill_index中
                      continue
                     end
%                      lazy_copy(:, l_index) = lazy_copy(:,min_index_1);
%                         u_all(:,l_index)=u_all(:,min_index_1);
                if phi_mod_2 == 0%判断u下标是否为奇数
                     C(1, 2 * l_index - 1) = p_2(u_p,l_index);%是奇数，将uk=0放在C[1]的左边一位
                    else
                    C(1, 2 * l_index) = p_2(u_p,l_index);%是偶数，将uk=0放在C[1]的右边一位
                end 
                    
               end
                
              
                 
              
           else%  1≤通过奇偶校验路径数＜Lp
%                for l_index = 1 : L
%                     if activepath(l_index) == 0
%                 continue;
%                     end
                    len = length(active_index);
                   min_index=find(activepath==1);%在存活的路径中将PM排序
                    [ia, ~] = setdiff(min_index, active_index);%在存活路径中将通过奇偶校验的路径去除
%                     res = min_index(sort(ia));
                     [~, min_index_]=sort(PM(ia));%在剩余路径中找PM最小的Lp-len条路径保留
               min_index_12=ia(min_index_(1:Lp-len));
               min_index_12=min_index_12';
               real_index = [active_index min_index_12];
                for l_index=1:L
                  if ismember(l_index, real_index)==0   %保留路径度量最小的Lp条通过路径
                      activepath(l_index) = 0;%对应SC译码器死亡
%                       kill_cnt = kill_cnt + 1;%堆栈，死亡译码器的个数
%                       kill_index(kill_cnt) = l_index;%将死亡的SC译码器序号l_index放入kill_index中
                      continue
%                         lazy_copy(:, l_index) = lazy_copy(:, PM_min_index);
%                         u_all(:,l_index)=u_all(:,PM_min_index);
                  end
%                    lazy_copy(:, l_index) = lazy_copy(:, active_index);
%                         u_all(:,l_index)=u_all(:,active_index);
                  if phi_mod_2 == 0%判断u下标是否为奇数
                     C(1, 2 * l_index - 1) = p_2(u_p,l_index);%是奇数，将uk=0放在C[1]的左边一位
                    else
                C(1, 2 * l_index) = p_2(u_p,l_index);%是偶数，将uk=0放在C[1]的右边一位
                  end 
           
               end
           
                
                  
                
                  
           end
        
              
              
              u_p=u_p+1;
          else %较可靠信息位不进行路径扩张，SC译码
             for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            u_phi=P(1,l_index)<0;
            u_all(cnt_u,l_index)=u_phi;
              C(1, 2 * l_index - 1 + phi_mod_2) = u_phi;
             end
               cnt_u=cnt_u+1;
         end
         
        end 
    end
   
    %%
    %partial-sum return返回值x
    for l_index = 1 : L
        if activepath(l_index) == 0
            continue
        end
        %判断u下标是否是偶数，为偶数时需要更新C向量
        if (phi_mod_2  == 1) && (phi ~= N - 1)
            layer = bit_layer_vec(phi + 1);
            %对应SC译码器的C向量第二列赋值
            for i_layer = 0 : layer - 1
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : 2 * index_1 - 1
                    C(beta + index_1, 2 * l_index) = mod(C(beta, 2 *  lazy_copy(i_layer + 1, l_index) - 1) + C(beta, 2 * l_index), 2);%Left Column lazy copy
                    C(beta + index_2, 2 * l_index) = C(beta, 2 * l_index);   
                end
            end
            %对应SC译码器的C向量第一列赋值，C更新结束
            index_1 = lambda_offset(layer + 1);
            index_2 = lambda_offset(layer + 2);
            for beta = index_1 : 2 * index_1 - 1
                C(beta + index_1, 2 * l_index - 1) = mod(...
                    C(beta, 2 * lazy_copy(layer + 1, l_index) - 1) + ...
                    C(beta, 2 * l_index), 2);
                %Left Column lazy copy
                C(beta + index_2, 2 * l_index - 1) = C(beta, 2 * l_index);
            end 
        end
    end
   
    %%
    %lazy copy更新
    %C数组计算完毕，重新判断激活路径对应的LC中的列
    if phi < N - 1
        %下一个比特二进制展开后低位有m个连续0，m+1
        for i_layer = 1 : llr_layer_vec(phi + 2) + 1
            %对应前m+1行元素修改为对应激活译码器的序号
            for l_index = 1 : L
                lazy_copy(i_layer, l_index) = l_index;
            end
        end
    end
end

%%
%path selection.
%选择最佳路径，即PM值最小的路径
% [~, min_index]=min(PM);%min_index为最佳路径的序号
% polar_info_esti = u_all(:,min_index);%最佳路径的信息位u取出（只有信息位的译码结果,与信源中的信息位比较即可）
final_index=find(activepath==1);
if length(final_index) ~=1
   [~, final_index1]=min(PM(final_index));
   final_index = final_index(final_index1);
end
% final_index=final_index+1;
polar_info_esti = u_all(:,final_index);
end
