function [polar_info_esti,u_all,active_length] = SCL_decoder(llr, L, K, ...
    frozen_bits,  lambda_offset, llr_layer_vec, bit_layer_vec,c,M)
%LLR-based SCL deocoder, һ��������û�е����Ӻ���


%����
N = length(llr);
m = log2(N);

%memory declared
%Lazy Copy ��¼���������������������ĸ������� 
% Here,data refer to LLRs and partial sums.
%If you do not understand such operation, you can directly copy data.
lazy_copy = zeros(m, L);
P = zeros(N - 1, L); %�������ŵ�llr��ʼֵ���N-1���㹻
 C = zeros(2*N-1,2*L);%������ (x1, x2, ... , xN), so N - 1 is enough.
%C = ones(N - 1, 2 * L)*6;
u_all = zeros(K-M, L);%��¼��Ϣ���������������飬���п���
activepath = zeros(L, 1);%��עÿ��SC�������Ƿ񱻼���. '1'������; '0' otherwise.
cnt_u = 1;%��Ϣλ����
u_p=1;%У��λ����
active_length=zeros(1,M);
p_2=zeros(M,L);
% error=0;
%initialize��SC1������������
PM=zeros(L,1);%·������ֵ
activepath(1) = 1;
lazy_copy(:, 1) = 1;%lazy_copy�е�ֵ��������·��������Ҫ�����ݴ������ȡ
Lp=4;

%%
%decoding starts
%default: in the case of path clone, the origianl path...
%always corresponds to bit 0,while the new path bit 1.
for phi = 0 : N - 1
    layer = llr_layer_vec(phi + 1);
    phi_mod_2 = mod(phi, 2);%���ж�u�±����ż��
    
    for l_index = 1 : L%SC���������
        if activepath(l_index) == 0%������ŵ�SC������δ�����������˴�ѭ��
            continue;
        end
        
        %Decoding bits u_0 and u_N/2 needs channel LLR,...
        %so the decoding of them is separated from other bits. 
        switch phi
            case 0
                index_1 = lambda_offset(m);
                %��llr��f����
                for beta = 0 : index_1 - 1
                    P(beta + index_1, l_index) = sign(llr(beta + 1)) *...
                        sign(llr(beta + index_1 + 1)) *...
                        min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                end
                for i_layer = m - 2 : -1 : 0%�ж�ʣ��f����Ĵ���
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
                        %��P������f����
                        P(beta + index_1, l_index) = sign(P(beta + index_2, l_index)) *...
                            sign(P(beta + index_1 + index_2, l_index)) *...
                            min(abs(P(beta + index_2, l_index)), abs(P(beta + index_1 + index_2, l_index)));
                    end
                end
                
            case N/2
                index_1 = lambda_offset(m);
                 %��llr��g����
                for beta = 0 : index_1 - 1
                    x_tmp = C(beta + index_1, 2 * l_index - 1);
                    P(beta + index_1, l_index) = (1 - 2 * x_tmp) *...
                        llr(beta + 1) + llr(beta + 1 + index_1);
                end
                for i_layer = m - 2 : -1 : 0
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    %��P��f����
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
                %��P������g����
                for beta = 0 : index_1 - 1
                    P(beta + index_1, l_index) = (1 - 2 *...
                        C(beta + index_1, 2 * l_index - 1)) *...
                        P(beta + index_2, lazy_copy(layer + 2, l_index)) +...
                        P(beta + index_1 + index_2, lazy_copy(layer + 2, l_index));
                end
                % ��P������f����                
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
    %����u�ǲ��ɿ���Ϣλʱ
    if frozen_bits(phi + 1) == 3%if now we decode an info bit
        %��Ϣ������Ҫ���п�¡������һ��SC������
        PM_pair = realmax * ones(2, L);%realmax����ָ���������������ܱ�ʾ���������ֵ��
        
        for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            %���㸴�Ƶ�PMֵ���ж����޼ӳͷ�ֵ
            if P(1, l_index) >= 0%LNk>=0ʱ
                %uk=0,����Ӳ�о���ԭPMֵ
                PM_pair(1, l_index) = PM(l_index);
                %uk=1,������Ӳ�о���+|LNk|
                PM_pair(2, l_index) = PM(l_index) + P(1, l_index);
            else%LNk<0ʱ
                %uk=0,������Ӳ�о���+|LNk|
                PM_pair(1, l_index) = PM(l_index) - P(1, l_index);
                %uk=1,����Ӳ�о���ԭPMֵ
                PM_pair(2, l_index) = PM(l_index);
            end
        end
        
        %%
        %���ܳ���L��PM_pairֵ����������Ҫɾѡȡ��С��L��PM_pairֵ
%         activepath
        middle = min(2 * sum(activepath), L);
        PM_sort = sort(PM_pair(:));%������PM_pairת����������������С��������
        PM_cv = PM_sort(middle);%�ɱ���·�������ֵ
        compare = PM_pair <= PM_cv; 
        
        %��¼��ɾ����·��
        kill_index = zeros(L, 1);%to record the index of the path that is killed
        %��ɾ��·�����ܸ���
        kill_cnt = 0;%the total number of killed path
        %��������������һ����ջ���
        
        for i = 1 : L
            %�ж��Ƿ��ж�Ӧ2��·����Ҫ��ɾ����SC������
            if (compare(1, i) == 0)&&(compare(2, i) == 0)
                activepath(i) = 0;%��ӦSC����������
                kill_cnt = kill_cnt + 1;%��ջ�������������ĸ���
                kill_index(kill_cnt) = i;%��������SC���������i����kill_index��
            end
        end
        
        %%
        %�жϱ���������SC������������·���Ƿ�ֱ�������������Ӧ��compareֵΪ1
        for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            %һ��*2��Ϊ����������·��,˵��2�Ǳ�����һ��·����1�Ǳ����ڶ���·��
            path_state = compare(1, l_index) * 2 + compare(2, l_index);
            
            %������path_state=0����Ϊ�������Ѿ����ǹ���
            switch path_state
                case 1%ֻ������2��·��
                    u_all(cnt_u, l_index) = 1;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 1;
                    PM(l_index) = PM_pair(2, l_index);
                case 2%ֻ������1��·��
                    u_all(cnt_u, l_index) = 0;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 0;
                    PM(l_index) = PM_pair(1, l_index);
                case 3%����2��·�����赯��һ���������洢uk=1
                    index = kill_index(kill_cnt);%����ջ�����һ����������������indexΪ�����
                    kill_cnt = kill_cnt - 1;%ʣ������SC����������
                    activepath(index) = 1;%��������������������
                    
                    %lazy copy����
                    %��ԭ������LC��һ�и��Ƹ�������������Ӧ����
                    lazy_copy(:, index) = lazy_copy(:, l_index);
                    u_all(:, index) = u_all(:, l_index);%u���Ƶ����������������Ӧu��
                    u_all(cnt_u, l_index) = 0;%ԭ����������u=0
                    u_all(cnt_u, index) = 1;%����������������u=1
                    C(1, 2 * l_index - 1 + phi_mod_2) = 0;%��ӦC��u=0
                    C(1, 2 * index - 1 + phi_mod_2) = 1;%��ӦC��u=1
                    PM(l_index) = PM_pair(1, l_index);%��PM_pairֵ�����ӦPM��
                    PM(index) = PM_pair(2, l_index);
            end
        end
        cnt_u = cnt_u + 1;
     
       %%
      % ���Ƕ���λʱ
    else if frozen_bits(phi + 1) == 1||frozen_bits(phi + 1) == 11%frozen bit operation
        for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            if frozen_bits(phi + 1) == 11
            if P(1, l_index) < 0%LNk<0��uk=0,������Ӳ�о���ӳͷ�
                PM(l_index) = PM(l_index) - P(1, l_index);%PM+|LNk|
            end
            end
            if phi_mod_2 == 0%�ж�u�±��Ƿ�Ϊ����
                C(1, 2 * l_index - 1) = 0;%����������uk=0����C[1]�����һλ
            else
                C(1, 2 * l_index) = 0;%��ż������uk=0����C[1]���ұ�һλ
            end 
        end
     else if frozen_bits(phi + 1) == 2%У����� 
             active_index=[];%��¼��ͨ����żУ����ص�·����
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
%                if p_2(u_p,l_index)~=p_1(u_p)%����ͨ����żУ���·��ɾ����������ͨ����żУ���·��
               if p_2(u_p,l_index)~=p_3%��ͨ����żУ���·���ӳͷ�ֵ
                    PM(l_index) = PM(l_index) +abs(P(1, l_index));
%                      activepath(l_index) = 0;%��ӦSC����������
                  continue
               end
               active_index=[active_index,l_index];%��¼ͨ����żУ���·����� 
             end
             active_length(u_p)=length(active_index);
             
           %����Lp��ͨ����żУ���·��     
           if length(active_index)>=Lp %�ж���·��ͨ����żУ��
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
                  if ismember(l_index,PM_min_index)==0   %����·��������С��Lp��ͨ��·��
                      activepath(l_index) = 0;%��ӦSC����������
%                       kill_cnt = kill_cnt + 1;%��ջ�������������ĸ���
%                       kill_index(kill_cnt) = l_index;%��������SC���������l_index����kill_index��
                      continue
%                         lazy_copy(:, l_index) = lazy_copy(:, PM_min_index);
%                         u_all(:,l_index)=u_all(:,PM_min_index);
                  end
                  if phi_mod_2 == 0%�ж�u�±��Ƿ�Ϊ����
                     C(1, 2 * l_index - 1) = p_2(u_p,l_index);%����������uk=0����C[1]�����һλ
                    else
                    C(1, 2 * l_index) = p_2(u_p,l_index);%��ż������uk=0����C[1]���ұ�һλ
                  end  
                     
                 
                  end
                  
              
           elseif length(active_index)==0%��żУ���δͨ��,����PM��С��Lp��
%                index=3
%                activepath
           
                 row=find(activepath==1);%�ڴ���·���н�PM����
%                  PM
%                   PM(row)
                 
     
               [~, min_index_11]=sort(PM(row));%min_index_1Ϊ���·�������
               min_index_1=row(min_index_11(1:Lp));
%                min_index_1
               for l_index = 1 : L
                    if activepath(l_index) == 0
                continue;
                    end
                     if ismember(l_index,min_index_1)==0
                     activepath(l_index) = 0;%��ӦSC����������
%                       kill_cnt = kill_cnt + 1;%��ջ�������������ĸ���
%                       kill_index(kill_cnt) = l_index;%��������SC���������l_index����kill_index��
                      continue
                     end
%                      lazy_copy(:, l_index) = lazy_copy(:,min_index_1);
%                         u_all(:,l_index)=u_all(:,min_index_1);
                if phi_mod_2 == 0%�ж�u�±��Ƿ�Ϊ����
                     C(1, 2 * l_index - 1) = p_2(u_p,l_index);%����������uk=0����C[1]�����һλ
                    else
                    C(1, 2 * l_index) = p_2(u_p,l_index);%��ż������uk=0����C[1]���ұ�һλ
                end 
                    
               end
                
              
                 
              
           else%  1��ͨ����żУ��·������Lp
%                for l_index = 1 : L
%                     if activepath(l_index) == 0
%                 continue;
%                     end
                    len = length(active_index);
                   min_index=find(activepath==1);%�ڴ���·���н�PM����
                    [ia, ~] = setdiff(min_index, active_index);%�ڴ��·���н�ͨ����żУ���·��ȥ��
%                     res = min_index(sort(ia));
                     [~, min_index_]=sort(PM(ia));%��ʣ��·������PM��С��Lp-len��·������
               min_index_12=ia(min_index_(1:Lp-len));
               min_index_12=min_index_12';
               real_index = [active_index min_index_12];
                for l_index=1:L
                  if ismember(l_index, real_index)==0   %����·��������С��Lp��ͨ��·��
                      activepath(l_index) = 0;%��ӦSC����������
%                       kill_cnt = kill_cnt + 1;%��ջ�������������ĸ���
%                       kill_index(kill_cnt) = l_index;%��������SC���������l_index����kill_index��
                      continue
%                         lazy_copy(:, l_index) = lazy_copy(:, PM_min_index);
%                         u_all(:,l_index)=u_all(:,PM_min_index);
                  end
%                    lazy_copy(:, l_index) = lazy_copy(:, active_index);
%                         u_all(:,l_index)=u_all(:,active_index);
                  if phi_mod_2 == 0%�ж�u�±��Ƿ�Ϊ����
                     C(1, 2 * l_index - 1) = p_2(u_p,l_index);%����������uk=0����C[1]�����һλ
                    else
                C(1, 2 * l_index) = p_2(u_p,l_index);%��ż������uk=0����C[1]���ұ�һλ
                  end 
           
               end
           
                
                  
                
                  
           end
        
              
              
              u_p=u_p+1;
          else %�Ͽɿ���Ϣλ������·�����ţ�SC����
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
    %partial-sum return����ֵx
    for l_index = 1 : L
        if activepath(l_index) == 0
            continue
        end
        %�ж�u�±��Ƿ���ż����Ϊż��ʱ��Ҫ����C����
        if (phi_mod_2  == 1) && (phi ~= N - 1)
            layer = bit_layer_vec(phi + 1);
            %��ӦSC��������C�����ڶ��и�ֵ
            for i_layer = 0 : layer - 1
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : 2 * index_1 - 1
                    C(beta + index_1, 2 * l_index) = mod(C(beta, 2 *  lazy_copy(i_layer + 1, l_index) - 1) + C(beta, 2 * l_index), 2);%Left Column lazy copy
                    C(beta + index_2, 2 * l_index) = C(beta, 2 * l_index);   
                end
            end
            %��ӦSC��������C������һ�и�ֵ��C���½���
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
    %lazy copy����
    %C���������ϣ������жϼ���·����Ӧ��LC�е���
    if phi < N - 1
        %��һ�����ض�����չ�����λ��m������0��m+1
        for i_layer = 1 : llr_layer_vec(phi + 2) + 1
            %��Ӧǰm+1��Ԫ���޸�Ϊ��Ӧ���������������
            for l_index = 1 : L
                lazy_copy(i_layer, l_index) = l_index;
            end
        end
    end
end

%%
%path selection.
%ѡ�����·������PMֵ��С��·��
% [~, min_index]=min(PM);%min_indexΪ���·�������
% polar_info_esti = u_all(:,min_index);%���·������Ϣλuȡ����ֻ����Ϣλ��������,����Դ�е���Ϣλ�Ƚϼ��ɣ�
final_index=find(activepath==1);
if length(final_index) ~=1
   [~, final_index1]=min(PM(final_index));
   final_index = final_index(final_index1);
end
% final_index=final_index+1;
polar_info_esti = u_all(:,final_index);
end
