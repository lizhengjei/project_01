clc
clear
close all
%%N=256,R=1/2,L=16
%%校验位在固定位上
%%译码采用SCL译码，保留所有路径，奇偶校验不通过加惩罚值
% 信噪比0.5，仿真次数15000,误码率0.1314472，误块率0.4854369
% 信噪比1.0，仿真次数15000,误码率0.0622210，误块率0.2747253
% 信噪比1.5，仿真次数15000,误码率0.0183833，误块率0.0933333
% 信噪比2.0，仿真次数15000,误码率0.0050229，误块率0.0341333
% 信噪比2.5，仿真次数15000,误码率0.0011667，误块率0.0112000
% 信噪比3.0，仿真次数15000,误码率0.0002031，误块率0.0018000
%%译码采用SCL，只保留奇偶校验通过的路径
% 信噪比0.5，仿真次数15000,误码率0.1128095，误块率0.4830918
% 信噪比1.0，仿真次数15000,误码率0.0631510，误块率0.2873563
% 信噪比1.5，仿真次数15000,误码率0.0168276，误块率0.0888000
% 信噪比2.0，仿真次数15000,误码率0.0048156，误块率0.0330000
% 信噪比2.5，仿真次数15000,误码率0.0011365，误块率0.0108667
% 信噪比3.0，仿真次数15000,误码率0.0001901，误块率0.0017333
%%较不可靠信息位采取SCL，其余信息位SC译码，只保留奇偶校验通过的路径
% 信噪比1.0，仿真次数15000,误码率0.1603936，误块率0.5524862
% 信噪比1.5，仿真次数15000,误码率0.0876202，误块率0.3076923
% 信噪比2.0，仿真次数15000,误码率0.0278374，误块率0.1108647
% 信噪比2.5，仿真次数15000,误码率0.0108865，误块率0.0450667
% 信噪比3.0，仿真次数15000,误码率0.0029021，误块率0.0125333
% 信噪比3.5，仿真次数15000,误码率0.0005302，误块率0.0025333
%%
%校验位在信息位上，较不可靠信息位采取SCL，其余信息位和冻结比特SC译码，保留所有路径，加惩罚值m=60,M=5
% 信噪比1.0，仿真次数15000,误码率0.0875884，误块率0.3731343
% 信噪比1.5，仿真次数15000,误码率0.0393264，误块率0.1536098
% 信噪比2.0，仿真次数15000,误码率0.0133268，误块率0.0644667
% 信噪比2.5，仿真次数15000,误码率0.0039584，误块率0.0208667
% 信噪比3.0，仿真次数15000,误码率0.0005564，误块率0.0036000
% 信噪比3.5，仿真次数15000,误码率0.0000391，误块率0.0002667
%m=50
% 信噪比1.0，仿真次数15000,误码率0.0911756，误块率0.3610108
% 信噪比1.5，仿真次数15000,误码率0.0403367，误块率0.1572327
% 信噪比2.0，仿真次数15000,误码率0.0138421，误块率0.0670000
% 信噪比2.5，仿真次数15000,误码率0.0043429，误块率0.0226667
% 信噪比3.0，仿真次数15000,误码率0.0009368，误块率0.0049333
% 信噪比3.5，仿真次数15000,误码率0.0001188，误块率0.0007333
%L=8
% 信噪比1.0，仿真次数15000,误码率0.0947311，误块率0.3816794
% 信噪比1.5，仿真次数15000,误码率0.0367224，误块率0.1545595
% 信噪比2.0，仿真次数15000,误码率0.0130115，误块率0.0642667
% 信噪比2.5，仿真次数15000,误码率0.0039353，误块率0.0213333
% 信噪比3.0，仿真次数15000,误码率0.0007830，误块率0.0044667
% 信噪比3.5，仿真次数15000,误码率0.0000787，误块率0.0003333
%%只保留通过奇偶校验的路径，m=60,M=4
% 信噪比1.0，仿真次数15000,误码率0.0941533，误块率0.3731343
% 信噪比1.5，仿真次数15000,误码率0.0445035，误块率0.1841621
% 信噪比2.0，仿真次数15000,误码率0.0148471，误块率0.0679333
% 信噪比2.5，仿真次数15000,误码率0.0035724，误块率0.0195333
% 信噪比3.0，仿真次数15000,误码率0.0007880，误块率0.0042667
% 信噪比3.5，仿真次数15000,误码率0.0001489，误块率0.0010667

%%
n=8;
N=2^n;
 M=5;%较可靠固定比特放置奇偶校验位
%  R=1/2;
R=(floor(112/256*N)+M)/N;
S=N*R; %信息位码长,包括奇偶校验
% M=5;%奇偶校验位
% m=55;%较不可靠信息比特

Se=S-M;%真正信息位
m=floor(0.42*Se);%包含SCL的信息位和奇偶校验位
F=N-S; %冻结位码长
L=16;
lambda_offset=2.^(0:log2(N));
 bit_layer_vec = get_bit_layer(N);
 llr_layer_vec = get_llr_layer(N);
 frozen_bits=zeros(1,N);%冻结比特标志位初始化
 u=zeros(1,N);%信源比特初始化
% x_outer=zeros(1,S);
% real_info_index=zeros(1,K);
 experimentnumber=15000;
 E_N=0.5:0.25:3;%信噪比，单位为dB
 IE=length(E_N);
 BER=zeros(1,IE);
FER=zeros(1,IE);
for i=1:IE
 sigma=(1/sqrt(2*R))*10.^((-E_N(i))/20);
 [g,~] = GA(sigma, N);
 pe=1/2*erfc(sqrt(g)/2);
%  [G_in_order,index] = sort(g,'descend');%将高斯子信道错误概率按从大到小排列，选取大的做信息位
 [~,index]=sort(pe);%从小到大

 signal_index =  sort(index( 1:S )) ;   %前S位作为信息位
 frozen_index =sort( index( S+1:end ));   %后面的作为冻结位
%  index1 = sort(index(S-m+1:S));
%  u_frozen_index = index1(m-M+1:end);
%  u_frozen_index = sort(index1(getIndex(1, m, M)));%平均分布，奇偶校验比特
%  u_signal_index = setdiff(index1,u_frozen_index);
% u_signal_index=sort(index(S-m+1:S-M+1));%较不可靠信息位，用来被校验
% u_frozen_index=sort(index(S-M+2:S));%最不可靠信息位放置校验比特
u_signal_index=sort(index(S-m+M:S));%较不可靠信息位，用来被校验，64-5
u_frozen_index=sort(index(S-m+1:S-m+M-1));%最可靠信息位放置M-1个校验比特,4
% u_frozen_index=sort(index(S-m-1:S-m+M-1));%n=8时
% u_frozen_index=u_frozen_index(1:end-2);%n=8时
u_frozen_index(M)=u_signal_index(end);%最后一位校验比特在最后一个,1
u_frozen_index=sort(u_frozen_index);
u_signal_index=u_signal_index(1:end-1);%信息比特除去最后一个校验比特，较不可靠信息位
u_info_index=setdiff(signal_index,u_frozen_index);%setdiff表示从A数组中去除B数组所剩的数，其中B数组包含于A数组

  u_x=sort([u_signal_index u_frozen_index]);%较不可靠信息位和校验位结合形成外码块
  for vv=1:length(u_frozen_index)
  u_x_frozen_index(vv)=find(u_x==u_frozen_index(vv));%校验比特在外码码字中的序号
  end

%在u_x每个外码段之间的冻结比特加惩罚值
u_x_frozen_index1=[];
for jj=1:length(frozen_index)
    for ii=1:M
        if ii==1&&frozen_index(jj)>u_x(1)&&frozen_index(jj)<u_frozen_index(1)
            u_x_frozen_index1=[u_x_frozen_index1,frozen_index(jj)];
        else if ii~=1&&frozen_index(jj)>u_x(u_x_frozen_index(ii-1)+1)&&frozen_index(jj)<u_frozen_index(ii)
            u_x_frozen_index1=[u_x_frozen_index1,frozen_index(jj)];
            end
         end
    end
    
end

        
            
for j=1:experimentnumber
signal=randi([0 1],1,Se);%生成1*S的0、1随机的矩阵,信息比特

u(1,frozen_index)=0;
u(1,u_info_index)=signal;
%计算校验方程以及校验比特
c_info=cell(1,M);%元胞数组，记录校验方程在u_signal_index中的序号 
for kk=1:M
    index_=0;
    index_=find(u_x==u_frozen_index(kk));
   c_info{kk}=[1:1:index_-kk]; 
end
c_check=cell(1,M);%信息比特根据概率需要校验的序号（在u_signal_index中）
for jj=1:M
%     c=[];
    c=calcu(c_info{jj});
    c_check{jj}=c;
end
%检验方程转换为在u_info_index中(signal)的序号
c_check_real=cell(1,M);
for kkk=1:M
    for jjj=1:length(c_check{kkk})
        c_check_real{kkk}(jjj)=find(u_info_index==u_signal_index(c_check{kkk}(jjj)));
    end
end
%计算校验比特
for iii=1:M
    x_p=0;
   for s=1:length(c_check_real{iii})
       x_p=x_p+signal(c_check_real{iii}(s));
   end
   u(1,u_frozen_index(iii))=mod(x_p,2);
end

 x= polar_encoder( u,lambda_offset,llr_layer_vec );%编码，不含置换
  x=x';
s=1-2*x; %BPSK调制
n=0+sigma.*randn(1,N);%高斯噪声，均值为0，方差为sigma^2，是*，不是^！！！！！！！！！！！！！！！！
y=s+n;
llr=2*y/(sigma^2);
frozen_bits(1,frozen_index)=1;%为1说明是冻结比特
frozen_bits(1,u_x_frozen_index1)=11;%校验区间内的冻结比特加惩罚值
frozen_bits(1,signal_index)=0;%为0是信息比特，较可靠信息位，SC译码
frozen_bits(1,u_frozen_index)=2;%校验比特，
frozen_bits(1,u_signal_index)=3;%较不可靠信息位，SCL译码并奇偶校验
%  p_1=u(1,u_frozen_index);%校验比特
[polar_info_esti,u_all,active_length] = SCL_decoder(llr, L, S, frozen_bits,  lambda_offset, llr_layer_vec, bit_layer_vec,c_check_real,M);


% for vv=1:S-M
%     real_info_index(vv)=find(signal_index==u_info_index(vv));
%    
% end
% polar_info_esti1=polar_info_esti(real_info_index);
errbits=length(find(polar_info_esti'~=signal));
if errbits~=0
            errblock=1;
        else
            errblock=0;
end
        BER(i)=BER(i)+errbits;
        FER(i)=FER(i)+errblock;%统计误块个数
        if(FER(i)>=100&&j<=1000)
            break;
        end
end
BER(i)=BER(i)/(Se*j);
FER(i)=FER(i)/j;
fprintf('信噪比%.1f，仿真次数%.f,误码率%.7f，误块率%.7f\n',E_N(i),experimentnumber,BER(i),FER(i));
end
% figure(1)
% semilogy(E_N,BER,'LineWidth',2,'MarkerSize',8);
% grid on;
% title('误码率曲线图','fontsize',12,'fontweight','b');
% xlabel('channelstate','fontsize',12,'fontweight','b');
% ylabel('ber','fontsize',12,'fontweight','b');
% figure(2)
% semilogy(E_N,FER,'LineWidth',2,'MarkerSize',8);
% grid on;
% title('误块率曲线图','fontsize',12,'fontweight','b');
% xlabel('channelstate','fontsize',12,'fontweight','b');
% ylabel('fer','fontsize',12,'fontweight','b');








