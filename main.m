clc
clear
close all
%%N=256,R=1/2,L=16
%%У��λ�ڹ̶�λ��
%%�������SCL���룬��������·������żУ�鲻ͨ���ӳͷ�ֵ
% �����0.5���������15000,������0.1314472�������0.4854369
% �����1.0���������15000,������0.0622210�������0.2747253
% �����1.5���������15000,������0.0183833�������0.0933333
% �����2.0���������15000,������0.0050229�������0.0341333
% �����2.5���������15000,������0.0011667�������0.0112000
% �����3.0���������15000,������0.0002031�������0.0018000
%%�������SCL��ֻ������żУ��ͨ����·��
% �����0.5���������15000,������0.1128095�������0.4830918
% �����1.0���������15000,������0.0631510�������0.2873563
% �����1.5���������15000,������0.0168276�������0.0888000
% �����2.0���������15000,������0.0048156�������0.0330000
% �����2.5���������15000,������0.0011365�������0.0108667
% �����3.0���������15000,������0.0001901�������0.0017333
%%�ϲ��ɿ���Ϣλ��ȡSCL��������ϢλSC���룬ֻ������żУ��ͨ����·��
% �����1.0���������15000,������0.1603936�������0.5524862
% �����1.5���������15000,������0.0876202�������0.3076923
% �����2.0���������15000,������0.0278374�������0.1108647
% �����2.5���������15000,������0.0108865�������0.0450667
% �����3.0���������15000,������0.0029021�������0.0125333
% �����3.5���������15000,������0.0005302�������0.0025333
%%
%У��λ����Ϣλ�ϣ��ϲ��ɿ���Ϣλ��ȡSCL��������Ϣλ�Ͷ������SC���룬��������·�����ӳͷ�ֵm=60,M=5
% �����1.0���������15000,������0.0875884�������0.3731343
% �����1.5���������15000,������0.0393264�������0.1536098
% �����2.0���������15000,������0.0133268�������0.0644667
% �����2.5���������15000,������0.0039584�������0.0208667
% �����3.0���������15000,������0.0005564�������0.0036000
% �����3.5���������15000,������0.0000391�������0.0002667
%m=50
% �����1.0���������15000,������0.0911756�������0.3610108
% �����1.5���������15000,������0.0403367�������0.1572327
% �����2.0���������15000,������0.0138421�������0.0670000
% �����2.5���������15000,������0.0043429�������0.0226667
% �����3.0���������15000,������0.0009368�������0.0049333
% �����3.5���������15000,������0.0001188�������0.0007333
%L=8
% �����1.0���������15000,������0.0947311�������0.3816794
% �����1.5���������15000,������0.0367224�������0.1545595
% �����2.0���������15000,������0.0130115�������0.0642667
% �����2.5���������15000,������0.0039353�������0.0213333
% �����3.0���������15000,������0.0007830�������0.0044667
% �����3.5���������15000,������0.0000787�������0.0003333
%%ֻ����ͨ����żУ���·����m=60,M=4
% �����1.0���������15000,������0.0941533�������0.3731343
% �����1.5���������15000,������0.0445035�������0.1841621
% �����2.0���������15000,������0.0148471�������0.0679333
% �����2.5���������15000,������0.0035724�������0.0195333
% �����3.0���������15000,������0.0007880�������0.0042667
% �����3.5���������15000,������0.0001489�������0.0010667

%%
n=8;
N=2^n;
 M=5;%�Ͽɿ��̶����ط�����żУ��λ
%  R=1/2;
R=(floor(112/256*N)+M)/N;
S=N*R; %��Ϣλ�볤,������żУ��
% M=5;%��żУ��λ
% m=55;%�ϲ��ɿ���Ϣ����

Se=S-M;%������Ϣλ
m=floor(0.42*Se);%����SCL����Ϣλ����żУ��λ
F=N-S; %����λ�볤
L=16;
lambda_offset=2.^(0:log2(N));
 bit_layer_vec = get_bit_layer(N);
 llr_layer_vec = get_llr_layer(N);
 frozen_bits=zeros(1,N);%������ر�־λ��ʼ��
 u=zeros(1,N);%��Դ���س�ʼ��
% x_outer=zeros(1,S);
% real_info_index=zeros(1,K);
 experimentnumber=15000;
 E_N=0.5:0.25:3;%����ȣ���λΪdB
 IE=length(E_N);
 BER=zeros(1,IE);
FER=zeros(1,IE);
for i=1:IE
 sigma=(1/sqrt(2*R))*10.^((-E_N(i))/20);
 [g,~] = GA(sigma, N);
 pe=1/2*erfc(sqrt(g)/2);
%  [G_in_order,index] = sort(g,'descend');%����˹���ŵ�������ʰ��Ӵ�С���У�ѡȡ�������Ϣλ
 [~,index]=sort(pe);%��С����

 signal_index =  sort(index( 1:S )) ;   %ǰSλ��Ϊ��Ϣλ
 frozen_index =sort( index( S+1:end ));   %�������Ϊ����λ
%  index1 = sort(index(S-m+1:S));
%  u_frozen_index = index1(m-M+1:end);
%  u_frozen_index = sort(index1(getIndex(1, m, M)));%ƽ���ֲ�����żУ�����
%  u_signal_index = setdiff(index1,u_frozen_index);
% u_signal_index=sort(index(S-m+1:S-M+1));%�ϲ��ɿ���Ϣλ��������У��
% u_frozen_index=sort(index(S-M+2:S));%��ɿ���Ϣλ����У�����
u_signal_index=sort(index(S-m+M:S));%�ϲ��ɿ���Ϣλ��������У�飬64-5
u_frozen_index=sort(index(S-m+1:S-m+M-1));%��ɿ���Ϣλ����M-1��У�����,4
% u_frozen_index=sort(index(S-m-1:S-m+M-1));%n=8ʱ
% u_frozen_index=u_frozen_index(1:end-2);%n=8ʱ
u_frozen_index(M)=u_signal_index(end);%���һλУ����������һ��,1
u_frozen_index=sort(u_frozen_index);
u_signal_index=u_signal_index(1:end-1);%��Ϣ���س�ȥ���һ��У����أ��ϲ��ɿ���Ϣλ
u_info_index=setdiff(signal_index,u_frozen_index);%setdiff��ʾ��A������ȥ��B������ʣ����������B���������A����

  u_x=sort([u_signal_index u_frozen_index]);%�ϲ��ɿ���Ϣλ��У��λ����γ������
  for vv=1:length(u_frozen_index)
  u_x_frozen_index(vv)=find(u_x==u_frozen_index(vv));%У����������������е����
  end

%��u_xÿ�������֮��Ķ�����ؼӳͷ�ֵ
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
signal=randi([0 1],1,Se);%����1*S��0��1����ľ���,��Ϣ����

u(1,frozen_index)=0;
u(1,u_info_index)=signal;
%����У�鷽���Լ�У�����
c_info=cell(1,M);%Ԫ�����飬��¼У�鷽����u_signal_index�е���� 
for kk=1:M
    index_=0;
    index_=find(u_x==u_frozen_index(kk));
   c_info{kk}=[1:1:index_-kk]; 
end
c_check=cell(1,M);%��Ϣ���ظ��ݸ�����ҪУ�����ţ���u_signal_index�У�
for jj=1:M
%     c=[];
    c=calcu(c_info{jj});
    c_check{jj}=c;
end
%���鷽��ת��Ϊ��u_info_index��(signal)�����
c_check_real=cell(1,M);
for kkk=1:M
    for jjj=1:length(c_check{kkk})
        c_check_real{kkk}(jjj)=find(u_info_index==u_signal_index(c_check{kkk}(jjj)));
    end
end
%����У�����
for iii=1:M
    x_p=0;
   for s=1:length(c_check_real{iii})
       x_p=x_p+signal(c_check_real{iii}(s));
   end
   u(1,u_frozen_index(iii))=mod(x_p,2);
end

 x= polar_encoder( u,lambda_offset,llr_layer_vec );%���룬�����û�
  x=x';
s=1-2*x; %BPSK����
n=0+sigma.*randn(1,N);%��˹��������ֵΪ0������Ϊsigma^2����*������^��������������������������������
y=s+n;
llr=2*y/(sigma^2);
frozen_bits(1,frozen_index)=1;%Ϊ1˵���Ƕ������
frozen_bits(1,u_x_frozen_index1)=11;%У�������ڵĶ�����ؼӳͷ�ֵ
frozen_bits(1,signal_index)=0;%Ϊ0����Ϣ���أ��Ͽɿ���Ϣλ��SC����
frozen_bits(1,u_frozen_index)=2;%У����أ�
frozen_bits(1,u_signal_index)=3;%�ϲ��ɿ���Ϣλ��SCL���벢��żУ��
%  p_1=u(1,u_frozen_index);%У�����
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
        FER(i)=FER(i)+errblock;%ͳ��������
        if(FER(i)>=100&&j<=1000)
            break;
        end
end
BER(i)=BER(i)/(Se*j);
FER(i)=FER(i)/j;
fprintf('�����%.1f���������%.f,������%.7f�������%.7f\n',E_N(i),experimentnumber,BER(i),FER(i));
end
% figure(1)
% semilogy(E_N,BER,'LineWidth',2,'MarkerSize',8);
% grid on;
% title('����������ͼ','fontsize',12,'fontweight','b');
% xlabel('channelstate','fontsize',12,'fontweight','b');
% ylabel('ber','fontsize',12,'fontweight','b');
% figure(2)
% semilogy(E_N,FER,'LineWidth',2,'MarkerSize',8);
% grid on;
% title('���������ͼ','fontsize',12,'fontweight','b');
% xlabel('channelstate','fontsize',12,'fontweight','b');
% ylabel('fer','fontsize',12,'fontweight','b');








