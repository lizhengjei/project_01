function x= polar_encoder( u,lambda_offset,layer_vec )
%����forѭ���Ĵ�����б��룬��Ҫ����lambda_offset��layer_vec
%encoding��x = u * Fn, δ�����û���������

%lambda_offset��ʾ���ֶΡ��õ�������������n+1
%lambda_offsetԪ��ȡֵ2.^(0:log2(N))
%����get_layer(N)������

N =length(u);
m = log2(N);
x_internal_value =zeros(2*N -  1,1);%�洢�м������������
x_internal_value(end-N+1:end)=u;
x =zeros(N,1);%�������

%����x1
for i_layer = m -  1  : - 1  : 0
    index_1 = lambda_offset(i_layer +  1);%��λ
    index_2 = lambda_offset(i_layer +  2);
    for beta = index_1 : index_2 -  1
        x_internal_value(beta)=...
            x_internal_value(index_2+beta-index_1)+...
             x_internal_value(index_2+beta);
    end
end
x(1)=x_internal_value(1);

%����x1~xN
for phi=1:N-1
    layer=layer_vec(phi+1);
    %ֱ�ӵõ�λ
    index_1 = lambda_offset(layer +  1);%��λ
    index_2 = lambda_offset(layer +  2);
     for beta = index_1 : index_2 -  1
        x_internal_value(beta)=...
            x_internal_value(index_2+beta);
     end
     %��ӵõ�λ
    for i_layer = layer -  1  : - 1  : 0
    index_1 = lambda_offset(i_layer +  1);%��λ
    index_2 = lambda_offset(i_layer +  2);
        for beta = index_1 : index_2 -  1
            x_internal_value(beta)=...
             x_internal_value(index_2+beta-index_1)+...
               x_internal_value(index_2+beta);
        end
    end
    x(phi +  1)= x_internal_value(1);
end

% ǰ����Ӷ�����ֱֵ����ӣ���ת����ģ2��
x = mod(x,2);
end

