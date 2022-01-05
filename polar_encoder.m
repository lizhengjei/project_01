function x= polar_encoder( u,lambda_offset,layer_vec )
%仅用for循环的代码进行编码，需要参数lambda_offset，layer_vec
%encoding：x = u * Fn, 未考虑置换！！！！

%lambda_offset表示“分段”用的向量，长度是n+1
%lambda_offset元素取值2.^(0:log2(N))
%调用get_layer(N)函数。

N =length(u);
m = log2(N);
x_internal_value =zeros(2*N -  1,1);%存储中间计算结果的向量
x_internal_value(end-N+1:end)=u;
x =zeros(N,1);%最终输出

%计算x1
for i_layer = m -  1  : - 1  : 0
    index_1 = lambda_offset(i_layer +  1);%定位
    index_2 = lambda_offset(i_layer +  2);
    for beta = index_1 : index_2 -  1
        x_internal_value(beta)=...
            x_internal_value(index_2+beta-index_1)+...
             x_internal_value(index_2+beta);
    end
end
x(1)=x_internal_value(1);

%计算x1~xN
for phi=1:N-1
    layer=layer_vec(phi+1);
    %直接得到位
    index_1 = lambda_offset(layer +  1);%定位
    index_2 = lambda_offset(layer +  2);
     for beta = index_1 : index_2 -  1
        x_internal_value(beta)=...
            x_internal_value(index_2+beta);
     end
     %相加得到位
    for i_layer = layer -  1  : - 1  : 0
    index_1 = lambda_offset(i_layer +  1);%定位
    index_2 = lambda_offset(i_layer +  2);
        for beta = index_1 : index_2 -  1
            x_internal_value(beta)=...
             x_internal_value(index_2+beta-index_1)+...
               x_internal_value(index_2+beta);
        end
    end
    x(phi +  1)= x_internal_value(1);
end

% 前面相加都是数值直接相加，需转换成模2加
x = mod(x,2);
end

