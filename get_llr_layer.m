function llr_layer_vec = get_llr_layer(N)
llr_layer_vec = zeros(N , 1);
for phi = 1 : N - 1
    psi = phi;
    layer = 0;
    while(mod(psi, 2) == 0)
        psi = floor(psi/2);
        layer = layer + 1;
    end
    llr_layer_vec(phi + 1) = layer;
end
end