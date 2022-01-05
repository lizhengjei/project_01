function bit_layer_vec = get_bit_layer(N)
bit_layer_vec = zeros(N, 1);
for phi = 0 : N - 1
    psi = floor(phi/2);
    layer = 0;
    while(mod(psi, 2) == 1)
        psi = floor(psi/2);
        layer = layer + 1;
    end
    bit_layer_vec(phi + 1) = layer;
end
end