function Ls = regularization(Km)

% *****************************
% JSheng
% *****************************
Ds0 = zeros(Km, Km);
for il = 1 : Km
    for jl = 1 : Km
        if 2 * (il - 1) == jl - 1
           Ds0(il, jl) = 1;
        end   
    end
end
Ds = kron(sparse(eye(Km)), sparse(Ds0));    % OK  
Qs_1d_1 = zeros(1, Km + 1);
Qs_1d_1(1, 1) = -2;
Qs_1d_1(1, 2) = 1;
Qs_1d_1(1, Km + 1) = 1;
Qs_1d_2 = zeros(1, Km + 1);
Qs_1d_2(1, Km + 1) = -2;
Qs_1d_2(1, Km) = 1;
Qs_1d_2(1, 1) = 1;
for il = 1 : Km * Km - Km
    Qs(il, il : il + Km) = sparse(Qs_1d_1(1, 1 : Km + 1));
end
for il = Km * Km - Km + 1 : Km * Km
    Qs(il, il - Km : il) = sparse(Qs_1d_2(1, 1 : Km + 1));
end
Ls = Qs' * Qs;   
