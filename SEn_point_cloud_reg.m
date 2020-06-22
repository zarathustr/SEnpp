clear all
close all
clc

len = 100000;
dim = 10;
d = 1e6;
num = 1;
err = 1e-4;

ll = [0, 0];

for j = 1 : num
    b = zeros(dim, len);
    r = randn(dim, len);
    T = randn(dim, 1);

    B = randn(dim, dim);
    R = orthonormalize(B);

    for i = 1 : len
        b(:, i) = R * r(:, i) + T + err * randn(dim, 1);
    end


    M = zeros(dim + 1, dim + 1);
    for i = 1 : len
        b_ = b(:, i);
        r_ = r(:, i);
    
        b_aug = [b_; d];
        r_aug = [r_; d];
    
        b_aug = b_aug ./ norm(b_aug);
        r_aug = r_aug ./ norm(r_aug);
    
        M = M + 1 / len * b_aug * r_aug';
    end

    [RR, t] = SOnp1_SEn(orthonormalize(M), d);


    L = 0;
    for i = 1 : len
        L = L + 1 / len * norm(b(:, i) - RR * r(:, i) - t)^2;
    end


    mean_b = zeros(dim, 1);
    mean_r = zeros(dim, 1);
    for i = 1 : len
        mean_b = mean_b + 1 / len * b(:, i);
        mean_r = mean_r + 1 / len * r(:, i);
    end

    M = zeros(dim, dim);
    for i = 1 : len
        b_ = b(:, i) - mean_b;
        r_ = r(:, i) - mean_b;
    
        M = M + 1 / len * b_ * r_';
    end
    R_ = orthonormalize(M);
    t_ = mean_b - R_ * mean_r;


    LL = 0;
    for i = 1 : len
        LL = LL + 1 / len * norm(b(:, i) - R_ * r(:, i) - t_)^2;
    end
    
    
    ll = ll + 1 / num * [L, LL];
end

ll
