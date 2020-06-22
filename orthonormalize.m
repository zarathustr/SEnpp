function R = orthonormalize(A)
ss = size(A);
dim = ss(1);
[u, ~, v] = svd(A);
s = eye(dim, dim);
s(dim, dim) = det(u * v);
R = u * s * v';
end