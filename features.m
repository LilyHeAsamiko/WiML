function S = features(Matrix)
[m, n] = size(Matrix);
S= 0.00000001*Matrix;
S=sum(sum(Matrix))/(m*n);
end