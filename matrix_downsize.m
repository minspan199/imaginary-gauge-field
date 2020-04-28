
%% input matrix A

[width, length] = size(A);

if width ~= length
    A = convert2Square(A);
end

time = 40;
[N,K] = size(A);
clear B;
for ii = 1:time: N
    for jj = 1:time:K
        B(floor(ii/time) + 1, floor(jj/time)+1) = A(ii, jj);
    end
end

[M,~] = size(B);

imwrite(B, 'LetterACompressed2.png')
figure;
imagesc(B)

function SMatrix = convert2Square(A)

    [width, length] = size(A);
    N = width + length;
    left = floor(length / 2);
    right = width - left;
    B = [zeros(width, left) A zeros(width, right)];
    top = floor(width / 2);
    bottom = length - top;
    SMatrix = [zeros(top, N); B; zeros(bottom, N)];
end

    
    