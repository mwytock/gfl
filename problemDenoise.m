function Y = problemDenoise(name,w)
Y = double(imresize(imread(name), [w nan]))/256;
Y = max(min(Y + 0.1*randn(size(Y)),1),0);
Y = permute(Y, [3 1 2]);