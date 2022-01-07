% mean square error
function y = mse(I)
y = mean(I(:).*I(:));
end