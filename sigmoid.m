function y = sigmoid(x,a) %a is GAIN
   y = 1./(1+exp(a*(0.5-x)));
end