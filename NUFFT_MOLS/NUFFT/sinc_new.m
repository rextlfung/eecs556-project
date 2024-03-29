function out = sinc_new(x)

indicator = (abs(x)<1e-6);
x = x+indicator;

out = sin(x)./x;
out = out.*not(indicator)+indicator;
