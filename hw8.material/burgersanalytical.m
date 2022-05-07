function [x,u]=burgersanalytical(x_,fh,t)
u=feval(fh,x_);
x=x_+u*t;


