function f = dist_minimizer(a,x)
Dx=nan(size(x,1),1);
for k=1:size(x,1)
  DD=(x(k,1) + a - x(:,2)).^2;
  Dx(k)=min(DD);
end
f=nansum(Dx);