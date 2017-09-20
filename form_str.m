function sss = form_str(vec,a)
if ~isempty(vec)
    s = [', ',num2str(vec(1))];
    vec = vec(2:end);
    sss = [s,form_str(vec,a)];
else
    sss = a;
end

