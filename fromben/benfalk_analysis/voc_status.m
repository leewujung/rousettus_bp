%for mark_dur.m and mark_good_dur.m
function voc_status(hh,offset,txt,col,pause_amnt)
axes(hh(1))
aa=axis;
text(aa(1)+offset,aa(3),txt,'fontsize',128,...
  'verticalalign','bottom','color',col)
axes(hh(2))
text(aa(1)+offset,aa(3),txt,'fontsize',128,...
  'verticalalign','bottom','color',col)
pause(pause_amnt)