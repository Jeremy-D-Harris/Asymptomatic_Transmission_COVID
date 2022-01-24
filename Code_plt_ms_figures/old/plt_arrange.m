% figs = [f1, f2, f3];   %as many as needed
nfig = 4; %length(figs);
frac = 1/nfig;
for K = 1 : nfig
%   old_pos = get(figs(K), 'Position');
  set(figs(K), 'Position', [(K-1)*frac, old_pos(2), frac, old_pos(4)]);
end