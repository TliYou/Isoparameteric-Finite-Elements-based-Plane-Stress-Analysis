function y = Q6Assemble(K, k, node)

lp = [];
for c=1:length(node)
    lp = [lp 2*node(c)-1 2*node(c)];
end

    
for r = 1:length(lp)
    for s = 1:length(lp)
        K(lp(r),lp(s)) = k(r,s);
    end
end

y = K;

end