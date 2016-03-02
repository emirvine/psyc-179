function x_avg = average_x_by_y_bin(x, y, y_edges)
[~, idx] = histc(y, y_edges);
x_avg = zeros(size(y_edges));
for bin = length(y_edges):-1:1
    if sum(idx == bin) ~= 0
        x_avg(bin) = nanmean(x(idx == bin));
    end
end
end