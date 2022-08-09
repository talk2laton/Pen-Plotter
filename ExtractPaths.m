function [Paths, xmin, xmax, ymin, ymax] = ExtractPaths(Paths, stp)
minxy = min(fliplr(Paths{1}));
maxxy = max(fliplr(Paths{1}));
for i = 1:numel(Paths)
    Paths{i} = fliplr(Paths{i});
    minxy = min(minxy, min(Paths{i}));
    maxxy = max(maxxy, max(Paths{i}));
end
xmin = minxy(1); xmax = maxxy(1);
ymin = minxy(2); ymax = maxxy(2);
for i = 1:numel(Paths)
    Paths{i}(:,1) = xmax - Paths{i}(:,1);
end
xmax = xmax - xmin; xmin = 0;
for i = 1:numel(Paths)
    path = Paths{i};
    s = size(path, 1);
    Paths{i} = path(1:stp:end, :);
end
