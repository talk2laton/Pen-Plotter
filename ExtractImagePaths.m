function [Paths, xmin, xmax, ymin, ymax] = ExtractImagePaths(filename, stp)
img = imread(filename);
img2 = rgb2gray(img);
img3 = ~imbinarize(img2, 'adaptive', 'ForegroundPolarity', 'dark' );
img4 = bwmorph(img3, 'clean');
img5 = bwmorph(img4, 'thin', inf);
coordsPix = getCoords(img5);
segmentsPix = coords2segments(coordsPix);
Paths = connectSegments(segmentsPix);
[Paths, xmin, xmax, ymin, ymax] = ExtractPaths(Paths, stp);