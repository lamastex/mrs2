% txt file importing
% for output from Exercise 11.35 using new SPnodes

[ImX1 ImX2 ImY1 ImY2] = textread('new_image.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);

figure(2);
fill([ImX1'; ImX2'; ImX2'; ImX1'],[ImY1'; ImY1'; ImY2'; ImY2'], 'b');

