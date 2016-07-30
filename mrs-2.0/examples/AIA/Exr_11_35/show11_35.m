% txt file importing
% for output from Exercise 11.35

[ImX1 ImX2 ImY1 ImY2] = textread('AIAimage.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);

figure(2);
subplot(1,2,1)
fill([-2 2 2 -2]',[-2 -2 2 2]', 'b')
title('Initial subpaving', 'FontName', 'Ariel', 'FontSize', 12, 'FontWeight', 'bold')
subplot(1,2,2)
fill([ImX1'; ImX2'; ImX2'; ImX1'],[ImY1'; ImY1'; ImY2'; ImY2'], 'b');
title('Image subpaving using ImageSp', 'FontName', 'Ariel', 'FontSize', 12, 'FontWeight', 'bold')

