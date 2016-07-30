% txt file importing
% for output from Example 3.3

[aX1 aX2 aY1 aY2] = textread('AIA3_3a.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);
[evalX1 evalX2 evalY1 evalY2] = textread('eval.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);
[dX1 dX2 dY1 dY2] = textread('AIA3_3d.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);


figure(3);
subplot(1,3,1);
fill([aX1'; aX2'; aX2'; aX1'],[aY1'; aY1'; aY2'; aY2'], 'b');
title('Initial subpaving', 'FontName', 'Ariel', 'FontSize', 12, 'FontWeight', 'bold')
subplot(1,3,2);
fill([evalX1'; evalX2'; evalX2'; evalX1'],[evalY1'; evalY1'; evalY2'; evalY2'], 'b');
title('After Evaluate(): Set of images of initial subpaving', 'FontName', 'Ariel', 'FontSize', 12, 'FontWeight', 'bold')
subplot(1,3,3);
fill([dX1'; dX2'; dX2'; dX1'],[dY1'; dY1'; dY2'; dY2'], 'b');
title('After Regularize(): A mininal subpaving covering the image', 'FontName', 'Ariel', 'FontSize', 12, 'FontWeight', 'bold')

