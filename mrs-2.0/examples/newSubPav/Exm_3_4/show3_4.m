% txt file importing
% for output from Example 3.4, using new SPnode

[aX1 aX2 aY1 aY2] = textread('New3_4a.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);
[bX1 bX2 bY1 bY2] = textread('New3_4b.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);
[cX1 cX2 cY1 cY2] = textread('New3_4c.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);
[evalX1 evalX2 evalY1 evalY2] = textread('eval.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);

figure(4);
subplot(2,2,1);
fill([aX1'; aX2'; aX2'; aX1'],[aY1'; aY1'; aY2'; aY2'], 'b');
title('Initial subpaving', 'FontName', 'Ariel', 'FontSize', 12, 'FontWeight', 'bold')

subplot(2,2,2);
fill([evalX1'; evalX2'; evalX2'; evalX1'],[evalY1'; evalY1'; evalY2'; evalY2'], 'b');
title('After Evaluate(): subpaving images', 'FontName', 'Ariel', 'FontSize', 12, 'FontWeight', 'bold')

subplot(2,2,4);
fill([bX1'; bX2'; bX2'; bX1'],[bY1'; bY1'; bY2'; bY2'], 'b');
title('Image subpaving', 'FontName', 'Ariel', 'FontSize', 12, 'FontWeight', 'bold')

subplot(2,2,3);
fill([cX1'; cX2'; cX2'; cX1'],[cY1'; cY1'; cY2'; cY2'], 'b');
title('Subpaving for reciprocal image of the image subpaving', 'FontName', 'Ariel', 'FontSize', 12, 'FontWeight', 'bold')


