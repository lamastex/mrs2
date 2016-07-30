% txt file importing
% for output from Example 3.3, new version for our SPnode class

[aX1 aX2 aY1 aY2] = textread('New3_3a.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);
[evalX1 evalX2 evalY1 evalY2] = textread('eval.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);
[dX1 dX2 dY1 dY2] = textread('New3_3d.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);


figure(3);
subplot(1,3,1);
fill([aX1'; aX2'; aX2'; aX1'],[aY1'; aY1'; aY2'; aY2'], 'b');
subplot(1,3,2);
fill([evalX1'; evalX2'; evalX2'; evalX1'],[evalY1'; evalY1'; evalY2'; evalY2'], 'b');
subplot(1,3,3);
fill([dX1'; dX2'; dX2'; dX1'],[dY1'; dY1'; dY2'; dY2'], 'b');

