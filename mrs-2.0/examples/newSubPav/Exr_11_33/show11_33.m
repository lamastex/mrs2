% txt file importing
% for output from Exercise 11.33, using new SPnodes

[AnnX1 AnnX2 AnnY1 AnnY2] = textread('new_annular.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);
[DirX1 DirX2 DirY1 DirY2] = textread('new_direct.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);
[InvX1 InvX2 InvY1 InvY2] = textread('new_inverse.txt', '%*s %f %*s %f %*s %*s %*s %f %*s %f %*s', 'headerlines', 3);

figure(1);
subplot(1,3,1);
fill([AnnX1'; AnnX2'; AnnX2'; AnnX1'],[AnnY1'; AnnY1'; AnnY2'; AnnY2'], 'b');
subplot(1,3,2);
fill([DirX1'; DirX2'; DirX2'; DirX1'],[DirY1'; DirY1'; DirY2'; DirY2'], 'b');
subplot(1,3,3);
fill([InvX1'; InvX2'; InvX2'; InvX1'],[InvY1'; InvY1'; InvY2'; InvY2'], 'b');
