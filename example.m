% EXAMPLE  Gram-Schmidt on three vectors in R^3.

clear; clc;

Y = [ 1  -1   0
      0   2   1
      1   1   2 ];          % columns are y1, y2, y3

Q = gramschmidt(Y);

disp('Orthonormal basis (columns):');
disp(Q);

disp('Check: Q''*Q (should be the identity):');
disp(Q' * Q);

fprintf('Max deviation from identity: %.2e\n', max(abs(Q'*Q - eye(3)), [], 'all'));
