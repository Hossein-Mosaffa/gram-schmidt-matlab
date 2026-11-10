function Q = gramschmidt(Y, tol)
%GRAMSCHMIDT  Orthonormalize the columns of a matrix (modified Gram-Schmidt).
%
%   Q = GRAMSCHMIDT(Y) returns a matrix Q whose columns form an orthonormal
%   basis for the span of the columns of Y. Columns of Y must be linearly
%   independent.
%
%   Q = GRAMSCHMIDT(Y, TOL) errors out if a column's norm falls below TOL
%   after projection, which signals linear dependence. Default: 1e-12.
%
%   Example:
%       Y = [1 -1 0; 0 2 1; 1 1 2];
%       Q = gramschmidt(Y);
%       disp(Q' * Q)        % identity matrix
%
%   See also QR, ORTH.

    if nargin < 2
        tol = 1e-12;
    end

    Q = Y;
    numVectors = size(Y, 2);

    for i = 1:numVectors
        % Remove the components along all previously computed directions.
        for k = 1:i-1
            Q(:, i) = Q(:, i) - (Q(:, k)' * Q(:, i)) * Q(:, k);
        end

        len = norm(Q(:, i));
        if len < tol
            error('gramschmidt:dependent', ...
                  'Column %d is linearly dependent on earlier columns.', i);
        end

        Q(:, i) = Q(:, i) / len;
    end
end
