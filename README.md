# Gram–Schmidt Orthogonalization

A short write-up and a MATLAB implementation of the Gram–Schmidt process, which turns a set of linearly independent vectors into an orthonormal basis of the space they span.

## Contents

| File | Description |
|---|---|
| `gschmit.m` | MATLAB function `gschmit(y)`; columns of `y` are the input vectors |
| `example.m` | Example script |
| `report.pdf` | Original notes |

## The idea

Given linearly independent vectors $y_1, \dots, y_n$, build orthogonal vectors $v_i$ by removing from each $y_i$ its components along the earlier $v_k$:

$$
v_i = y_i - \sum_{k=1}^{i-1} \frac{\langle v_k, y_i\rangle}{\langle v_k, v_k\rangle}\, v_k
$$

Then normalize:

$$
\hat v_i = \frac{v_i}{\lVert v_i \rVert}, \qquad i = 1, \dots, n
$$

### Step-by-step (modified form)

1. $v_1 = y_1$ and $\hat v_1 = v_1 / \lVert v_1\rVert$
2. $y_j' = y_j - \langle \hat v_1, y_j\rangle \hat v_1$ for $j = 2,\dots,n$
3. $v_2 = y_2'$ and $\hat v_2 = v_2 / \lVert v_2\rVert$
4. $y_j'' = y_j' - \langle \hat v_2, y_j'\rangle \hat v_2$ for $j = 3,\dots,n$
5. Repeat until all vectors are processed.

## Worked example

$$
y_1=\begin{bmatrix}1\\0\\1\end{bmatrix},\quad
y_2=\begin{bmatrix}-1\\2\\1\end{bmatrix},\quad
y_3=\begin{bmatrix}0\\1\\2\end{bmatrix}
$$

**Vector 1.** $\hat v_1 = \tfrac{1}{\sqrt2}(1,0,1)^T$.

**Vector 2.** $\langle \hat v_1, y_2\rangle = 0$, so $y_2' = y_2$ and $\hat v_2 = \tfrac{1}{\sqrt6}(-1,2,1)^T$.

**Vector 3.**
- $\langle \hat v_1, y_3\rangle = \sqrt2$, so $y_3' = (-1,1,1)^T$
- $\langle \hat v_2, y_3'\rangle = 4/\sqrt6$, so $y_3'' = \left(-\tfrac13,-\tfrac13,\tfrac13\right)^T$
- $\hat v_3 = \tfrac{1}{\sqrt3}(-1,-1,1)^T$

**Result:** $\hat v_1,\hat v_2,\hat v_3$ are orthonormal, i.e. $Q^TQ = I$.

## MATLAB usage

```matlab
y = [1 -1 0; 0 2 1; 1 1 2];   % columns are y1, y2, y3
v = gschmit(y);
disp(v' * v)                  % should be the identity
```

## Known issue in `gschmit.m`

The inner loop overwrites `v(:,i)` on every pass of `k`, so only the projection on the *last* earlier vector is removed. For the example above, the third output vector is not orthogonal to the first. The fix is to accumulate the projections:

```matlab
function v = gschmit(y)
[~, M] = size(y);
v = y;
for i = 2:M
    for k = 1:i-1
        v(:,i) = v(:,i) - (v(:,k)' * y(:,i)) / (v(:,k)' * v(:,k)) * v(:,k);
    end
end
for i = 1:M
    v(:,i) = v(:,i) / norm(v(:,i));
end
end
```

Using `v(:,i)` instead of `y(:,i)` in the coefficient gives the numerically more stable modified Gram–Schmidt.

## License

MIT

[Hossein Mosaffa's personal website](https://hossein-mosaffa.github.io)
