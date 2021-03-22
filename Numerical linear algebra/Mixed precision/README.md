# Projects on Mixed precision

Mixed precision and numerical linear algebra is currently an interesting field to pursuit. Zhihan Li is now working in this field with many interesting problems.

## lowfrequency-roundoff.ipynb
In this project, we formulate the round-off of a float sequence $\alpha$ into a $0, 1$ sequence $\beta$ to a combinatorial optimization problems where we want to find an algorithm
which has $O(N)$ running time and returns $\beta$ s.t. $\alpha - \beta$ has first several Fourier coefficient of order $O(N^{-1})$.

Since this problem is a hard combinatorial optimization problem, another possible direction is to implement the cutting hyperplane algorithm for this combinatorial optimization.
