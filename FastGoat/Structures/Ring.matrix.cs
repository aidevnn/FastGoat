using FastGoat.Commons;
using FastGoat.Structures.VecSpace;

namespace FastGoat.Structures;

public static partial class Ring
{
    public static T[,] Matrix<T>(int rows, params T[] coefs)
    {
        if (coefs.Length % rows != 0)
            throw new ArgumentException();

        var cols = coefs.Length / rows;
        var A = new T[rows, cols];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                A[i, j] = coefs[i * cols + j];
            }
        }

        return A;
    }

    public static T[,] Matrix<T>(int rows, T t0, params int[] coefs) where T : IElt<T>, IRingElt<T>, IFieldElt<T>
    {
        return Matrix(rows, coefs.Select(c => t0.One.Mul(c)).ToArray());
    }

    public static T[,] Diagonal<T>(T v, int n) where T : IElt<T>, IRingElt<T>
    {
        var coefs = n.Range().Grid2D().Select(e => e.t1 == e.t2 ? v : v.Zero).ToArray();
        return Matrix(n, coefs);
    }

    public static T[,] Dot<T>(T[,] A, T[,] B, T t0) where T : IElt<T>, IRingElt<T>
    {
        var rowsA = A.GetLength(0);
        var colsAB = A.GetLength(1);
        if (colsAB != B.GetLength(0))
            throw new ArgumentException();

        var colsB = B.GetLength(1);
        var C = new T[rowsA, colsB];

        for (int i = 0; i < rowsA; i++)
        {
            for (int j = 0; j < colsAB; j++)
            {
                var sum = t0.Zero;
                for (int k = 0; k < colsAB; k++)
                    sum = sum.Add(A[i, k].Mul(B[k, j]));

                C[i, j] = sum;
            }
        }

        return C;
    }

    public static void Dot<T>(T[,] A, T[,] B, T[,] C, T t0) where T : IElt<T>, IRingElt<T>, IFieldElt<T>
    {
        if (A.GetLength(1) != B.GetLength(0) || A.GetLength(0) != C.GetLength(0) || B.GetLength(1) != C.GetLength(1))
            throw new ArgumentException();

        var rowsA = C.GetLength(0);
        var colsAB = A.GetLength(1);
        for (int i = 0; i < rowsA; i++)
        {
            for (int j = 0; j < colsAB; j++)
            {
                var sum = t0.Zero;
                for (int k = 0; k < colsAB; k++)
                    sum = sum.Add(A[i, k].Mul(B[k, j]));

                C[i, j] = sum;
            }
        }
    }

    public static T[,] Cofactor<T>(T[,] A, int row, int col) where T : IElt<T>, IRingElt<T>
    {
        var n = A.GetLength(0);
        if (n != A.GetLength(1))
            throw new ArgumentException();

        var B = new T[n - 1, n - 1];
        for (int i = 0; i < n; i++)
        {
            if (i == row) continue;
            int i0 = i < row ? i : i - 1;
            for (int j = 0; j < n; j++)
            {
                if (j == col) continue;
                int j0 = j < col ? j : j - 1;
                B[i0, j0] = A[i, j];
            }
        }

        return B;
    }

    public static T Determinant<T>(T[,] A, T t0) where T : IElt<T>, IRingElt<T>
    {
        var n = A.GetLength(0);
        if (n != A.GetLength(1))
            throw new ArgumentException();

        if (n == 1)
            return A[0, 0];

        var det = t0.Zero;
        var s = 1;
        for (int i = 0; i < n; i++)
        {
            var cof = Cofactor(A, i, 0);
            var det0 = Determinant(cof, t0);
            det = det.Add(det0.Mul(A[i, 0]).Mul(s));
            s *= -1;
        }

        return det;
    }

    public static T[,] Transpose<T>(T[,] A) where T : IElt<T>, IRingElt<T>
    {
        var rows = A.GetLength(0);
        var cols = A.GetLength(1);
        var B = new T[cols, rows];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                B[j, i] = A[i, j];
            }
        }

        return B;
    }

    public static T[,] CoMatrix<T>(T[,] A, T t0) where T : IElt<T>, IRingElt<T>
    {
        var n = A.GetLength(0);
        if (n != A.GetLength(1))
            throw new ArgumentException();

        var B = new T[n, n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                var s = (int)Math.Pow(-1, i + j);
                B[i, j] = Determinant(Cofactor(A, i, j), t0).Mul(s);
            }
        }

        return B;
    }

    private static void SwapRows<T>(int i, int j, T[,] A) where T : IElt<T>, IRingElt<T>
    {
        var cols = A.GetLength(1);
        for (int k = 0; k < cols; k++)
            (A[i, k], A[j, k]) = (A[j, k].Opp(), A[i, k]);
    }

    private static void MulRows<K>(int i, K a, K[,] A) where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var cols = A.GetLength(1);
        for (int k = 0; k < cols; k++)
            A[i, k] *= a;
    }

    private static void CombineRows<K>(int i, int j, K a, K[,] A) where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var cols = A.GetLength(1);
        for (int k = 0; k < cols; k++)
            A[i, k] = A[i, k].Add(A[j, k].Mul(a));
    }

    private static T[,] CloneMatrix<T>(T[,] A)
    {
        var m = A.GetLength(0);
        var n = A.GetLength(1);
        var B = new T[m, n];
        for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            B[i, j] = A[i, j];

        return B;
    }

    public static (K[,] P, K[,] A0) ReducedRowsEchelonForm<K>(K[,] A, bool details = true) where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var m = A.GetLength(0);
        var n = A.GetLength(1);
        if (m == 0 || n == 0)
            throw new ArgumentException();

        var k0 = A[0, 0];

        var A0 = CloneMatrix(A);
        var P = Diagonal(k0.One, n);

        for (int k = 0; k < n; k++)
        {
            var pivotFound = true;
            if (A0[k, k].IsZero())
            {
                pivotFound = false;
                var i = 0;
                for (i = k + 1; i < n; i++)
                {
                    if (!A0[i, k].IsZero())
                    {
                        SwapRows(i, k, A0);
                        SwapRows(i, k, P);
                        pivotFound = true;
                        
                        Console.WriteLine($"Swap {i} {k}");
                        Console.WriteLine("Matrix A'");
                        DisplayMatrix(A0);
                        Console.WriteLine("Matrix P");
                        DisplayMatrix(P);
                        
                        break;
                    }
                }
            }
            
            if(!pivotFound)
                continue;
            
            if (details)
            {
                Console.WriteLine($"Pivot L{k}");
                var a = A0[k, k];
                MulRows(k, a.Inv(), A0);
                MulRows(k, a.Inv(), P);
                Console.WriteLine("Matrix A'");
                DisplayMatrix(A0);
                Console.WriteLine("Matrix P");
                DisplayMatrix(P);
            }

            for (int i = 0; i < n; i++)
            {
                if (i == k)
                    continue;
                
                var b = A0[i, k];
                if (b.IsZero())
                    continue;

                var a0 = b.Opp();
                CombineRows(i, k, a0, A0);
                CombineRows(i, k, a0, P);
                if (details)
                {
                    Console.WriteLine($"L{i} <- L{i} + {a0} L{k}");
                    Console.WriteLine("Matrix A'");
                    DisplayMatrix(A0);
                    Console.WriteLine("Matrix P");
                    DisplayMatrix(P);
                }
            }
        }

        if (details)
        {
            Console.WriteLine("Final Matrix A'");
            DisplayMatrix(A0);
            Console.WriteLine("Final Matrix P");
            DisplayMatrix(P);
            
        }

        return (P, A0);
    }

    public static K DeterminantByPivot<K>(K[,] A, bool details = false)
        where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var dim = A.GetLength(0);
        if (dim == 0)
            throw new ArgumentException();

        if (dim != A.GetLength(1))
            return A[0, 0].Zero;

        var A0 = new K[dim, dim];
        for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            A0[i, j] = A[i, j];

        for (int k = 0; k < dim; k++)
        {
            if (A0[k, k].IsZero())
            {
                var i = 0;
                for (i = k + 1; i < dim; i++)
                {
                    if (!A0[i, k].IsZero())
                        break;
                }

                if (i == dim)
                    return A0[k, k];

                SwapRows(i, k, A0);
                if (details)
                {
                    Console.WriteLine($"Swap {i} {k}");
                    DisplayMatrix(A0);
                }
            }

            var a = A0[k, k];
            if (details)
                Console.WriteLine($"Pivot L{k}");

            for (int i = k + 1; i < dim; i++)
            {
                var b = A0[i, k];
                if (b.IsZero())
                    continue;

                var a0 = b.Mul(a.Inv()).Opp();
                CombineRows(i, k, a0, A0);
                if (details)
                {
                    Console.WriteLine($"L{i} <- L{i} + {a0} L{k}");
                    DisplayMatrix(A0);
                }
            }
        }

        return dim.Range().Select(k => A0[k, k]).Aggregate((a, b) => a.Mul(b));
    }

    public static void DisplayMatrix<T>(T[,] mat, string sep = ", ")
    {
        var rows = mat.GetLength(0);
        var cols = mat.GetLength(1);
        var rgCols = cols.Range();
        var digits = rgCols.Select(j => rows.Range().Max(i => $"{mat[i, j]}".Length)).Max();
        var fmt = $"{{0,{digits}}}";
        for (int i = 0; i < rows; i++)
        {
            var i0 = i;
            Console.WriteLine("[{0}]", rgCols.Select(j => mat[i0, j]).Glue(sep, fmt));
        }

        Console.WriteLine();
    }

    public static Polynomial<K, T>[] PolyBase<K, T>(K zero, T t, int n)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var ft = new Polynomial<K, T>(t, zero.One);
        var seq = new Stack<Polynomial<K, T>>(new[] { ft.One });
        while (seq.Count <= n)
        {
            var f0 = seq.Peek().Mul(ft);
            seq.Push(f0);
        }

        return seq.ToArray();
    }

    public static (SortedList<Polynomial<K, T>, Polynomial<K, T>>, SortedList<int, Polynomial<K, T>>)
        Decompose<K, T>(Polynomial<K, T> f, T ft)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var n = f.DegreeOf(ft);
        var baseFt = PolyBase(f.KZero, ft, n);
        var rem = f;
        var coefs = new SortedList<Polynomial<K, T>, Polynomial<K, T>>();
        var bs = new SortedList<int, Polynomial<K, T>>();
        foreach (var b in baseFt)
        {
            bs[b.DegreeOf(ft)] = b;
            if (rem.IsZero())
            {
                coefs[b] = f.Zero;
                continue;
            }

            var qr = rem.Div(b);
            if (!qr.quo.IsZero())
                coefs[b] = qr.quo;
            else
                coefs[b] = f.Zero;

            rem = qr.rem;
        }

        return (coefs, bs);
    }

    public static Polynomial<K, T>[,] SylvesterMatrix<K, T>(Polynomial<K, T> f, T ft, Polynomial<K, T> g, T gt)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        if (!f.Zero.Equals(g.Zero))
            throw new GroupException(GroupExceptionType.GroupDef);

        var (decF, bsF) = Decompose(f, ft);
        var (decG, bsG) = Decompose(g, gt);
        var m = decF.First().Key.DegreeOf(ft);
        var n = decG.First().Key.DegreeOf(gt);
        var S = new Polynomial<K, T>[m + n, m + n];
        for (int i = 0; i < m + n; i++)
        {
            for (int j = 0; j < m + n; j++)
            {
                if (i < n)
                {
                    if (j < i || j > i + m)
                        S[i, j] = f.Zero;
                    else
                    {
                        var mnm = bsF[m - j + i];
                        S[i, j] = decF[mnm];
                    }
                }
                else
                {
                    var i0 = i - n;
                    if (j < i0 || j > i0 + n)
                        S[i, j] = f.Zero;
                    else
                    {
                        var mnm = bsG[n - j + i0];
                        S[i, j] = decG[mnm];
                    }
                }
            }
        }

        return S;
    }

    public static Polynomial<K, T> Discriminant<K, T>(Polynomial<K, T> f, T X)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var g = f.D(X);
        if (g.IsZero())
            return f.Zero;

        var S = SylvesterMatrix(f, X, g, X);
        var am = f.CoefMax(X);
        var n = f.DegreeOf(X);
        var n0 = (int)(n * (n - 1) / 2);
        var s = n0 % 2 == 0 ? 1 : -1;
        var det = Determinant(S, f.Zero);
        return det.Div(am).quo.Mul(s);
    }

    public static K[,] SylvesterMatrix<K>(KPoly<K> f, KPoly<K> g)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        if (!f.Zero.Equals(g.Zero))
            throw new GroupException(GroupExceptionType.GroupDef);

        var m = f.Degree;
        var n = g.Degree;
        var S = new K[m + n, m + n];
        for (int i = 0; i < m + n; i++)
        {
            for (int j = 0; j < m + n; j++)
            {
                if (i < n)
                {
                    if (j < i || j > i + m)
                        S[i, j] = f.KZero;
                    else
                    {
                        S[i, j] = f[m - j + i];
                    }
                }
                else
                {
                    var i0 = i - n;
                    if (j < i0 || j > i0 + n)
                        S[i, j] = f.KZero;
                    else
                        S[i, j] = g[n - j + i0];
                }
            }
        }

        return S;
    }

    public static K Discriminant<K>(KPoly<K> f) where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var g = f.Derivative;
        if (g.IsZero())
            return f.KZero;

        var S = SylvesterMatrix(f, g);
        var am = f.Coefs.Last();
        var n = f.Degree;
        var n0 = (int)(n * (n - 1) / 2);
        var s = n0 % 2 == 0 ? 1 : -1;
        var det = DeterminantByPivot(S);
        return det.Mul(s).Mul(am.Inv());
    }
}