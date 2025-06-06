using System.Text;
using FastGoat.Commons;
using FastGoat.Structures.VecSpace;

namespace FastGoat.Structures;

public static partial class Ring
{
    public static T[,] Matrix<T>(int rows, T[] coefs)
    {
        if (coefs.Length % rows != 0)
            throw new ArgumentException($"nb coefs = {coefs.Length} and rows = {rows}");

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

    public static T[,] Matrix<T>(int rows, T t0, params dynamic[] coefs) where T : IElt<T>, IRingElt<T>, IFieldElt<T>
    {
        return Matrix(rows, coefs.Select(c => c is int c0 ? t0.One.Mul(c0) : c is T c1 ? c1 : t0.Zero).Prepend(t0).ToArray());
    }

    public static T[,] Diagonal<T>(T v, int n) where T : IElt<T>, IRingElt<T>
    {
        var coefs = n.Range().Grid2D().Select(e => e.t1 == e.t2 ? v : v.Zero).ToArray();
        return Matrix(n, coefs);
    }

    public static T[,] Matrix<T>(T v, int nbRows, int nbCols)
    {
        return Matrix(nbRows, Enumerable.Repeat(v, nbRows * nbCols).ToArray());
    }

    public static T[,] Dot<T>(T[,] A, T[,] B) where T : IElt<T>, IRingElt<T>
    {
        var rowsA = A.GetLength(0);
        var colsAB = A.GetLength(1);
        if (colsAB != B.GetLength(0) || rowsA * colsAB == 0)
            throw new ArgumentException();

        var colsB = B.GetLength(1);
        var t0 = A[0, 0];
        var C = new T[rowsA, colsB];

        for (int i = 0; i < rowsA; i++)
        {
            for (int j = 0; j < colsB; j++)
            {
                var sum = t0.Zero;
                for (int k = 0; k < colsAB; k++)
                    sum = sum.Add(A[i, k].Mul(B[k, j]));

                C[i, j] = sum;
            }
        }

        return C;
    }

    public static void Dot<T>(T[,] A, T[,] B, T[,] C) where T : IElt<T>, IRingElt<T>
    {
        if (A.GetLength(1) != B.GetLength(0) || A.GetLength(0) != C.GetLength(0) || B.GetLength(1) != C.GetLength(1))
            throw new ArgumentException();

        var rowsA = A.GetLength(0);
        var colsAB = B.GetLength(1);
        var t0 = A[0, 0];
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

    public static void SwapCols<T>(int i, int j, T[,] A) where T : IElt<T>, IRingElt<T>
    {
        var rows = A.GetLength(0);
        for (int k = 0; k < rows; k++)
            (A[k, i], A[k, j]) = (A[k, j].Opp(), A[k, i]);
    }

    public static void MulCols<K>(int i, K a, K[,] A) where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var rows = A.GetLength(0);
        for (int k = 0; k < rows; k++)
            A[k, i] *= a;
    }

    public static void CombineCols<K>(int i, int j, K a, K[,] A) where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var rows = A.GetLength(0);
        for (int k = 0; k < rows; k++)
            A[k, i] = A[k, i].Add(A[k, j].Mul(a));
    }

    public static void SwapRows<T>(int i, int j, T[,] A) where T : IElt<T>, IRingElt<T>
    {
        var cols = A.GetLength(1);
        for (int k = 0; k < cols; k++)
            (A[i, k], A[j, k]) = (A[j, k].Opp(), A[i, k]);
    }

    public static void MulRows<K>(int i, K a, K[,] A) where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var cols = A.GetLength(1);
        for (int k = 0; k < cols; k++)
            A[i, k] *= a;
    }

    public static void CombineRows<K>(int i, int j, K a, K[,] A) where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var cols = A.GetLength(1);
        for (int k = 0; k < cols; k++)
            A[i, k] = A[i, k].Add(A[j, k].Mul(a));
    }

    public static T[,] CloneMatrix<T>(T[,] A)
    {
        var m = A.GetLength(0);
        var n = A.GetLength(1);
        var B = new T[m, n];
        for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            B[i, j] = A[i, j];

        return B;
    }

    public static KMatrix<K> ToKMatrix<K>(this K[,] mat)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(mat);
    }

    public static KMatrix<K> ToKMatrix<K>(this IEnumerable<K> mat)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(Matrix(1, mat.ToArray()));
    }

    public static KMatrix<K> ToKMatrix<K>(this IEnumerable<K> mat, int m)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(Matrix(m, mat.ToArray()));
    }

    public static KMatrix<K> ToKMatrix<K>(this IEnumerable<int> mat, int m, K scalar)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(Matrix(m, mat.Select(i => i * scalar.One).ToArray()));
    }

    public static K SquareNorm2<K>(KMatrix<K> mat)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (mat.M == 1)
        {
            var sq = mat.KZero;
            for (int i = 0; i < mat.N; i++)
                sq += mat.Coefs[0, i] * mat.Coefs[0, i];

            return sq;
        }

        if (mat.N == 1)
        {
            var sq = mat.KZero;
            for (int i = 0; i < mat.M; i++)
                sq += mat.Coefs[i, 0] * mat.Coefs[i, 0];

            return sq;
        }

        throw new("SquareNorm2");
    }

    public static K SquareNorm2f<K>(KMatrix<K> mat) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
    {
        if (mat.M == 1)
        {
            var sq = mat.KZero;
            for (int i = 0; i < mat.N; i++)
                sq += mat.Coefs[0, i].Absolute2;

            return sq;
        }

        if (mat.N == 1)
        {
            var sq = mat.KZero;
            for (int i = 0; i < mat.M; i++)
                sq += mat.Coefs[i, 0].Absolute2;

            return sq;
        }

        throw new("SquareNorm2");
    }

    public static K Norm2f<K>(KMatrix<K> mat) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
    {
        return K.Sqrt(SquareNorm2f(mat));
    }

    public static (K[,] P, K[,] A0) ReducedRowsEchelonForm<K>(K[,] A)
        where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (m, n) = (A.GetLength(0), A.GetLength(1));
        if (m == 0 || n == 0)
            throw new ArgumentException();

        var A0 = CloneMatrix(A);
        var P = Diagonal(A[0, 0].One, m);

        int i = 0, j = 0;
        while (i < m && j < n)
        {
            var pivotFound = true;
            if (A0[i, j].IsZero())
            {
                pivotFound = false;
                for (int i0 = i + 1; i0 < m; i0++)
                {
                    if (!A0[i0, j].IsZero())
                    {
                        SwapRows(i0, i, A0);
                        SwapRows(i0, i, P);
                        pivotFound = true;
                        break;
                    }
                }
            }

            if (!pivotFound)
            {
                ++j;
                continue;
            }

            var a = A0[i, j];
            MulRows(i, a.Inv(), A0);
            MulRows(i, a.Inv(), P);
            for (int i0 = 0; i0 < m; i0++)
            {
                if (i0 == i)
                    continue;

                var b = A0[i0, j];
                if (b.IsZero())
                    continue;

                var a0 = b.Opp();
                CombineRows(i0, i, a0, A0);
                CombineRows(i0, i, a0, P);
            }

            ++i;
            ++j;
        }

        return (P, A0);
    }

    // wikipedia
    public static K[,] ToReducedRowEchelonForm<K>(K[,] M0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var M = CloneMatrix(M0);
        var lead = 0;
        var rowCount = M.GetLength(0);
        var columnCount = M.GetLength(1);
        for (int r = 0; r < rowCount; r++)
        {
            if (columnCount <= lead)
                return M;

            var i = r;
            while (M[i, lead].IsZero())
            {
                ++i;
                if (rowCount == i)
                {
                    i = r;
                    ++lead;
                    if (columnCount == lead)
                        return M;
                }
            }

            if (i != r) SwapRows(i, r, M);
            MulRows(r, M[r, lead].Inv(), M);
            for (int j = 0; j < rowCount; j++)
            {
                if (j != r)
                    CombineRows(j, r, M[j, lead].Opp(), M);
            }

            ++lead;
        }

        return M;
    }


    public static (KMatrix<K> P, KMatrix<K> A0) ReducedRowsEchelonForm<K>(KMatrix<K> A)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var e = ReducedRowsEchelonForm(A.Coefs);
        return (new(e.P), new(e.A0));
    }

    public static (K[,] P, K[,] R) RowsEchelonForm<K>(K[,] A, bool details = false)
        where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (m, n) = (A.GetLength(0), A.GetLength(1));
        if (m == 0 || n == 0 || m > n)
            throw new ArgumentException();


        // Console.WriteLine(new { dim });
        var A0 = CloneMatrix(A);
        var P = Diagonal(A[0, 0].One, m);

        for (int k = 0; k < m; k++)
        {
            // Console.WriteLine(new { k, dim });
            if (A0[k, k].IsZero())
            {
                var i = 0;
                for (i = k + 1; i < m; i++)
                {
                    if (!A0[i, k].IsZero())
                        break;
                }

                if (i == m)
                    continue;

                SwapRows(i, k, A0);
                SwapRows(i, k, P);
                if (details)
                {
                    Console.WriteLine($"Swap {i} {k}");
                    DisplayMatrix(A0);
                }
            }

            var a = A0[k, k];
            if (details)
                Console.WriteLine($"Pivot L{k}");

            for (int i = k + 1; i < m; i++)
            {
                var b = A0[i, k];
                if (b.IsZero())
                    continue;

                var a0 = b.Mul(a.Inv()).Opp();
                CombineRows(i, k, a0, A0);
                CombineRows(i, k, a0, P);
                if (details)
                {
                    Console.WriteLine($"L{i} <- L{i} + {a0} L{k}");
                    DisplayMatrix(A0);
                }
            }
        }

        // DisplayMatrix(A0);
        return (P, A0);
    }

    public static (KMatrix<K> P, KMatrix<K> R) RowsEchelonForm<K>(KMatrix<K> mat)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (P, R) = RowsEchelonForm(mat.Coefs);
        return (new(P), new(R));
    }

    public static K DeterminantByPivot<K>(K[,] A, bool details = false)
        where K : IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var dim = A.GetLength(0);
        if (dim == 0)
            throw new ArgumentException();

        if (dim != A.GetLength(1))
            return A[0, 0].Zero;

        // Console.WriteLine(new { dim });
        var A0 = new K[dim, dim];
        for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            A0[i, j] = A[i, j];

        for (int k = 0; k < dim; k++)
        {
            // Console.WriteLine(new { k, dim });
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

        // DisplayMatrix(A0);
        return dim.Range().Select(k => A0[k, k]).Aggregate((a, b) => a.Mul(b));
    }

    public static string Matrix2String<T>(T[,] mat, string sep = ", ")
    {
        var cursorPos = Enumerable.Repeat(' ', Console.CursorLeft).Glue();
        var rows = mat.GetLength(0);
        var cols = mat.GetLength(1);
        var rgCols = cols.Range();
        var rgRows = rows.Range();
        if (MatrixDisplayForm == MatrixDisplay.OneLineArray)
            return $"[{rgRows.Grid2D(rgCols).Select(e => mat[e.t1, e.t2]).Glue(sep)}]";
        
        var digitsByCols = rgCols.Select(j => rgRows.Max(i => $"{mat[i, j]}".Length))
            .Select(d => MatrixDisplayForm == MatrixDisplay.TableLeft ? $"{{0,-{d}}}" : $"{{0,{d}}}").ToArray();
        if (MatrixDisplayForm == MatrixDisplay.CurlyBracketNoFmt || MatrixDisplayForm == MatrixDisplay.SquareBracketNoFmt)
            digitsByCols = digitsByCols.Select(d => "{0}").ToArray();

        var lt = new List<string>();
        for (int i = 0; i < rows; i++)
        {
            var i0 = i;
            var s0 = rgCols.Select(j => string.Format(digitsByCols[j], mat[i0, j])).Glue(sep);
            if (MatrixDisplayForm == MatrixDisplay.CurlyBracket)
                lt.Add($"{{ {s0} }}");
            else
                lt.Add($"[{s0}]");
        }

        if (MatrixDisplayForm == MatrixDisplay.Table || MatrixDisplayForm == MatrixDisplay.TableLeft)
        {
            return lt.Select((l, i) => i == 0 ? l : cursorPos + l).Glue("\n");
        }
        else
        {
            if (MatrixDisplayForm == MatrixDisplay.CurlyBracket)
                return $"{{ {lt.Glue(",")} }}";
            else
                return $"[{lt.Glue(",")}]";
        }
    }

    public enum MatrixDisplay
    {
        CurlyBracket,
        SquareBracket,
        CurlyBracketNoFmt,
        SquareBracketNoFmt,
        OneLineArray,
        Table,
        TableLeft
    }

    public static MatrixDisplay MatrixDisplayForm { get; set; } = MatrixDisplay.Table;

    public static void DisplayMatrix<T>(T[,] mat, MatrixDisplay display, string sep = ", ")
    {
        var prev = MatrixDisplayForm;
        MatrixDisplayForm = display;
        Console.WriteLine(Matrix2String(mat, sep));
        MatrixDisplayForm = prev;
    }

    public static void DisplayMatrix<T>(T[,] mat, string sep = ", ")
    {
        Console.WriteLine(Matrix2String(mat, sep));
    }

    public static Polynomial<K, T>[] PolyBase<K, T>(Polynomial<K, T> f, T t, int n)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var ft = new Polynomial<K, T>(new Monom<T>(f.Indeterminates, t), f.KOne);
        var seq = new Stack<Polynomial<K, T>>(new[] { ft.One });
        while (seq.Count <= n)
        {
            var f0 = seq.Peek().Mul(ft);
            seq.Push(f0);
        }

        return seq.ToArray();
    }

    public static Polynomial<K, T>[] PolyBase<K, T>(Polynomial<K, T> f, T t)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return PolyBase(f, t, f.DegreeOf(t));
    }

    public static Polynomial<K, T>[] PolyBase<K, T>(Polynomial<K, T> f, string x)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var t = f.Indeterminates[x];
        return PolyBase(f, t, f.DegreeOf(t));
    }

    public static (SortedList<Polynomial<K, T>, Polynomial<K, T>>, SortedList<int, Polynomial<K, T>>)
        Decompose<K, T>(Polynomial<K, T> f, T ft)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var n = f.DegreeOf(ft);
        var baseFt = PolyBase(f, ft, n);
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

    public static (SortedList<Polynomial<K, T>, Polynomial<K, T>>, SortedList<int, Polynomial<K, T>>)
        Decompose<K, T>(Polynomial<K, T> f, string x)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var t = f.Indeterminates[x];
        return Decompose(f, t);
    }

    public static Polynomial<K, T>[,] SylvesterMatrix<K, T>(Polynomial<K, T> f, T ft, Polynomial<K, T> g, T gt)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        if (!f.Zero.Equals(g.Zero))
            throw new GroupException(GroupExceptionType.GroupDef);

        var (decF, bsF) = Decompose(f, ft);
        var (decG, bsG) = Decompose(g, gt);
        var m = decF.Last().Key.DegreeOf(ft);
        var n = decG.Last().Key.DegreeOf(gt);
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
        var n0 = n * (n - 1) / 2;
        var s = n0 % 2 == 0 ? 1 : -1;
        
        var nbInd = f.NbIndeterminates;
        if (nbInd == 2)
        {
            var Y = f.ExtractAllIndeterminates.First(xi => !xi.Equals(X));
            var (dim1, dim2) = (S.GetLength(0), S.GetLength(1));
            var S1 = new FracPoly<K>[dim1, dim2];
            for (int i = 0; i < dim1; i++)
            {
                for (int j = 0; j < dim2; j++)
                    S1[i, j] = new(S[i, j].ToKPoly(Y));
            }
        
            var det = DeterminantByPivot(S1);
            return det.Num.ToPolynomial(f.Indeterminates, Y).Div(am).quo.Mul(s);
        }
        else
        {
            var det = Determinant(S, f.Zero);
            return det.Div(am).quo.Mul(s);
        }
    }

    public static Polynomial<K, T> Discriminant<K, T>(Polynomial<K, T> f, Polynomial<K, T> X)
        where T : struct, IElt<T>
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Discriminant(f, X.ExtractIndeterminate);
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

    public static K Resultant<K>(KPoly<K> f, KPoly<K> g)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return DeterminantByPivot(SylvesterMatrix(f, g));
    }

    public static K Discriminant<K>(KPoly<K> f) where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var g = f.Derivative;
        if (g.IsZero())
            return f.KZero;

        var am = f.Coefs.Last();
        var n = f.Degree;
        var n0 = (int)(n * (n - 1) / 2);
        var s = n0 % 2 == 0 ? 1 : -1;
        var det = FastResultant(f, g);
        return det.Mul(s).Mul(am.Inv());
    }

    public static K FastResultant<K>(KPoly<K> A, KPoly<K> B)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (f, g, s) = (A.KOne, A.KOne, A.KOne);
        while (B.Degree > 0)
        {
            var d = A.Degree - B.Degree;
            var bDeg = B.Coefs.Last().Pow(d + 1);
            var R = (bDeg * A).Div(B).rem; // Pseudo remainders, AECF
            if ((A.Degree * B.Degree) % 2 == 1)
                s = -s;

            A = B;
            B = R / (f * g.Pow(d));
            f = A.Coefs.Last();
            g = f * (f / g).Pow(d - 1);
        }

        var dA = A.Degree;
        return ((s * B.Pow(dA)) / g.Pow(dA - 1))[0];
    }

    static (KPoly<K>, KPoly<K>) FastBezoutInternal<K>(KPoly<K> A, KPoly<K> B, K w0, K y0, int d0, int i = 0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (B.IsZero())
            return A.CompareTo(A.Opp()) == -1 ? (A.One.Opp(), A.One) : (A.One, A.Zero);

        var d1 = A.Degree - B.Degree;
        var y1 = B.LT;
        var w1 = i == 0 ? -A.KOne : (-y0).Pow(d0) * w0.Pow(d0 - 1).Inv();
        var b = i == 0 ? ((d0 + 1) % 2 == 0 ? A.KOne : -A.KOne) : -y0 * w1.Pow(d1);
        var bi = b.Inv();
        var cf = y1.Pow(d1 + 1);
        var (q, r) = (cf * A).Div(B);
        var (X0, Y0) = FastBezoutInternal(B, r * bi, w1, y1, d1, i + 1);
        return (Y0 * bi * cf, X0 - (q * Y0 * bi));
    }

    public static (KPoly<K>, KPoly<K>) FastBezout<K>(KPoly<K> A, KPoly<K> B)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (A.Degree < B.Degree)
        {
            var (y, x) = FastBezout(B, A);
            return (x, y);
        }

        return FastBezoutInternal(A, B, A.KOne, A.KOne, 0);
    }

    public static KPoly<K> FastInv<K>(KPoly<K> A, KPoly<K> mod)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (U, V) = FastBezout(A, mod);
        var gcd = (A * U + mod * V);
        return (U / gcd).Div(mod).rem;
    }

    // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Pseudo-remainder_sequences
    static KPoly<K> FastGCDInternal<K>(KPoly<K> A, KPoly<K> B, K w0, K y0, int d0, int i = 0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (B.IsZero())
            return A;

        var d1 = A.Degree - B.Degree;
        var y1 = B.LT;
        var w1 = i == 0 ? -A.KOne : FastPow(-y0, d0) / FastPow(w0, d0 - 1);
        var b = i == 0 ? ((d0 + 1) % 2 == 0 ? A.KOne : -A.KOne) : -y0 * FastPow(w1, d1);
        var r = (FastPow(y1, d1 + 1) * A).Div(B).rem / b;
        return FastGCDInternal(B, r, w1, y1, d1, i + 1);
    }

    public static KPoly<K> FastGCD<K>(KPoly<K> A, KPoly<K> B) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (A.Degree < B.Degree)
            return FastGCD(B, A);

        return FastGCDInternal(A, B, A.KOne, A.KOne, 0);
    }

    public static KPoly<K> FastGCD<K>(KPoly<K>[] arr) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (arr.Length == 0)
            throw new ArgumentException();

        if (arr.Length == 1)
            return arr[0];

        return FastGCD(arr.First(), FastGCD(arr.Skip(1).ToArray()));
    }
    
    public static KPoly<K> FastLCM<K>(KPoly<K> A, KPoly<K> B) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return (A * B) / FastGCD(A, B);
    }
    
    public static KPoly<K> FastLCM<K>(KPoly<K>[] arr) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (arr.Length == 0)
            throw new ArgumentException();

        if (arr.Length == 1)
            return arr[0];

        return FastLCM(arr.First(), FastLCM(arr.Skip(1).ToArray()));
    }

    public static (KMatrix<K> O, KMatrix<K> U) GramSchmidt<K>(KMatrix<K> A)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var n = A.M;
        var vs = new KMatrix<K>[n];
        var u = Diagonal(A.KOne, n);
        for (int i = 0; i < n; i++)
        {
            vs[i] = A.GetCol(i);
            for (int j = 0; j < i; j++)
            {
                var uij = u[i, j] = ((A.GetCol(i).T * vs[j]) / (vs[j].T * vs[j]))[0, 0];
                vs[i] -= uij * vs[j];
            }
        }

        return (KMatrix<K>.MergeSameRows(vs), new(u));
    }

    public static (KMatrix<K> O, KMatrix<K> U) GramSchmidt2<K>(KMatrix<K> A, bool details = false)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (m, n) = A.Dim;
        var vs = new KMatrix<K>[n];
        var u = Diagonal(A.KOne, n);
        var allCols = A.Cols;
        for (int i = 0; i < n; i++)
        {
            vs[i] = allCols[i];
            for (int j = 0; j < i; j++)
            {
                var sum0 = A.KZero;
                var sum1 = A.KZero;
                for (int k = 0; k < m; k++)
                {
                    sum0 += vs[i].Coefs[k, 0] * vs[j].Coefs[k, 0];
                    sum1 += vs[j].Coefs[k, 0] * vs[j].Coefs[k, 0];
                }

                // var uij = u[i, j] = (ScalarProduct(vs[i], vs[j]) / ScalarProduct(vs[j], vs[j]));
                var uij = u[i, j] = sum0 / sum1;

                // vs[i] -= uij * vs[j];
                for (int k = 0; k < m; k++)
                    vs[i].Coefs[k, 0] -= uij * vs[j].Coefs[k, 0];
            }
        }

        return (KMatrix<K>.MergeSameRows(vs), new(u));
    }

    public static KMatrix<K> ToMatrix<K>(this KPoly<K> f, int n) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (n < f.Degree)
            throw new ArgumentException();

        var mat = new KMatrix<K>(f.KZero, n + 1, n + 1);
        for (int i = 0; i <= n; i++)
        {
            mat.Coefs[i, n] = f[i];
            if (i < n)
                mat.Coefs[1 + i, i] = f.KOne;
        }

        return mat;
    }

    public static KMatrix<K> ToHMatrix<K>(this KPoly<K> f, int n) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (n < f.Degree)
            throw new ArgumentException();

        var mat = new KMatrix<K>(f.KZero, 1, n + 1);
        for (int i = 0; i <= n; i++)
            mat.Coefs[0, i] = f[i];

        return mat;
    }

    public static KMatrix<K> ToVMatrix<K>(this KPoly<K> f, int n) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (n < f.Degree)
            throw new ArgumentException();

        var mat = new KMatrix<K>(f.KZero, n + 1, 1);
        for (int i = 0; i <= n; i++)
            mat.Coefs[i, 0] = f[i];

        return mat;
    }

    public static KMatrix<K> ToHMatrix<K>(this KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return f.ToHMatrix(f.Degree);
    }

    public static KMatrix<K> ToVMatrix<K>(this KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return f.ToVMatrix(f.Degree);
    }

    public static KMatrix<K> ToMatrix<K>(this KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return f.ToMatrix(f.Degree);
    }

    public static KMatrix<FracPoly<K>> ToKPolyMatrix<K>(this FracPoly<K> f0, int n)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (!f0.Denom.Equals(f0.One.Num) || n < f0.Num.Degree)
            throw new ArgumentException();

        var f = f0.Num;
        var mat = new KMatrix<FracPoly<K>>(f0.Zero, n + 1, n + 1);
        for (int i = 0; i <= n; i++)
        {
            mat.Coefs[i, n] = f[i] * f0.One;
            if (i < n)
                mat.Coefs[1 + i, i] = f0.One;
        }

        return mat;
    }
    
    public static KMatrix<K> Conj<K>(this KMatrix<K> A)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
    {
        return A.Select(e => e.Conj).ToKMatrix(A.M).T;
    }

}