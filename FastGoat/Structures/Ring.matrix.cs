using FastGoat.Commons;

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

    public static T[,] Matrix<T>(int rows, T t0, params int[] coefs) where T : IRingElt<T>, IElt<T>
    {
        return Matrix(rows, coefs.Select(c => t0.One.Mul(c)).ToArray());
    }

    public static T[,] Diagonal<T>(T v, int n) where T : IElt<T>, IRingElt<T>
    {
        var coefs = n.Range().Grid2D().Select(e => e.t1 == e.t2 ? v : v.Zero).ToArray();
        return Matrix(n, coefs);
    }

    public static T[,] Dot<T>(T[,] A, T[,] B, T t0) where T : IRingElt<T>, IElt<T>
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
}