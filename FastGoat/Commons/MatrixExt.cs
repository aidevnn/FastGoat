namespace FastGoat.Commons;

/// <summary>
/// Extension class that provides methods for manipulating matrices. 
/// </summary>
public static class MatrixExt
{
    /// <summary>
    /// Creates an matrix of size nb x nb and fills it with the value v on the diagonal.
    /// </summary>
    /// <param name="nb">Size of the matrix.</param>
    /// <param name="v">Value to fill the diagonal.</param>
    /// <returns>An array of size nb filled with v on the diagonal.</returns>
    public static int[] Diagonal(int nb, int v)
    {
        int[] mat0 = new int[nb * nb];
        for (int i = 0; i < nb; i++)
        {
            mat0[i * (nb + 1)] = v;
        }

        return mat0;
    }

    /// <summary>
    /// Creates an identity matrix of size nb x nb. 
    /// </summary>
    /// <param name="nb">The size of the matrix.</param>
    /// <returns>An identity matrix of size nb.</returns>
    public static int[] Identity(int nb) => Diagonal(nb, 1);

    /// <summary>
    /// Calculates the cofactor of a matrix at a given row and column.
    /// </summary>
    /// <param name="mat">The matrix to calculate the cofactor of.</param>
    /// <param name="row">The row index of the element to calculate the cofactor for.</param>
    /// <param name="col">The column index of the element to calculate the cofactor for.</param>
    /// <returns>An array containing the cofactor of the given element.</returns>
    public static int[] Cofactor(int[] mat, int row, int col)
    {
        int nb = (int)Math.Sqrt(mat.Length);
        int nb0 = nb - 1;
        int[] mat0 = new int[nb0 * nb0];

        for (int i = 0; i <= nb0; ++i)
        {
            if (i == row) continue;
            int i0 = i < row ? i : i - 1;
            for (int j = 0; j <= nb0; ++j)
            {
                if (j == col) continue;
                int j0 = j < col ? j : j - 1;
                mat0[i0 * nb0 + j0] = mat[i * nb + j];
            }
        }

        return mat0;
    }

    /// <summary>
    /// Calculates the determinant of a given matrix.
    /// </summary>
    /// <param name="mat">The matrix to calculate the determinant of.</param>
    /// <returns>The determinant of the given matrix.</returns>
    public static int ComputeDeterminant(int[] mat)
    {
        int nb = (int)Math.Sqrt(mat.Length);
        if (nb == 1) return mat[0];

        int det = 0;
        int s = 1;
        for (int k = 0; k < nb; ++k)
        {
            var mat0 = Cofactor(mat, k, 0);
            det += s * mat[k * nb] * ComputeDeterminant(mat0);
            s *= -1;
        }

        return det;
    }

    /// <summary>
    /// Transposes the given matrix.
    /// </summary>
    /// <param name="mat">The matrix to transpose.</param>
    /// <returns>The transposed matrix.</returns>
    public static int[] Transpose(int[] mat)
    {
        int nb = (int)Math.Sqrt(mat.Length);
        int[] mat0 = new int[mat.Length];
        for (int i = 0; i < nb; i++)
        {
            for (int j = 0; j < nb; j++)
            {
                mat0[j * nb + i] = mat[i * nb + j];
            }
        }

        return mat0;
    }

    /// <summary>
    /// Calculates the comatrix of the given matrix.
    /// </summary>
    /// <param name="mat">The matrix to calculate the comatrix of.</param>
    /// <returns>The comatrix of the given matrix.</returns>
    public static int[] Comatrix(int[] mat)
    {
        int nb = (int)Math.Sqrt(mat.Length);
        int[] mat0 = new int[mat.Length];
        for (int i = 0; i < nb; i++)
        {
            for (int j = 0; j < nb; j++)
            {
                var s = (int)Math.Pow(-1, i + j);
                mat0[i * nb + j] = s * ComputeDeterminant(Cofactor(mat, i, j));
            }
        }

        return mat0;
    }

    /// <summary>
    /// Calculates the dot product of two matrices.
    /// </summary>
    /// <param name="matA">The first matrix.</param>
    /// <param name="matB">The second matrix.</param>
    /// <returns>The dot product of the two matrices.</returns>
    public static int[] Dot(int[] matA, int[] matB)
    {
        int nb = (int)Math.Sqrt(matA.Length);
        int[] matC = new int[matA.Length];
        for (int i = 0; i < nb; i++)
        {
            for (int j = 0; j < nb; j++)
            {
                int sum = 0;
                for (int k = 0; k < nb; k++)
                {
                    sum += matA[i * nb + k] * matB[k * nb + j];
                }

                matC[i * nb + j] = sum;
            }
        }

        return matC;
    }

    /// <summary>
    /// Displays a matrix of integers.
    /// </summary>
    /// <param name="mat"></param>
    public static void DisplayMatrix(int[] mat)
    {
        int nb = (int)Math.Sqrt(mat.Length);
        for (int i = 0; i < nb; i++)
        {
            Console.WriteLine(mat.Skip(i * nb).Take(nb).Glue(" "));
        }
    }

    /// <summary>
    /// Reduce a matrix modulo p. 
    /// </summary>
    /// <param name="p">An integer representing the modulus.</param>
    /// <param name="mat">A matrix of integers.</param>
    /// <returns>An array of integers.</returns>
    public static int[] ModP(int p, int[] mat)
    {
        return mat.Select(e => IntExt.AmodP(e, p)).ToArray();
    }

    /// <summary>
    /// Calculates the dot product of two matrices and reduce the result modulo p.
    /// </summary>
    /// <param name="p">The number to use for the modulo operation.</param>
    /// <param name="matA">The first matrix to use for the dot product.</param>
    /// <param name="matB">The second matrix to use for the dot product.</param>
    /// <returns>An array containing the result of the dot product modulo p.</returns>
    public static int[] DotModP(int p, int[] matA, int[] matB) => ModP(p, Dot(matA, matB));

    /// <summary>
    /// Inverts matrix and reduce it modulo p.
    /// </summary>
    /// <param name="p">The modulo p.</param>
    /// <param name="mat">The matrix to invert.</param>
    /// <returns>The inverted matrix.</returns>
    public static int[] InvertModP(int p, int[] mat)
    {
        int n = (int)Math.Sqrt(mat.Length);
        var det = ComputeDeterminant(mat);
        var invDet = IntExt.InvModP(det, p);
        var diag = MatrixExt.Diagonal(n, invDet);

        var com = MatrixExt.Comatrix(mat);
        var tcom = MatrixExt.Transpose(com);
        var inv = MatrixExt.Dot(diag, tcom);
        var invP = MatrixExt.ModP(p, inv);

        return invP;
    }

    /// <summary>
    /// This method calculates the determinant of an NxN matrix with N from 1 to 5. 
    /// </summary>
    /// <param name="mat">The matrix to calculate the determinant of.</param>
    /// <returns>The determinant of the matrix.</returns>
    public static int DetNxN(int[] mat)
    {
        if (mat.Length == 1)
            return mat[0];
        else if (mat.Length == 4)
            return Det2x2(mat);
        else if (mat.Length == 9)
            return Det3x3(mat);
        else if (mat.Length == 16)
            return Det4x4(mat);
        else if (mat.Length == 25)
            return Det5x5(mat);

        throw new();
    }

    /// <summary>
    /// Calculates the determinant of a 2x2 matrix. 
    /// </summary>
    /// <param name="mat">The 2x2 matrix to calculate the determinant of.</param>
    /// <returns>The determinant of the given matrix.</returns>
    public static int Det2x2(int[] mat)
    {
        var (a, b) = (mat[0], mat[1]);
        var (c, d) = (mat[2], mat[3]);
        var det = (a * d - c * b);
        return det;
    }

    /// <summary>
    /// Calculates the determinant of a 3x3 matrix. 
    /// </summary>
    /// <param name="mat">The 3x3 matrix to calculate the determinant of.</param>
    /// <returns>The determinant of the given matrix.</returns>
    public static int Det3x3(int[] mat)
    {
        var (a, b, c) = (mat[0], mat[1], mat[2]);
        var (d, e, f) = (mat[3], mat[4], mat[5]);
        var (g, h, i) = (mat[6], mat[7], mat[8]);
        var det = (a * e * i + d * h * c + g * b * f - a * f * h - b * d * i - c * e * g);
        return det;
    }

    /// <summary>
    /// Calculates the determinant of a 4x4 matrix. 
    /// </summary>
    /// <param name="mat">The 4x4 matrix to calculate the determinant of.</param>
    /// <returns>The determinant of the given matrix.</returns>
    public static int Det4x4(int[] mat)
    {
        var (a, b, c, d) = (mat[0], mat[1], mat[2], mat[3]);
        var (e, f, g, h) = (mat[4], mat[5], mat[6], mat[7]);
        var (i, j, k, l) = (mat[8], mat[9], mat[10], mat[11]);
        var (m, n, o, p) = (mat[12], mat[13], mat[14], mat[15]);

        // mathematica.wolframcloud.com
        var det = (d * g * j * m - c * h * j * m - d * f * k * m + b * h * k * m + c * f * l * m - b * g * l * m -
            d * g * i * n + c * h * i * n + d * e * k * n - a * h * k * n - c * e * l * n + a * g * l * n +
            d * f * i * o - b * h * i * o - d * e * j * o + a * h * j * o + b * e * l * o - a * f * l * o -
            c * f * i * p + b * g * i * p + c * e * j * p - a * g * j * p - b * e * k * p + a * f * k * p);
        return det;
    }

    /// <summary>
    /// Calculates the determinant of a 5x5 matrix. 
    /// </summary>
    /// <param name="mat">The 5x5 matrix to calculate the determinant of.</param>
    /// <returns>The determinant of the given matrix.</returns>
    public static int Det5x5(int[] mat)
    {
        var (a, b, c, d, e) = (mat[0], mat[1], mat[2], mat[3], mat[4]);
        var (f, g, h, i, j) = (mat[5], mat[6], mat[7], mat[8], mat[9]);
        var (k, l, m, n, o) = (mat[10], mat[11], mat[12], mat[13], mat[14]);
        var (p, q, r, s, t) = (mat[15], mat[16], mat[17], mat[18], mat[19]);
        var (u, v, w, x, y) = (mat[20], mat[21], mat[22], mat[23], mat[24]);

        // mathematica.wolframcloud.com
        var det = e * i * m * q * u - d * j * m * q * u - e * h * n * q * u + c * j * n * q * u + d * h * o * q * u -
            c * i * o * q * u - e * i * l * r * u + d * j * l * r * u + e * g * n * r * u -
            b * j * n * r * u - d * g * o * r * u + b * i * o * r * u + e * h * l * s * u - c * j * l * s * u -
            e * g * m * s * u + b * j * m * s * u + c * g * o * s * u - b * h * o * s * u -
            d * h * l * t * u + c * i * l * t * u + d * g * m * t * u - b * i * m * t * u - c * g * n * t * u +
            b * h * n * t * u - e * i * m * p * v + d * j * m * p * v + e * h * n * p * v -
            c * j * n * p * v - d * h * o * p * v + c * i * o * p * v + e * i * k * r * v - d * j * k * r * v -
            e * f * n * r * v + a * j * n * r * v + d * f * o * r * v - a * i * o * r * v -
            e * h * k * s * v + c * j * k * s * v + e * f * m * s * v - a * j * m * s * v - c * f * o * s * v +
            a * h * o * s * v + d * h * k * t * v - c * i * k * t * v - d * f * m * t * v +
            a * i * m * t * v + c * f * n * t * v - a * h * n * t * v + e * i * l * p * w - d * j * l * p * w -
            e * g * n * p * w + b * j * n * p * w + d * g * o * p * w - b * i * o * p * w -
            e * i * k * q * w + d * j * k * q * w + e * f * n * q * w - a * j * n * q * w - d * f * o * q * w +
            a * i * o * q * w + e * g * k * s * w - b * j * k * s * w - e * f * l * s * w +
            a * j * l * s * w + b * f * o * s * w - a * g * o * s * w - d * g * k * t * w + b * i * k * t * w +
            d * f * l * t * w - a * i * l * t * w - b * f * n * t * w + a * g * n * t * w - e * h * l * p * x +
            c * j * l * p * x + e * g * m * p * x - b * j * m * p * x - c * g * o * p * x + b * h * o * p * x +
            e * h * k * q * x -
            c * j * k * q * x - e * f * m * q * x + a * j * m * q * x + c * f * o * q * x - a * h * o * q * x -
            e * g * k * r * x + b * j * k * r * x + e * f * l * r * x -
            a * j * l * r * x - b * f * o * r * x + a * g * o * r * x + c * g * k * t * x - b * h * k * t * x -
            c * f * l * t * x + a * h * l * t * x + b * f * m * t * x -
            a * g * m * t * x + d * h * l * p * y - c * i * l * p * y - d * g * m * p * y + b * i * m * p * y +
            c * g * n * p * y - b * h * n * p * y - d * h * k * q * y +
            c * i * k * q * y + d * f * m * q * y - a * i * m * q * y - c * f * n * q * y + a * h * n * q * y +
            d * g * k * r * y - b * i * k * r * y - d * f * l * r * y +
            a * i * l * r * y + b * f * n * r * y - a * g * n * r * y - c * g * k * s * y + b * h * k * s * y +
            c * f * l * s * y - a * h * l * s * y - b * f * m * s * y + a * g * m * s * y;
        return det;
    }
}