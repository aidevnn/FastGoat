using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class SymbolicMatrix
{
    static SymbolicMatrix()
    {
        Monom.ShowStar = true;
    }

    public static void Mat2x2()
    {
        var coefs = Ring.Polynomial(ZnInt.KZero(), "abcd".ToArray());
        var z0 = Ring.PolynomialZero(ZnInt.KZero());
        var mat = Ring.Matrix(2, coefs);

        Console.WriteLine("Matrix");
        Ring.DisplayMatrix(mat);

        Console.WriteLine("CoMatrix");
        Ring.DisplayMatrix(Ring.CoMatrix(mat, z0));
        Console.WriteLine("Det = {0}", Ring.Determinant(mat, z0));
        Console.WriteLine();
    }

    public static void Mat2x2bis()
    {
        var coefs = Ring.Polynomial(ZnInt.KZero(), "abba".ToArray());
        var z0 = Ring.PolynomialZero(ZnInt.KZero());
        var mat = Ring.Matrix(2, coefs);

        Console.WriteLine("Matrix");
        Ring.DisplayMatrix(mat);

        Console.WriteLine("CoMatrix");
        Ring.DisplayMatrix(Ring.CoMatrix(mat, z0));
        Console.WriteLine("Det = {0}", Ring.Determinant(mat, z0));
        Console.WriteLine();
    }

    public static void Mat3x3()
    {
        var coefs = Ring.Polynomial(ZnInt.KZero(), "abcdefghi".ToArray());
        var z0 = Ring.PolynomialZero(ZnInt.KZero());
        var mat = Ring.Matrix(3, coefs);

        Console.WriteLine("Matrix");
        Ring.DisplayMatrix(mat);

        Console.WriteLine("CoMatrix");
        Ring.DisplayMatrix(Ring.CoMatrix(mat, z0));
        Console.WriteLine("Det = {0}", Ring.Determinant(mat, z0));
        Console.WriteLine();
    }

    public static void Mat3x3bis()
    {
        var coefs = Ring.Polynomial(ZnInt.KZero(), "abcbaecea".ToArray());
        var z0 = Ring.PolynomialZero(ZnInt.KZero());
        var mat = Ring.Matrix(3, coefs);

        Console.WriteLine("Matrix");
        Ring.DisplayMatrix(mat);

        Console.WriteLine("CoMatrix");
        Ring.DisplayMatrix(Ring.CoMatrix(mat, z0));
        Console.WriteLine("Det = {0}", Ring.Determinant(mat, z0));
        Console.WriteLine();
    }

    public static void Mat4x4()
    {
        var coefs = Ring.Polynomial(ZnInt.KZero(), "abcdefghijklmnop".ToArray());
        var z0 = Ring.PolynomialZero(ZnInt.KZero());
        var mat = Ring.Matrix(4, coefs);

        Console.WriteLine("Matrix");
        Ring.DisplayMatrix(mat);

        Console.WriteLine("CoMatrix");
        Ring.DisplayMatrix(Ring.CoMatrix(mat, z0), ",\n    ");
        Console.WriteLine("Det = {0}", Ring.Determinant(mat, z0));
        Console.WriteLine();
    }

    public static void Mat4x4bis()
    {
        var coefs = Ring.Polynomial(ZnInt.KZero(), "abcdbaefceagdfga".ToArray());
        var z0 = Ring.PolynomialZero(ZnInt.KZero());
        var mat = Ring.Matrix(4, coefs);

        Console.WriteLine("Matrix");
        Ring.DisplayMatrix(mat);

        Console.WriteLine("CoMatrix");
        Ring.DisplayMatrix(Ring.CoMatrix(mat, z0), ",\n    ");
        Console.WriteLine("Det = {0}", Ring.Determinant(mat, z0));
        Console.WriteLine();
    }

    public static void Mat5x5()
    {
        var coefs = Ring.Polynomial(ZnInt.KZero(), "abcdefghijklmnopqrstuvwxy".ToArray());
        var z0 = Ring.PolynomialZero(ZnInt.KZero());
        var mat = Ring.Matrix(5, coefs);

        Console.WriteLine("Matrix");
        Ring.DisplayMatrix(mat);

        Console.WriteLine("CoMatrix");
        Ring.DisplayMatrix(Ring.CoMatrix(mat, z0), ",\n    ");
        Console.WriteLine("Det = {0}", Ring.Determinant(mat, z0));
        Console.WriteLine();
    }

    public static void SylvesterMatrix()
    {
        var x = Ring.Polynomial(ZnInt.KZero());
        var X = x.Indeterminates.First();
        var f = 2 * x.Pow(2) + 3 * x + 1;
        var g = 7 * x.Pow(2) + x + 3;
        var S = Ring.SylvesterMatrix(f, X, g, X);
        Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero));
        Ring.DisplayMatrix(S);
    }

    public static void QuadraticDiscriminant()
    {
        var x = Ring.Polynomial(ZnInt.KZero());
        var (a, b, c) = Ring.Polynomial('a', 'b', 'c', ZnInt.KZero());
        var X = x.Indeterminates.First();
        var f = a * x.Pow(2) + b * x + c;
        var g = f.D(X);
        var S = Ring.SylvesterMatrix(f, X, g, X);
        var am = f.CoefMax(X);
        var n = f.DegreeOf(X);
        var s = (int)Math.Pow(-1, n * (n - 1) / 2);
        
        Console.WriteLine($"f({X})  = {f}");
        Console.WriteLine($"f'({X}) = {g}");
        Console.WriteLine("Sylvester Matrix f, f'");
        Ring.DisplayMatrix(S);
        Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero).Div(am).quo.Mul(s));
        
        // Direct method
        Console.WriteLine("Disc = {0}", Ring.Discriminant(f, X));
    }

    public static void CubicDiscriminant()
    {
        var x = Ring.Polynomial(ZnInt.KZero());
        var (p, q) = Ring.Polynomial('p', 'q', ZnInt.KZero());
        var X = x.Indeterminates.First();
        var f = x.Pow(3) + p * x + q;
        var g = f.D(X);
        var S = Ring.SylvesterMatrix(f, X, g, X);
        var am = f.CoefMax(X);
        var n = f.DegreeOf(X);
        var s = (int)Math.Pow(-1, n * (n - 1) / 2);
        
        Console.WriteLine($"f({X})  = {f}");
        Console.WriteLine($"f'({X}) = {g}");
        Console.WriteLine("Sylvester Matrix f, f'");
        Ring.DisplayMatrix(S);
        Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero).Div(am).quo.Mul(s));
        
        // Direct method
        Console.WriteLine("Disc = {0}", Ring.Discriminant(f, X));
    }
}