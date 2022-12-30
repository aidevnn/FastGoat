using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class SymbolicMatrix
{
    static SymbolicMatrix()
    {
        Monom.Display = MonomDisplay.StarSuperscript;
    }

    public static void Mat2x2()
    {
        var coefs = Ring.Polynomial(ZnInt.KZero(), "abcd".ToArray());
        var z0 = coefs[0].Zero;
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
        var z0 = coefs[0].Zero;
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
        var z0 = coefs[0].Zero;
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
        var z0 = coefs[0].Zero;
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
        var z0 = coefs[0].Zero;
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
        var z0 = coefs[0].Zero;
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
        var z0 = coefs[0].Zero;
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
        var f = 2 * x.Pow(2) + 3 * x + 1;
        var g = 7 * x.Pow(2) + x + 3;
        Console.WriteLine(new { f, g });
        var X = x.Indeterminates[0];
        var S = Ring.SylvesterMatrix(f, X, g, X);
        Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero));
        Ring.DisplayMatrix(S);
    }

    public static void QuadraticDiscriminant()
    {
        var (x, a, b, c) = Ring.Polynomial("x", "a", "b", "c", ZnInt.KZero());
        var X = x.Indeterminates.First();
        var f = a * x.Pow(2) + b * x + c;
        var g = f.D(X);
        var S = Ring.SylvesterMatrix(f, X, g, X);
        var am = f.CoefMax(X);
        var n = f.DegreeOf(X);
        var s = (n * (n - 1) / 2) % 2 == 0 ? 1 : -1;

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
        var (x, p, q) = Ring.Polynomial("x", "p", "q", ZnInt.KZero());
        var X = x.Indeterminates.First();
        var f = x.Pow(3) + p * x + q;
        var g = f.D(X);
        var S = Ring.SylvesterMatrix(f, X, g, X);
        var am = f.CoefMax(X);
        var n = f.DegreeOf(X);
        var s = (n * (n - 1) / 2) % 2 == 0 ? 1 : -1;

        Console.WriteLine($"f({X})  = {f}");
        Console.WriteLine($"f'({X}) = {g}");
        Console.WriteLine("Sylvester Matrix f, f'");
        Ring.DisplayMatrix(S);
        Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero).Div(am).quo.Mul(s));

        // Direct method
        Console.WriteLine("Disc = {0}", Ring.Discriminant(f, X));
        // Disc = -4*p³ + -27*q²
    }

    public static void CubicDiscriminantLong()
    {
        var ps = Ring.Polynomial(ZnInt.KZero(), "xabcd".ToArray());
        var (x, a, b, c, d) = (ps[0], ps[1], ps[2], ps[3], ps[4]);
        var X = x.Indeterminates[0];
        var f = a * x.Pow(3) + b * x.Pow(2) + c * x + d;
        var g = f.D(X);
        var S = Ring.SylvesterMatrix(f, X, g, X);
        var am = f.CoefMax(X);
        var n = f.DegreeOf(X);
        var s = (n * (n - 1) / 2) % 2 == 0 ? 1 : -1;

        Console.WriteLine($"f({X})  = {f}");
        Console.WriteLine($"f'({X}) = {g}");
        Console.WriteLine("Sylvester Matrix f, f'");
        Ring.DisplayMatrix(S);
        Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero).Div(am).quo.Mul(s));

        // Direct method
        Console.WriteLine("Disc = {0}", Ring.Discriminant(f, X));

        // https://en.wikipedia.org/wiki/Discriminant#Degree_3
        // Disc = -27·a²*d² + 18·a*b*c*d + -4·a*c³ + -4·b³*d + b²*c²
    }

    public static void QuarticDiscriminant()
    {
        var ps = Ring.Polynomial(ZnInt.KZero(), "xabcde".ToArray());
        var (x, a, b, c, d, e) = (ps[0], ps[1], ps[2], ps[3], ps[4], ps[5]);
        var X = x.Indeterminates.First();
        var f = a * x.Pow(4) + b * x.Pow(3) + c * x.Pow(2) + d * x + e;
        var g = f.D(X);
        var S = Ring.SylvesterMatrix(f, X, g, X);
        var am = f.CoefMax(X);
        var n = f.DegreeOf(X);
        var s = (n * (n - 1) / 2) % 2 == 0 ? 1 : -1;

        Console.WriteLine($"f({X})  = {f}");
        Console.WriteLine($"f'({X}) = {g}");
        Console.WriteLine("Sylvester Matrix f, f'");
        Ring.DisplayMatrix(S);
        Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero).Div(am).quo.Mul(s));

        // Direct method
        Console.WriteLine("Disc = {0}", Ring.Discriminant(f, X));

        // https://en.wikipedia.org/wiki/Discriminant#Degree_4
        // Disc = 256·a³*e³ + -192·a²*b*d*e² + -128·a²*c²*e² + 144·a²*c*d²*e +
        // -27·a²*d⁴ + 144·a*b²*c*e² + -6·a*b²*d²*e + -80·a*b*c²*d*e +
        // 18·a*b*c*d³ + 16·a*c⁴*e + -4·a*c³*d² + -27·b⁴*e² + 18·b³*c*d*e +
        // -4·b³*d³ + -4·b²*c³*e + b²*c²*d²
    }
}