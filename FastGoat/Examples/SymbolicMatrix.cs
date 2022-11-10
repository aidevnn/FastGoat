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
}