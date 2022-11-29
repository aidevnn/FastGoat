using System.Diagnostics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using System.Numerics;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

for (int n = 2; n <= 5; ++n)
{
    var alphabet0 = "abcdefghijklmnopqrstuvwxyz";
    var alphabet1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    var z0 = Ring.PolynomialZero(ZnInt.KZero());
    var coefs0 = Ring.Polynomial(ZnInt.KZero(), alphabet0.Take(n * n).ToArray());
    var coefs1 = Ring.Polynomial(ZnInt.KZero(), alphabet1.Take(n * n).ToArray());
    var mat0 = Ring.Matrix(n, coefs0);
    var mat1 = Ring.Matrix(n, coefs1);
    var mat = Ring.Dot(mat0, mat1, z0);
    Monom.Display = MonomDisplay.Star;

    Console.WriteLine($"EPoly<ZnInt>[] Dot{n}x{n}(EPoly<ZnInt>[] mat0, EPoly<ZnInt>[] mat1)\n{{");

    for (int i = 0; i < n; i++)
    {
        var i0 = i;
        var rg0 = n.Range().Select(j => mat0[i0, j]).Glue(", ");
        var rg1 = n.Range().Select(j => $"mat0[{n * i0 + j}]").Glue(", ");
        Console.WriteLine($"    var ({rg0}) = ({rg1});");
    }

    Console.WriteLine();
    for (int i = 0; i < n; i++)
    {
        var i0 = i;
        var rg0 = n.Range().Select(j => mat1[i0, j]).Glue(", ");
        var rg1 = n.Range().Select(j => $"mat1[{n * i0 + j}]").Glue(", ");
        Console.WriteLine($"    var ({rg0}) = ({rg1});");
    }

    Console.WriteLine();
    Console.WriteLine($"    var mat = new EPoly<ZnInt>[{n * n}];");

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Console.WriteLine("    mat[{0}] = {1};", n * i + j, mat[i, j]);
        }
    }

    Console.WriteLine();
    Console.WriteLine("    return mat;\n}");
    Console.WriteLine();
}
