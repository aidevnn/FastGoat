using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Linq.Expressions;
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
using System.Threading.Channels;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Random rnd = new Random();

{
    var m = Ring.Matrix(2, Rational.KZero(), 3, 2, 1, 2);
    var A = new KMatrix<Rational>(m);
    var (O, U) = Ring.GramSchmidt(A);
    Console.WriteLine("Matrix A");
    Console.WriteLine(A);
    Console.WriteLine("Matrix O");
    Console.WriteLine(O);
    Console.WriteLine("Matrix U");
    Console.WriteLine(U);
    Console.WriteLine("(U * O.T).T");
    Console.WriteLine((U * O.T).T);
}

{
    var m = Ring.Matrix(2, Rational.KZero(), 6, 8, 4, 4);
    var A = new KMatrix<Rational>(m);
    var (O, U) = Ring.GramSchmidt(A);
    Console.WriteLine("Matrix A");
    Console.WriteLine(A);
    Console.WriteLine("Matrix O");
    Console.WriteLine(O);
    Console.WriteLine("Matrix U");
    Console.WriteLine(U);
    Console.WriteLine("(U * O.T).T");
    Console.WriteLine((U * O.T).T);
}

{
    var n = 3;
    var m = Ring.Matrix(n, Rational.KZero(), 1, 2, 3, 4, 7, 6, 2, 1, 4);
    var A = new KMatrix<Rational>(m);
    var (O, U) = Ring.GramSchmidt(A);
    Console.WriteLine("Matrix A");
    Console.WriteLine(A);
    Console.WriteLine("Matrix O");
    Console.WriteLine(O);
    Console.WriteLine("Matrix U");
    Console.WriteLine(U);
    Console.WriteLine("(U * O.T).T");
    Console.WriteLine((U * O.T).T);
}

KMatrix<Rational> RandMatrixRational(int n)
{
    var m = new KMatrix<Rational>(Rational.KZero(), n, n);
    while (true)
    {
        var m0 = Ring.Matrix(n, Rational.KZero(), (n * n).Range().Select(i => rnd.Next(n + 1)).ToArray());
        m = new(m0);
        if (!m.Det.IsZero())
            return m;
    }
}

for (int k = 0; k < 15; ++k)
{
    var n = 2 + rnd.Next(15);
    GlobalStopWatch.Restart();
    var A = RandMatrixRational(n);
    var (O, U) = Ring.GramSchmidt(A);
    Console.WriteLine($"Dim {A.Dim}");
    Console.WriteLine("U * O.T = A.T : {0}", (U * O.T).Equals(A.T));
    Console.WriteLine("All Orthogonals : {0}", n.Range().Grid2D(n.Range())
        .All(e => e.t1.Equals(e.t2) || (O.GetCol(e.t1).T * O.GetCol(e.t2)).IsZero()));

    GlobalStopWatch.Show("GS2");
    Console.WriteLine();
}
