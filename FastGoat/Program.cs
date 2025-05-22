using System.Numerics;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using static FastGoat.Commons.IntExt;


//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

{
    var (p, n) = (2, 4);
    var q = BigInteger.Pow(p, n);
    var a = FG.FqX(q, 'a');
    var g = NumberTheory.PrimitiveRoot(a);
    foreach (var k in DividorsBigInt(q - 1).Select(k => (int)k).Order().ToArray())
    {
        for (int i = 1; i < q; i++)
        {
            var e = g.Pow(i);
            if (!e.FastPow((q - 1) / k).IsOne())
                continue;
            var kthRoot = NumberTheory.NthRoot(e, k, g);
            var et = kthRoot.Pow(k);
            var test = et.Equals(e);
            Console.WriteLine(new { i, k, e, kthRoot });
            if (!test)
                throw new();
        }

        Console.WriteLine();
    }
}