using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using Microsoft.VisualBasic.CompilerServices;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly();
    
    var P = x.Pow(6) + 12;
    var rootsK = IntFactorisation.AlgebraicRoots(P);
    var gal = GaloisTheory.GaloisGroup(rootsK, details: true);
    GaloisApplicationsPart2.CheckChebotarev(P, gal, 100);
    
    var rootsC = FG.NRoots(P.ToCPoly()).ToList();
    rootsC.Println();
    var r0 = rootsC[0]; // any index
    rootsC = rootsK.Select(r => r.Poly.Substitute(r0)).ToList(); // mandatory
    
    Console.WriteLine(rootsC.Aggregate(Cplx.X.One, (prod, c1) => prod * (Cplx.X - c1)));
    
    Console.WriteLine();
    foreach (var perm in gal)
    {
        var n = perm.Sn.N;
        var r = rootsK.First(r => perm.Table.Select((j, i) => (i, j)).All(e => r.Substitute(rootsK[e.i]).Equals(rootsK[e.j])));
        perm.Table.Select((j, i) => (rootsC[i], rootsC[j], r.Poly.Substitute(rootsC[i]).Equals(rootsC[j])))
            .Println($"r = {r} perm = {perm}");
    
        var mat0 = perm.Table.Select((j, i) => n.Range().Select(k => rootsC[i].Pow(k)).Append(rootsC[j]).ToKMatrix()).ToArray();
        var mat = KMatrix<Cplx>.MergeSameCols(mat0);
        var (MatP, A0) = Ring.ReducedRowsEchelonForm(mat);
        var rc = A0.Cols.Last().ToKPoly('y');
        perm.Table.Select((j, i) => (rootsC[i], rootsC[j], rc.Substitute(rootsC[i]), rc.Substitute(rootsC[i]).Equals(rootsC[j])))
            .Println($"rc = {rc}");
        Console.WriteLine();
    }
}
