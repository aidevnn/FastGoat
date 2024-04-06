using System.Collections;
using System.ComponentModel;
using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

MatFq TransposeConj(MatFq m, EPoly<ZnInt> ax)
{
    var table = m.Table.ToArray();
    var n = m.GLnq.N;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    {
        var u = m.Table[i * n + j];
        var v = u.Substitute(ax);
        table[j * n + i] = v;
        // Console.WriteLine(new { i, j, u, v });
    }
    
    return new(m.GLnq, table);
}

ConcreteGroup<MatFq> U2p(int p)
{
    if (p < 2 || !Primes10000.Contains(p))
        throw new();
    
    var n = 2;
    var p2 = p * p;
    var Glnq = new GLnq(2, p2);
    
    var a = FG.FqX(p2, 'a');
    var X = FG.KPoly('X', a);
    var eq = a.F.Substitute(X) / (X - a);
    var ax = -eq[0];

    var (_, x1, x2, x3, x4, x5, x6, x7, x8) = Ring.Polynomial(a, MonomOrder.Lex, (9, "x")).Deconstruct();
    var m = new[] { x1 + a * x2, x3 + a * x4, x5 + a * x6, x7 + a * x8 }.ToKMatrix(2);
    var mc = new[] { x1 + ax * x2, x5 + ax * x6, x3 + ax * x4, x7 + ax * x8 }.ToKMatrix(2);
    var j = new[] { x1.Zero, x1.One, x1.One, x1.Zero }.ToKMatrix(n);
    var P = mc * j * m - j;

    var gens = new List<MatFq>();
    if (p == 2)
        gens.Add(Glnq[0, 1, 1, 1]);
    
    {
        var subs = new[]
        {
            (x1.ExtractIndeterminate, x1.Zero),
            (x2.ExtractIndeterminate, x1.One),
            (x3.ExtractIndeterminate, x1.Zero),
            (x4.ExtractIndeterminate, x1.Zero),
            (x5.ExtractIndeterminate, x1.Zero),
            (x6.ExtractIndeterminate, x1.Zero),
            (x7.ExtractIndeterminate, x1.Zero)
        }.ToList();
        var P1 = P.Select(c => c.Substitute(subs)).ToKMatrix(n);
        var arrP1 = P1.Where(c => !c.IsZero()).ToArray();
        var sol = Ring.ReducedGrobnerBasis(arrP1)[0];
    
        var c0 = -sol.ConstTerm;
        var gi = Glnq[a, 0, 0, c0 * a];
        gens.Add(gi);
    }

    {
        var subs = new[]
        {
            (x1.ExtractIndeterminate, x1.Zero),
            (x2.ExtractIndeterminate, x1.Zero),
            (x3.ExtractIndeterminate, x1.One),
            (x4.ExtractIndeterminate, x1.Zero),
            (x5.ExtractIndeterminate, x1.One),
            (x6.ExtractIndeterminate, x1.Zero)
        }.ToList();
        var P1 = P.Select(c => c.Substitute(subs)).ToKMatrix(n);
        var arrP1 = P1.Where(c => !c.IsZero()).ToArray();
        var sol = Ring.ReducedGrobnerBasis(arrP1)[0];
        
        var e = sol[x7.ExtractMonom];
        var f = sol[x8.ExtractMonom];
        var c0 = -(p - 1) * e / f;
        var gi = Glnq[0, 1, 1, (p - 1) * e + c0 * a];
        gens.Add(gi);
    }
    
    var unp = Group.Generate($"U(2,{p})", Glnq, gens.ToArray());
    // var J = Glnq[0, 1, 1, 0];
    // Console.WriteLine(unp.All(A => J.Equals(Glnq.Op(TransposeConj(A, ax), Glnq.Op(J, A)))));
    return unp;
}

{
    GlobalStopWatch.Restart();
    foreach (var p in Primes10000.Take(7)) // p in { 2,3,5,7,11,13,17 }
    {
        GlobalStopWatch.AddLap();
        var Unp = U2p(p);
        Console.WriteLine(Unp.ShortName);
        if (Unp.Count() < 100)
        {
            var sg = Unp.AllSubgroups();
            var idGr = FG.FindIdGroup(Unp, sg.Infos)[0];
            Console.WriteLine(idGr.FullName);
        }
        
        GlobalStopWatch.Show(Unp.Name);
        Console.WriteLine();
    }
    
    Console.Beep();
}
/* |U(2,2)| = 18
   SmallGroup(18,3) Name:C3 x S3
   # U(2,2) Time:20ms
   
   |U(2,3)| = 96
   SmallGroup(96,67) Name:SL(2,3) : C4
   # U(2,3) Time:351ms
   
   |U(2,5)| = 720
   # U(2,5) Time:66ms
   
   |U(2,7)| = 2688
   # U(2,7) Time:541ms
   
   |U(2,11)| = 15840
   # U(2,11) Time:14.788s
   
   |U(2,13)| = 30576
   # U(2,13) Time:1m6s
   
   |U(2,17)| = 88128
   # U(2,17) Time:10m17s

*/