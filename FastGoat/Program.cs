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

MatFq Transpose(MatFq m, EPoly<ZnInt> ax)
{
    var table = m.Table.ToArray();
    var n = m.GLnq.N;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        table[j * n + i] = m.Table[i * n + j];

    return new(m.GLnq, table);
}

MatFq Adjoint(MatFq m, EPoly<ZnInt> ax)
{
    var table = m.Table.ToArray();
    var n = m.GLnq.N;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        table[j * n + i] = m.Table[i * n + j].Substitute(ax);

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
    var ax = a.F[0] / a;

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
    var J = Glnq[0, 1, 1, 0];
    Console.WriteLine(unp.All(A => J.Equals(Glnq.Op(Adjoint(A, ax), Glnq.Op(J, A)))));
    return unp;
}

void TestUnp()
{
    GlobalStopWatch.Restart();
    foreach (var p in Primes10000.Take(4)) // p in { 2,3,5,7,11,13,17 }
    {
        GlobalStopWatch.AddLap();
        var Unp = U2p(p);
        Console.WriteLine(Unp.ShortName);
        // if (Unp.Count() < 100)
        // {
        //     var sg = Unp.AllSubgroups();
        //     var idGr = FG.FindIdGroup(Unp, sg.Infos)[0];
        //     Console.WriteLine(idGr.FullName);
        // }
        //
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

ConcreteGroup<MatFq> O3p(int p)
{
    if (p < 2 || PrimesDec(p).Count != 1)
        throw new();

    var n = 3;
    var Glnq = new GLnq(n, p);
    var a = Glnq.Fq.X;
    var gens = new List<MatFq>();
    var randGen1 = () =>
    {
        while (true)
        {
            var (r0, r1, r2, z) = (Rng.Next(0, p) * a.One, Rng.Next(0, p) * a.One, Rng.Next(0, p) * a.One, a.Zero);
            var m = new[] { r0, z, z, z, z, r1, z, r2, z }.ToKMatrix(3);
            var m0 = new MatFq(Glnq, m.ToArray());
            if ((m * m.T).Equals(m.One)) // (m.T * m).Equals(m.One) && 
                return m0;
        }
    };
    var randGen2 = (bool t) =>
    {
        var t0 = t ? 1 : 0;
        var t1 = t ? 0 : 1;
        while (true)
        {
            var m = p > 5
                ? 6.Range().Select(_ => a.One * Rng.Next(0, p)).Concat([a.One * t0, a.One * t1, a.Zero]).ToKMatrix(3)
                : 9.Range().Select(_ => a.One * Rng.Next(0, p)).ToKMatrix(3);
            var m0 = new MatFq(Glnq, m.ToArray());
            if ((m * m.T).Equals(m.One)) // (m.T * m).Equals(m.One) && 
                return m0;
        }
    };
    
    gens.AddRange([randGen1(), randGen2(true), randGen2(false)]);

    if (!gens.Select(m => m.Table.ToKMatrix(3)).All(m => (m.T * m).Equals(m.One) && (m * m.T).Equals(m.One)))
        throw new();

    var o3p = Group.Generate($"O(3,{p})", Glnq, gens.ToArray());
    var check = o3p.Select(m => m.Table.ToKMatrix(3)).All(m => (m.T * m).Equals(m.One) && (m * m.T).Equals(m.One));
    Console.WriteLine("All A in {0}, A*AT=I {1}", o3p.Name, check);
    o3p.GetGenerators().Select(m => $"Glnq[{m.Table.Select(e=>e[0].K).Glue(",")}]").Println($"Generators of {o3p.ShortName}");
    Console.WriteLine();

    return o3p;
}

{
    // TestUnp();
    GlobalStopWatch.Restart();
    var ps = Primes10000.Take(7).ToArray();
    var ords = new[] { 6, 48, 240, 672, 2640, 4368, 9792 };
    for (int k = 0; k < 7; k++)
    {
        for (int i = 0; i < 5; i++)
        {
            GlobalStopWatch.AddLap();
            var o3p = O3p(ps[k]);
            if (o3p.Count() == ords[k])
            {
                GlobalStopWatch.Show(o3p.Name);
                Console.WriteLine();
                break;
            }
            else
            {
                GlobalStopWatch.Show($"###### |{o3p}| expected:{ords[k]} actual:{o3p.Count()} ######");
                Console.WriteLine();
            }
        }
    }

    Console.Beep();
}

/* All A in O(3,2), A*AT=I True
   Generators of |O(3,2)| = 6
       Glnq[0,0,1,0,1,0,1,0,0]
       Glnq[0,1,0,0,0,1,1,0,0]
   
   # O(3,2) Time:131ms
   
   All A in O(3,3), A*AT=I True
   Generators of |O(3,3)| = 48
       Glnq[0,0,1,0,1,0,1,0,0]
       Glnq[0,0,1,0,2,0,1,0,0]
       Glnq[1,0,0,0,0,2,0,1,0]
   
   # O(3,3) Time:40ms
   
   Generators of |O(3,5)| = 240
       Glnq[1,0,0,0,0,4,0,1,0]
       Glnq[2,1,4,1,2,4,4,4,2]
   
   # O(3,5) Time:760ms
   
   All A in O(3,7), A*AT=I True
   Generators of |O(3,7)| = 672
       Glnq[0,2,2,0,2,5,1,0,0]
       Glnq[0,6,0,0,0,1,1,0,0]
   
   # O(3,7) Time:1.012s
   
   All A in O(3,11), A*AT=I True
   Generators of |O(3,11)| = 2640
       Glnq[0,8,6,0,5,8,1,0,0]
       Glnq[3,0,6,5,0,3,0,1,0]
   
   # O(3,11) Time:6.356s
   
   All A in O(3,13), A*AT=I True
   Generators of |O(3,13)| = 4368
       Glnq[0,0,12,0,12,0,1,0,0]
       Glnq[0,2,6,0,7,2,1,0,0]
   
   # O(3,13) Time:12.084s
   
   
   All A in O(3,17), A*AT=I True
   Generators of |O(3,17)| = 9792
       Glnq[0,6,4,0,13,6,1,0,0]
       Glnq[0,11,4,0,4,6,1,0,0]
   
   # O(3,17) Time:48.947s
*/