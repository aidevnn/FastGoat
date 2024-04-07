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

MatFq Transpose(MatFq m)
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

MatFq MulTransp(MatFq m) => m.GLnq.Op(m, Transpose(m));

MatFq TranspMul(MatFq m) => m.GLnq.Op(Transpose(m), m);

ConcreteGroup<MatFq> UnitaryGroup(int p, int ord, bool special = false)
{
    if (!Primes10000.Contains(p) || p > 17)
        throw new();

    var n = 2;
    var q2 = p * p;
    var Glnq = new GLnq(n, q2);
    var a = Glnq.Fq.X;
    var ax = a.F[0] / a;
    var arrFq = Group.MulGroup($"F{q2}", a).Prepend(a.Zero).ToArray();
    var J = Glnq[0, 1, 1, 0];
    
    MatFq Prod(MatFq m) => Glnq.Op(Adjoint(m, ax), Glnq.Op(J, m));
    var gen0 = arrFq.Grid2D().Where(x => !x.t1.IsZero() && !x.t2.IsZero())
        .Select(x => !special ? Glnq[a, 0, 0, x.t1] : Glnq[0, x.t1, x.t2, 0])
        .First(m => !m.Equals(Glnq.Neutral()) && Prod(m).Equals(J) && (!special || Glnq.Determinant(m).Equals(a.One)));
    var gen1 = arrFq.Grid2D().Where(x => !x.t1.IsZero() && !x.t2.IsZero())
        .Select(x => !special ? Glnq[0, 1, 1, x.t1] : Glnq[0, x.t1, x.t2, 1])
        .First(m => Prod(m).Equals(J) && (!special || Glnq.Determinant(m).Equals(a.One)));

    var gens = new HashSet<MatFq>() { gen0, gen1 };
    if (p == 2)
        gens.Add(Glnq[1, 0, 1, 1]);

    var name = special ? $"SU(2,{p})" : $"GU(2,{p})";
    var group = Group.Generate(name, Glnq, [..gens]);
    
    var check1 = group.GetGenerators().All(m => Prod(m).Equals(J));
    if (special)
    {
        var check2 = group.GetGenerators().All(m => Glnq.Determinant(m).Equals(a.One));
        Console.WriteLine("All A in {0}, bA*J*A=J {1}, Det A=1 {2}", group.Name, check1, check2);
    }
    else
        Console.WriteLine("All A in {0}, bA*J*A=J {1}", group.Name, check1);
    
    group.GetGenerators().Select(m => $"Glnq[{m.Table.Glue(",")}]").Println($"Generators of {group.ShortName}");
    Console.WriteLine();
    
    if (group.Count() != ord)
        throw new ($"############ expected {ord}");
    
    return group;
}

ConcreteGroup<MatFq> OrthogonalGroup(int q, int ord, bool special = false)
{
    if (q < 2 || PrimesDec(q).Count != 1)
        throw new();

    var n = 3;
    var Glnq = new GLnq(n, q);
    var a = Glnq.Fq.X;
    var arrFq = Group.MulGroup($"F{q}", a).Prepend(a.Zero).ToArray();
    var randGen = (int t0 = 2) =>
    {
        while (true)
        {
            EPoly<ZnInt>[] m;
            var (z, o) = (a.Zero, a.One);
            if (t0 == 0 || t0 == 1)
            {
                var t1 = 1 - t0;
                m = q > 7
                    ? 6.Range().Select(_ => arrFq[Rng.Next(0, q)]).Concat([o * t0, o * t1, z]).ToArray()
                    : 9.Range().Select(_ => arrFq[Rng.Next(0, q)]).ToArray();
            }
            else
            {
                var (r0, r1, r2) = (arrFq[Rng.Next(1, q)], arrFq[Rng.Next(1, q)], arrFq[Rng.Next(1, q)]);
                m = new[] { r0, z, z, z, z, r1, z, r2, z }.ToArray();
            }
            
            var m0 = new MatFq(Glnq, m);
            if (MulTransp(m0).Equals(Glnq.Neutral()))
            {
                if (!special || Glnq.Determinant(m0).Equals(a.One))
                    return m0;
            }
        }
    };

    var gens0 = new List<MatFq>() { randGen(), randGen(0), randGen(1) };
    var set0 = Group.GenerateElements(Glnq, gens0.ToArray());
    var set = new GroupSubset<MatFq>(gens0.ToHashSet(), set0);
    while (set.Count < ord)
    {
        var sz = set.Count;
        gens0 = set.Generators.Concat([randGen(), randGen(0), randGen(1)]).ToList();
        var elts = Group.GenerateElements(Glnq, set.Elements, gens0);
        if (sz < elts.Count)
            set = new(gens0.ToHashSet(), elts);
    }

    if (!set.Generators.All(m => TranspMul(m).Equals(Glnq.Neutral()) && MulTransp(m).Equals(Glnq.Neutral())))
        throw new();

    var name = special ? $"SO(3,{q})" : $"O(3,{q})";
    var group = Group.Generate(name, Glnq, set.Generators.ToArray());
    var check1 = group.All(m => TranspMul(m).Equals(Glnq.Neutral()) && MulTransp(m).Equals(Glnq.Neutral()));
    if (special)
    {
        var check2 = group.All(m => Glnq.Determinant(m).Equals(a.One));
        Console.WriteLine("All A in {0}, A*AT=I {1}, Det A=1 {2}", group.Name, check1, check2);
    }
    else
        Console.WriteLine("All A in {0}, A*AT=I {1}", group.Name, check1);
    
    group.GetGenerators().Select(m => $"Glnq[{m.Table.Glue(",")}]").Println($"Generators of {group.ShortName}");
    Console.WriteLine();

    return group;
}

void Test_GU_SU()
{
    var infos = new[]
    {
        (2, 18, false), (2, 6, true),
        (3, 96, false), (3, 24, true),
        (5, 720, false), (5, 120, true),
        (7, 2688, false), (7, 336, true),
        (11, 15840, false), (11, 1320, true),
        (13, 30576, false), (13, 2184, true),
        // (17, 88128, false), (17, 4896, true), // 10mn
    };
    foreach (var (q, ord, isSpecial) in infos)
    {
        GlobalStopWatch.AddLap();
        var group = UnitaryGroup(q, ord, isSpecial);
        GlobalStopWatch.Show(group.Name);
        Console.WriteLine();
    }

    Console.Beep();
}

void Test_GO_SO()
{
    GlobalStopWatch.Restart();
    var infos = new[]
    {
        (2, 6, false), (2, 6, true),
        (4, 60, false), (4, 60, true),
        (8, 504, false), (8, 504, true),
        (3, 48, false), (3, 24, true),
        (9, 1440, false), (9, 720, true),
        (5, 240, false), (5, 120, true),
        (7, 672, false), (7, 336, true),
        (11, 2640, false), (11, 1320, true),
        (13, 4368, false), (13, 2184, true),
        (17, 9792, false), (17, 4896, true),
    };
    foreach (var (q, ord, isSpecial) in infos)
    {
        GlobalStopWatch.AddLap();
        var group = OrthogonalGroup(q, ord, isSpecial);
        GlobalStopWatch.Show(group.Name);
        Console.WriteLine();
    }

    Console.Beep();
}

{
    Test_GU_SU();
    Test_GO_SO();
}