using System.Collections;
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

MatFq MulTransp(MatFq m) => m.GLnq.Op(m, Transpose(m));

MatFq TranspMul(MatFq m) => m.GLnq.Op(Transpose(m), m);

MatFq SelfAdjoint(MatFq m, EPoly<ZnInt> ax)
{
    var table = m.Table.ToArray();
    var n = m.GLnq.N;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        table[j * n + i] = m.Table[i * n + j].Substitute(ax);

    return new(m.GLnq, table);
}

ConcreteGroup<MatFq> UnitaryGroup(int q, bool special = false)
{
    var dec = PrimesDec(q);
    if (dec.Count > 2 || q > 17)
        throw new();

    var n = 2;
    var q2 = q * q;
    var Glnq = new GLnq(n, q2);
    var a = Glnq.Fq.X;
    var arrFq = Group.MulGroup($"F{q2}", a).Prepend(a.Zero).ToArray();

    // conj(conj(x))=x then ax(ax(a)) = a
    var ax = arrFq.First(e => !e.IsZero() && !e.Equals(a) && a.F.Substitute(e).IsZero() && e.Substitute(e).Equals(a));

    var J = Glnq[0, 1, 1, 0];
    MatFq Prod(MatFq m) => Glnq.Op(SelfAdjoint(m, ax), Glnq.Op(J, m));

    var gens = arrFq.Grid2D().SelectMany(x => new[]
        {
            Glnq[0, 1, 1, x.t2],
            Glnq[0, x.t1, x.t2, 1],
            Glnq[a, 0, 0, x.t2]
        })
        .Where(m => Prod(m).Equals(J) && (!special || Glnq.Determinant(m).Equals(a.One)))
        .ToHashSet();

    var name = special ? $"SU(2,{q})" : $"GU(2,{q})";
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

    return group;
}

ConcreteGroup<MatFq> OrthogonalGroup(int q, bool special = false)
{
    if (q < 2 || PrimesDec(q).Count != 1 || q > 19)
        throw new();

    var n = 3;
    var Glnq = new GLnq(n, q);
    var a = Glnq.Fq.X;
    var arrFq = Group.MulGroup($"F{q}", a).Prepend(a.Zero).ToArray();

    var gen0 = arrFq.Grid3D().Where(x => !x.t1.IsZero() && !x.t2.IsZero() && !x.t3.IsZero())
        .SelectMany(x => new[]
        {
            Glnq[x.t1, 0, 0, 0, 0, x.t2, 0, x.t3, 0],
            Glnq[0, x.t1, 0, 0, 0, x.t2, x.t3, 0, 0],
            Glnq[0, x.t1, x.t2, 0, x.t3, x.t1, 1, 0, 0],
            Glnq[x.t1, x.t2, 0, 0, x.t3, x.t1, 1, 0, 0],
            Glnq[x.t1, x.t2, 1, 1, 1, x.t3, x.t2, x.t1, 1]
        })
        .Where(m => MulTransp(m).Equals(Glnq.Neutral()) && (!special || Glnq.Determinant(m).Equals(a.One)))
        .OrderBy(m => m, Comparer<MatFq>.Create((m0, m1) => m0.Table.SequenceCompareTo(m1.Table)))
        .ToHashSet();

    if (!gen0.All(m => TranspMul(m).Equals(Glnq.Neutral()) && MulTransp(m).Equals(Glnq.Neutral())))
        throw new();

    var name = special ? $"SO(3,{q})" : $"O(3,{q})";
    var group = Group.Generate(name, Glnq, gen0.ToArray());
    var check1 = group.GetGenerators().All(m => TranspMul(m).Equals(Glnq.Neutral()) && MulTransp(m).Equals(Glnq.Neutral()));
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

void Test_GO_SO()
{
    GlobalStopWatch.Restart();
    foreach (var q in new[] { 2, 4, 8, 16, 3, 9, 5, 7, 11, 13, 17 })
    {
        foreach (var isSpecial in new[] { false, true })
        {
            GlobalStopWatch.AddLap();
            var group = OrthogonalGroup(q, isSpecial);
            GlobalStopWatch.Show(group.Name);
            Console.WriteLine();
        }
    }

    Console.Beep();
}

void Test_GU_SU()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    foreach (var q in new[] { 2, 4, 8, 3, 9, 5, 7, 11 })
    {
        foreach (var isSpecial in new[] { false, true })
        {
            GlobalStopWatch.AddLap();
            var group = UnitaryGroup(q, isSpecial);
            GlobalStopWatch.Show(group.Name);
            Console.WriteLine();
        }
    }

    Console.Beep();
}

{
    Test_GO_SO();
}