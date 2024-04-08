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

MatFq SelfAdjoint(MatFq m, EPoly<ZnInt> ax)
{
    var table = m.Table.ToArray();
    var n = m.GLnq.N;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        table[j * n + i] = m.Table[i * n + j].Substitute(ax);

    return new(m.GLnq, table);
}

ConcreteGroup<MatFq> UnitaryGroup(int q, int ord, bool special = false)
{
    var dec = PrimesDec(q);
    if (dec.Count > 2 || q > 17)
        throw new();

    var n = 2;
    var q2 = q * q;
    var Glnq = new GLnq(n, q2);
    var a = Glnq.Fq.X;
    var arrFq = Group.MulGroup($"F{q2}", a).Prepend(a.Zero).ToArray();
    Console.WriteLine($"arrFq^2:{arrFq.Length.Pow(2)}");
    var sols = arrFq.Where(e => !e.IsZero() && !e.Equals(a) && a.F.Substitute(e).IsZero()).ToArray();
    
    foreach (var ax in sols)
    {
        var J = Glnq[0, 1, 1, 0];
        MatFq Prod(MatFq m) => Glnq.Op(SelfAdjoint(m, ax), Glnq.Op(J, m));
        try
        {
            var gens = arrFq.Grid2D()
                .SelectMany(x => new[]
                {
                    Glnq[0, 1, x.t1, x.t2],
                    Glnq[a, x.t1, 0, x.t2],
                    Glnq[0, x.t1, x.t2, 0],
                    Glnq[0, x.t1, x.t2, 1],
                    Glnq[x.t1, 0, 0, x.t2]
                })
                .Where(m => Prod(m).Equals(J) && (!special || Glnq.Determinant(m).Equals(a.One)))
                .ToHashSet();
            
            if (q == 2)
                gens.Add(Glnq[1, 0, 1, 1]);

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

            if (group.Count() != ord)
                throw new($"############ actual {group.Count()} expected {ord}");
        
            Console.WriteLine($"q={q} q^2={q2} a={a} conj(a)={ax}");
            return group;
        }
        catch (Exception e)
        {
            Console.WriteLine(e);
        }
    }

    throw new();
}

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var infos = new[]
    {
        (2, 18, false), (2, 6, true),
        (4, 300, false), (4, 60, true),
        (8, 4536, false), (8, 504, true),
        (3, 96, false), (3, 24, true),
        (9, 7200, false), (9, 720, true),
        (5, 720, false), (5, 120, true),
        (7, 2688, false), (7, 336, true),
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
