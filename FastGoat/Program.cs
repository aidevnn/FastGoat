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
    
    // conj(conj(x))=x then ax(ax) = a
    var ax = arrFq.First(e => !e.IsZero() && !e.Equals(a) && a.F.Substitute(e).IsZero() && e.Substitute(e).Equals(a));

    var J = Glnq[0, 1, 1, 0];
    MatFq Prod(MatFq m) => Glnq.Op(SelfAdjoint(m, ax), Glnq.Op(J, m));

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

    return group;
}

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    foreach (var q in new[] { 2, 4, 8, 16, 3, 9, 5, 7, 11, 13 })
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
/* All A in GU(2,2), bA*J*A=J True
   Generators of |GU(2,2)| = 18
       Glnq[0,1,1,0]
       Glnq[0,1,1,1]
       Glnq[0,x,x,0]
   
   # GU(2,2) Time:536ms
   
   All A in SU(2,2), bA*J*A=J True, Det A=1 True
   Generators of |SU(2,2)| = 6
       Glnq[0,1,1,0]
       Glnq[0,1,1,1]
   
   # SU(2,2) Time:8ms
   
   All A in GU(2,4), bA*J*A=J True
   Generators of |GU(2,4)| = 300
       Glnq[0,1,1,0]
       Glnq[0,1,1,1]
       Glnq[0,1,1,x^2 + x]
       Glnq[0,x,x^3 + x^2 + x,0]
   
   # GU(2,4) Time:289ms
   
   All A in SU(2,4), bA*J*A=J True, Det A=1 True
   Generators of |SU(2,4)| = 60
       Glnq[0,1,1,0]
       Glnq[0,1,1,1]
       Glnq[0,1,1,x^2 + x]
   
   # SU(2,4) Time:92ms
   
   All A in GU(2,8), bA*J*A=J True
   Generators of |GU(2,8)| = 4536
       Glnq[0,1,1,0]
       Glnq[0,1,1,1]
       Glnq[0,1,1,x^4 + x^2 + x]
       Glnq[0,x,x^5 + x^3 + x^2,0]
   
   # GU(2,8) Time:3.758s
   
   All A in SU(2,8), bA*J*A=J True, Det A=1 True
   Generators of |SU(2,8)| = 504
       Glnq[0,1,1,0]
       Glnq[0,1,1,1]
       Glnq[0,1,1,x^4 + x^2 + x]
   
   # SU(2,8) Time:1.664s
   
   All A in GU(2,16), bA*J*A=J True
   Generators of |GU(2,16)| = 69360
       Glnq[0,1,1,0]
       Glnq[0,1,1,1]
       Glnq[0,1,1,x^3 + x]
       Glnq[0,x,x^4 + x^2 + x,0]
   
   # GU(2,16) Time:3m40s
   
   All A in SU(2,16), bA*J*A=J True, Det A=1 True
   Generators of |SU(2,16)| = 4080
       Glnq[0,1,1,0]
       Glnq[0,1,1,1]
       Glnq[0,1,1,x^3 + x]
   
   # SU(2,16) Time:38.470s
   
   All A in GU(2,3), bA*J*A=J True
   Generators of |GU(2,3)| = 96
       Glnq[0,1,1,0]
       Glnq[0,1,1,x + 1]
       Glnq[0,x,2*x,0]
   
   # GU(2,3) Time:15ms
   
   All A in SU(2,3), bA*J*A=J True, Det A=1 True
   Generators of |SU(2,3)| = 24
       Glnq[0,x + 1,x + 1,0]
       Glnq[0,x + 1,x + 1,1]
   
   # SU(2,3) Time:11ms
   
   All A in GU(2,9), bA*J*A=J True
   Generators of |GU(2,9)| = 7200
       Glnq[0,1,1,0]
       Glnq[0,1,1,x^2 + 2]
       Glnq[0,x,x^3 + 2,0]
   
   # GU(2,9) Time:3.934s
   
   All A in SU(2,9), bA*J*A=J True, Det A=1 True
   Generators of |SU(2,9)| = 720
       Glnq[0,x^2 + 2,2*x^3 + x^2 + 2*x + 1,0]
       Glnq[0,x^2 + 2,2*x^3 + x^2 + 2*x + 1,1]
       Glnq[0,x^3 + x + 1,2*x^3 + 2*x^2 + 2*x,0]
       Glnq[0,x^3 + x + 1,2*x^3 + 2*x^2 + 2*x,1]
   
   # SU(2,9) Time:1.336s
   
   All A in GU(2,5), bA*J*A=J True
   Generators of |GU(2,5)| = 720
       Glnq[0,1,1,0]
       Glnq[0,1,1,x + 2]
       Glnq[0,x,3*x,0]
   
   # GU(2,5) Time:140ms
   
   All A in SU(2,5), bA*J*A=J True, Det A=1 True
   Generators of |SU(2,5)| = 120
       Glnq[0,x + 2,2*x + 4,0]
       Glnq[0,x + 2,2*x + 4,1]
   
   # SU(2,5) Time:70ms
   
   All A in GU(2,7), bA*J*A=J True
   Generators of |GU(2,7)| = 2688
       Glnq[0,1,1,0]
       Glnq[0,1,1,x + 3]
       Glnq[0,x,5*x,0]
   
   # GU(2,7) Time:808ms
   
   All A in SU(2,7), bA*J*A=J True, Det A=1 True
   Generators of |SU(2,7)| = 336
       Glnq[0,x + 3,x + 3,0]
       Glnq[0,x + 3,x + 3,1]
   
   # SU(2,7) Time:275ms
   
   All A in GU(2,11), bA*J*A=J True
   Generators of |GU(2,11)| = 15840
       Glnq[0,1,1,0]
       Glnq[0,1,1,x +  9]
       Glnq[0,x, 6*x,0]
   
   # GU(2,11) Time:17.813s
   
   All A in SU(2,11), bA*J*A=J True, Det A=1 True
   Generators of |SU(2,11)| = 1320
       Glnq[0,x +  9, 5*x + 1,0]
       Glnq[0,x +  9, 5*x + 1,1]
   
   # SU(2,11) Time:1.978s
   
   All A in GU(2,13), bA*J*A=J True
   Generators of |GU(2,13)| = 30576
       Glnq[0,1,1,0]
       Glnq[0,1,1,x +  6]
       Glnq[0,x, 7*x,0]
   
   # GU(2,13) Time:1m6s
   
   All A in SU(2,13), bA*J*A=J True, Det A=1 True
   Generators of |SU(2,13)| = 2184
       Glnq[0,x +  6, 8*x +  9,0]
       Glnq[0,x +  6, 8*x +  9,1]
   
   # SU(2,13) Time:3.483s
*/