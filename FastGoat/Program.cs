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
Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracketNoFmt;

ConcreteGroup<MatFq> GenGO3q(int q, bool special)
{
    var Glnq = new GLnq(3, q);
    var a = Glnq.Fq.X;
    var arrFq = Group.MulGroup($"F{q}", a);
    var pows = q.Range(-1).ToDictionary(k => k == -1 ? a.Zero : a.Pow(k), k => k);
    var S = Glnq.Neutral();
    var possibles = arrFq.Append(a.Zero)
        .Select(x0 => (x: x0, yList: arrFq.Append(a.Zero).Where(y0 => !x0.Equals(a.One) && (y0 * y0).Equals(1 - x0 * x0)).ToArray()))
        .Where(e => e.yList.Length != 0)
        .SelectMany(e => e.yList.Select(y0 => (e.x, y: y0)))
        .Distinct()
        .Select(e => (e.x, e.y, mat: Glnq[1, 0, 0, 0, e.x, e.y, 0, -e.y, e.x]))
        .Where(m => Glnq.Op(m.mat, m.mat.T).Equals(S) && (!special || m.mat.Det.Equals(a.One)))
        .Select(e => (e, Group.Cycle(Glnq, e.mat).Count))
        .OrderByDescending(e => e.Item2)
        .ThenByDescending(e => pows[e.e.x])
        .ToArray();
    
    var e0 = special
        ? a.One
        : arrFq.Where(e => e.Inv().Equals(e)).OrderByDescending(e => arrFq.ElementsOrders[e]).FirstOrDefault(a.One);
    var ide = Glnq[e0, 0, 0, 0, e0, 0, 0, 0, e0];
    
    var m0 = possibles[0].e.mat;
    var m1 = q != 5 ? Glnq[0, 1, 0, 0, 0, 1, 1, 0, 0] : Glnq[3, 1, 1, 1, 4, 3, 4, 2, 1];
    var gens = new List<MatFq>() { m0, Glnq.Op(ide, m1) };
    
    var name = special ? $"SO(3,{q})" : $"GO(3,{q})";
    var group = Group.Generate(name, Glnq, [..gens]);
    var nb = special ? FG.SOnqOrder(3, q) : FG.GOnqOrder(3, q);
    if (nb != group.Count())
    {
        gens.Println(m => $"[{m.Table.Glue(",").Replace(" ", "")}] ord:{group.ElementsOrders[m]}",
            $"Generators of {group.ShortName}");
        throw new($"actual:{group.ShortName} expected:{nb}");
    }
    
    return group;
}

void OrthogonalGroup3()
{
    foreach (var q in new[] { 2, 4, 8, 16, 32, 3, 9, 27, 5, 25, 7, 11, 13, 17, 19, 23, 29 })
    {
        var GO = GenGO3q(q, special: false);
        GO.GetGenerators().Println(m => $"[{m.Table.Glue(",").Replace(" ", "")}] ord:{GO.ElementsOrders[m]}",
            $"Generators of {GO.ShortName}");

        Console.WriteLine();
        
        var SO = GenGO3q(q, special: true);
        SO.GetGenerators().Println(m => $"[{m.Table.Glue(",").Replace(" ", "")}] ord:{SO.ElementsOrders[m]}",
            $"Generators of {SO.ShortName}");

        Console.WriteLine();
    }
}

{
    OrthogonalGroup3();
    Console.Beep();
}

/* Generators of |GO(3,2)| = 6
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,0,1,0,1,0] ord:2
   
   Generators of |SO(3,2)| = 6
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,0,1,0,1,0] ord:2
   
   Generators of |GO(3,4)| = 60
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,x+1,x,0,x,x+1] ord:2
   
   Generators of |SO(3,4)| = 60
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,x+1,x,0,x,x+1] ord:2
   
   Generators of |GO(3,8)| = 504
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,x²+1,x²,0,x²,x²+1] ord:2
   
   Generators of |SO(3,8)| = 504
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,x²+1,x²,0,x²,x²+1] ord:2
   
   Generators of |GO(3,16)| = 4080
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,x³+1,x³,0,x³,x³+1] ord:2
   
   Generators of |SO(3,16)| = 4080
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,x³+1,x³,0,x³,x³+1] ord:2
   
   Generators of |GO(3,32)| = 32736
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,x⁴+x,x⁴+x+1,0,x⁴+x+1,x⁴+x] ord:2
   
   Generators of |SO(3,32)| = 32736
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,x⁴+x,x⁴+x+1,0,x⁴+x+1,x⁴+x] ord:2
   
   Generators of |GO(3,3)| = 48
       [0,2,0,0,0,2,2,0,0] ord:6
       [1,0,0,0,0,1,0,2,0] ord:4
   
   Generators of |SO(3,3)| = 24
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,0,1,0,2,0] ord:4
   
   Generators of |GO(3,9)| = 1440
       [0,2,0,0,0,2,2,0,0] ord:6
       [1,0,0,0,2·x+2,x+1,0,2·x+2,2·x+2] ord:8
   
   Generators of |SO(3,9)| = 720
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,2·x+2,x+1,0,2·x+2,2·x+2] ord:8
   
   Generators of |GO(3,27)| = 39312
       [0,2,0,0,0,2,2,0,0] ord:6
       [1,0,0,0,2·x²+1,2·x²+x+2,0,x²+2·x+1,2·x²+1] ord:28
   
   Generators of |SO(3,27)| = 19656
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,2·x²+1,2·x²+x+2,0,x²+2·x+1,2·x²+1] ord:28
   
   Generators of |GO(3,5)| = 240
       [1,0,0,0,0,1,0,4,0] ord:4
       [2,4,4,4,1,2,1,3,4] ord:10
   
   Generators of |SO(3,5)| = 120
       [1,0,0,0,0,1,0,4,0] ord:4
       [3,1,1,1,4,3,4,2,1] ord:5
   
   Generators of |GO(3,25)| = 31200
       [0,4,0,0,0,4,4,0,0] ord:6
       [1,0,0,0,x+1,x+3,0,4·x+2,x+1] ord:24
   
   Generators of |SO(3,25)| = 15600
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,x+1,x+3,0,4·x+2,x+1] ord:24
   
   Generators of |GO(3,7)| = 672
       [0,6,0,0,0,6,6,0,0] ord:6
       [1,0,0,0,5,2,0,5,5] ord:8
   
   Generators of |SO(3,7)| = 336
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,5,2,0,5,5] ord:8
   
   Generators of |GO(3,11)| = 2640
       [0,10,0,0,0,10,10,0,0] ord:6
       [1,0,0,0,3,5,0,6,3] ord:12
   
   Generators of |SO(3,11)| = 1320
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,3,5,0,6,3] ord:12
   
   Generators of |GO(3,13)| = 4368
       [0,12,0,0,0,12,12,0,0] ord:6
       [1,0,0,0,11,6,0,7,11] ord:12
   
   Generators of |SO(3,13)| = 2184
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,11,6,0,7,11] ord:12
   
   Generators of |GO(3,17)| = 9792
       [0,16,0,0,0,16,16,0,0] ord:6
       [1,0,0,0,6,13,0,4,6] ord:16
   
   Generators of |SO(3,17)| = 4896
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,6,13,0,4,6] ord:16
   
   Generators of |GO(3,19)| = 13680
       [0,18,0,0,0,18,18,0,0] ord:6
       [1,0,0,0,3,7,0,12,3] ord:20
   
   Generators of |SO(3,19)| = 6840
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,3,7,0,12,3] ord:20
   
   Generators of |GO(3,23)| = 24288
       [0,22,0,0,0,22,22,0,0] ord:6
       [1,0,0,0,19,10,0,13,19] ord:24
   
   Generators of |SO(3,23)| = 12144
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,19,10,0,13,19] ord:24
   
   Generators of |GO(3,29)| = 48720
       [0,28,0,0,0,28,28,0,0] ord:6
       [1,0,0,0,5,18,0,11,5] ord:28
   
   Generators of |SO(3,29)| = 24360
       [0,1,0,0,0,1,1,0,0] ord:3
       [1,0,0,0,5,18,0,11,5] ord:28
*/