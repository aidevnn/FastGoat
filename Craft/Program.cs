using System.Text;
using Craft;
using Craft.Craft;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.Tools;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Perm.Style = DisplayPerm.CyclesComma;

Perm AutZnToPerm(Automorphism<ZnInt> aut)
{
    var su = new Sn(aut.Domain.Order);
    return su.CreateElementTable(aut.AutMap.Select(e => (e.Key.K, e.Value.K)).OrderBy(e => e.Item1)
        .Select(e => e.Item2).ToArray());
}

ConcreteGroup<Perm> UnPerm(Un u)
{
    var su = new Sn(u.Cn.Order);
    return Group.Generate($"{u.Name}", su, u.GetGenerators().Select(aut => AutZnToPerm(aut)).ToArray());
}

(Perm a, Perm b, Perm[] gensUi) AutDnPerm(int n)
{
    var Ui = IntExt.PrimesDec(n).Select(e => new Un(e.Key.Pow(e.Value))).OrderBy(u => u.Cn.Order).ToArray();
    var a = FG.Cycles(n);
    var b = FG.ConcatPerm(Ui.Select(u => AutZnToPerm(u[(u.Cn[1], -u.Cn[1])])).ToArray());
    var gpUi = Product.Gp(Ui.Select(u => UnPerm(u)).Cast<IGroup<Perm>>().ToArray());
    var gensUi = gpUi.GetGenerators().Select(e => FG.ConcatPerm(e.Ei)).ToArray();
    return (a, b, gensUi);
}

void RunAutomorphismDihedrals()
{
    for (int n = 3; n <= 128; n++)
    {
        Console.WriteLine($"################################## Aut[D{2 * n}]");
        var (a, b, gensUi) = AutDnPerm(n);
        var d2n = Group.Generate($"D{2 * n}", a.Sn, a, b);
        var autD2n2 = Group.Generate($"Aut[{d2n}]", a.Sn, gensUi.Prepend(a).ToArray());
        DisplayGroup.HeadOrdersGenerators(d2n);
        DisplayGroup.HeadOrdersGenerators(autD2n2); // Aut(D2n) ~ Hol(Cn)
        Console.WriteLine();
        
        if (!UGCraft.InnerAutCn(FG.Cycles(n)).ToHashSet().SetEquals(autD2n2))
            throw new();
    }
}

{
    RunAutomorphismDihedrals();
}