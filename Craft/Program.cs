using System.Numerics;
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

(Perm a, Perm b, Perm[] gensUi) AutD2nPerm(int n)
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
        var (a, b, gensUi) = AutD2nPerm(n);
        var d2n = Group.Generate($"D{2 * n}", a.Sn, a, b);
        var autD2n2 = Group.Generate($"Aut[{d2n}]", a.Sn, gensUi.Prepend(a).ToArray());
        DisplayGroup.HeadOrdersGenerators(d2n);
        DisplayGroup.HeadOrdersGenerators(autD2n2); // Aut(D2n) ~ Hol(Cn)
        Console.WriteLine();

        if (!UGCraft.HolCn(FG.Cycles(n)).ToHashSet().SetEquals(autD2n2))
            throw new();
    }
}

void TestHolCn()
{
    for (int i = 1; i <= 128; i++)
    {
        var a = FG.Cycles(i);
        var holCaOrd = i * IntExt.Phi(i);
        var holCa1 = UGCraft.HolCn(a).ToArray();
        if (i > 16)
        {
            Console.WriteLine($"|Hol(C{i})| = {holCa1.Length}/{holCaOrd}");
            continue;
        }

        var holCa2 = UGCraft.HolCnMatrix(a).ToArray();
        Console.WriteLine($"|Hol(C{i})| = {holCa1.Length}/{holCa2.Length}/{holCaOrd}");
    }

    Console.WriteLine();
}

void SearchHolDn(int n)
{
    var d2nOrd = 2 * n;
    var autD2nOrd = n * IntExt.Phi(n);
    var holOrd = d2nOrd * autD2nOrd;
    var snOrder = n.SeqLazy(1).Aggregate(BigInteger.One, (acc, e) => e * acc);
    Console.WriteLine(
        $"{$"Hol[D{2 * n}]",-10} in  {new Sn(n),-4}:{(snOrder % holOrd == 0 ? "Unknown" : "Impossible")}");
}

(Perm a, Perm b) PrepareDihedral(int n)
{
    var pType = IntExt.PrimesDec(n).Select(e => e.Key.Pow(e.Value)).ToArray();
    var a0 = IntExt.PermAndCyclesFromType(pType);
    var n1 = a0.perm.Length;
    var sn = new Sn(2 * n1);
    var bCycles = a0.cycles.Zip(a0.cycles)
        .SelectMany(e => e.First.Zip(e.Second.Select(f => f + n1).Reverse().ToArray()))
        .ToArray();
    var a = FG.ConcatPerm(FG.Cycles(n), FG.Cycles(n));
    var b = sn.OpSeq(bCycles.Select(c => sn.CycleP1([c.First, c.Second])));
    return (a, b);
}

Dictionary<T2, HashSet<T2>> AllOrbits<T1, T2>(ConcreteGroup<T1> gr, T2[] set, GroupAction<T1, T2> act)
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
    var xsets = new HashSet<T2>();
    var allOrbits = new Dictionary<T2, HashSet<T2>>();
    var gens = gr.GetGenerators().ToHashSet();
    foreach (var x in set)
    {
        if (xsets.Contains(x))
            continue;

        var orbx = Group.Orbits(gens, act, x);
        if (!xsets.Overlaps(orbx))
            allOrbits[orbx.Min()] = orbx;

        xsets.UnionWith(orbx);
    }

    return allOrbits;
}

ConcreteGroup<Perm> ProductPermGroup(params ConcreteGroup<Perm>[] Gs)
{
    var HnBase = Product.Gp(Gs.Cast<IGroup<Perm>>().ToArray());
    var HnGens = HnBase.GetGenerators().Select(e => FG.ConcatPerm(e.Ei)).ToArray();
    return Group.Generate(HnBase.Name, HnGens[0].Sn, HnGens);
}

(ConcreteGroup<Perm> d2n, ConcreteGroup<Perm> autd2n, ConcreteGroup<Perm> hold2n) HolomorphD2n(int n)
{
    Console.WriteLine($"#### Start D{2 * n,-6} ####");
    if (IntExt.PrimesDec(n).Count == 1)
    {
        var _d2n = FG.Dihedral(n);
        var (_d2nReg, _autD2nReg, _) = FG.RegPermAutGroup(_d2n);
        var _holD2n = Group.DirectProduct($"Hol[{_d2n}]", _d2nReg, _autD2nReg);
        DisplayGroup.HeadOrdersGenerators(_d2n);
        DisplayGroup.HeadOrdersGenerators(_autD2nReg);
        DisplayGroup.HeadOrdersGenerators(_holD2n);
        return (_d2n, _autD2nReg, _holD2n);
    }

    var phi = IntExt.Phi(n);
    var phiType = Group.AbelianGroupType(FG.UnInt(n));
    var (a, b) = PrepareDihedral(n);
    Console.WriteLine(new { a });
    Console.WriteLine(new { b });

    var d2n = Group.Generate($"D{2 * n}", a.Sn, a, b);
    DisplayGroup.HeadOrdersGenerators(d2n);
    var autD2nReg = FG.RegPermAutGroup(d2n).autG;
    Console.WriteLine($"Aut[{d2n}] = C{n} x: {phiType.ToAbString().WithParenthesis()}");
    var holCa = UGCraft.HolCn(a).ToArray();
    var hol = holCa.GroupBy(e => e.Order).ToDictionary(e => e.Key, e => e.ToArray());
    var d2nSet = d2n.ToSet();
    var act = Group.ByConjugate(a.Sn);
    var actSet = Group.ByConjugateSet(a.Sn);
    hol.ToDictionary(e => e.Key, e => e.Value.Length).AscendingByKey()
        .Println($"Hol:{hol.Sum(e => e.Value.Length)}/{2 * n * n * phi}");

    var gHolCa = Group.Generate("InnAuts", a.Sn, holCa.OrderByDescending(e => e.Order).ThenBy(e => e.Orbits.Length).ToArray());
    var allOrbx = AllOrbits(gHolCa, holCa, act);
    var holOrd = 2 * n * n * phi;

    var CnOrd = allOrbx.Keys.Where(e => !d2n.Contains(e) && e.Order == n && act.IsInvariant(e, a))
        .OrderByDescending(e => e.Orbits.Length).ThenBy(e => e).ToArray();
    var UnOrd = phiType.Select(t =>
            allOrbx.Keys.Where(e => !d2n.Contains(e) && e.Order == t && actSet.IsInvariant(e, d2nSet))
                .OrderByDescending(e => e.Orbits.Length).ThenBy(e => e).ToArray())
        .ToArray();
    Console.WriteLine($"CnOrd:{CnOrd.Length} UnOrd:[{UnOrd.Select(e => e.Length).Glue(", ")}]");

    foreach (var g0 in CnOrd)
    {
        Console.WriteLine($"g0:{g0}");
        var G0 = Group.Generate("G", a.Sn, g0).ToSet();
        var filterUnOrd = UnOrd
            .MultiLoop().Select(l => l.ToArray())
            .Where(l => l.Distinct().Count() == l.Length &&
                        l.Grid2D().All(g => act.IsInvariant(g.t1, g.t2)) &&
                        l.All(e => actSet.IsInvariant(e, G0)))
            .ToArray();
        var gsAll = filterUnOrd.Where(l => l.All(e => actSet.IsInvariant(e, G0))).ToArray();
        Console.WriteLine($"   filterUnOrd:{filterUnOrd.Length} gs:{gsAll.Length}");
        var gs = gsAll.FirstOrDefault(l =>
        {
            var elts = GroupCraft.GenerateElementsLimited(a.Sn, l.Prepend(g0).ToArray(), n * phi);
            return elts.Count == n * phi &&
                   elts.Intersect(d2n).Count() == 1 &&
                   Group.Generate(a.Sn, l.Prepend(g0).ToArray()).IsIsomorphicTo(autD2nReg) &&
                   GroupCraft.GenerateElementsLimited(a.Sn, l.Concat([a, b, g0]).ToArray(), holOrd).Count == holOrd;
        }, [a.Sn.Neutral()]);

        if (gs.All(e => e.Order != 1))
        {
            var autD2n = Group.Generate($"Aut[{d2n}]", a.Sn, gs.Prepend(g0).ToArray());
            var holD2n = Group.Generate($"Hol[{d2n}]", a.Sn, gs.Concat([a, b, g0]).ToArray());
            DisplayGroup.HeadOrdersGenerators(autD2n);
            DisplayGroup.HeadOrdersGenerators(holD2n);
            return (d2n, autD2n, holD2n);
        }
    }

    throw new("bad");
}

void HolD2n(int n)
{
    {
        var (d2n, autD2n, holD2n) = HolomorphD2n(n);
        var sn = holD2n.Neutral().Sn;
        var (d2nReg, autD2nReg, _) = FG.RegPermAutGroup(d2n);
        // var holReg = Group.DirectProduct(d2nReg, autD2nReg);
        var snReg = d2nReg.Neutral().Sn;
        var phi = IntExt.Phi(n);
        var holOrd = 2 * n * n * phi;
        Console.WriteLine(
            $"{d2n.ShortName,-20} {autD2n.ShortName,-20} autD2nOrd:{n * phi,-10} {holD2n.ShortName,-20} holOrd:{holOrd,-10} {sn} {snReg}");
        if (!autD2n.IsIsomorphicTo(autD2nReg))
            throw new();
    }

    Console.WriteLine();
}

{
    GlobalStopWatch.Restart();
    for (int n = 3; n <= 32; n++)
        HolD2n(n);
    GlobalStopWatch.Show(); // Time:35.480s
    // |D60| = 60           |Aut[D60]| = 240     autD2nOrd:240        |Hol[D60]| = 14400   holOrd:14400      S20 S60
}