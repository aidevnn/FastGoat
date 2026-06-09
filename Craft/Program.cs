using System.Numerics;
using System.Reflection;
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
Perm[] QuaternionGens(int n)
{
    if (!int.IsPow2(n))
        throw new();

    var k1 = n / 2;
    var k2 = n / 4;
    var sn = new Sn(n);
    var a = sn.Cycle(k1.Range(1)) * sn.Cycle(k1.Range(k1 + 1));
    var b = Group.OpSeq(sn, k2.SeqLazy().Select(i => sn.Cycle(i + 1, n - i, i + 1 + k2, n - i - k2)));
    return [a, b];
}
(Perm a, Perm b) OpsCnByUm(int n, Perm a, Perm b)
{
    var n0 = IntExt.DividorsInt(n).Order().First(x => IntExt.Lcm(x, b.Order) == n);
    var b0 = FG.Cycles(n0);
    return (FG.PaddingRight(a, b0.Sn.N), FG.ConcatPerm(b, b0));
}

int total, errors;

void HolDic(int n)
{
    var dicn = FG.DicyclicPg(n);
    var sn = dicn.Neutral().Sn;
    var autDicn = Group.AutomorphismGroup(dicn);
    DisplayGroup.HeadOrders(dicn);
    DisplayGroup.Generators(dicn);
    DisplayGroup.HeadOrders(autDicn);
    var (reg, autReg, _) = FG.RegPermAutGroup(dicn);
    var hol0 = Group.DirectProduct($"Hol[{dicn}]", reg, autReg);
    if (sn.N == dicn.Order)
    {
        DisplayGroup.HeadOrders(autReg);
        Console.WriteLine($"Hol[{dicn}] from Regular Permutation");
        DisplayGroup.HeadOrders(hol0);
    }
    else
    {
        var typeMatching = autDicn.GetGenerators()
            .All(aut => dicn.GetGenerators().All(g => Perm.TypeEquals(g, aut[g])));
        Console.WriteLine($"TypeMatch:{typeMatching}");
        var autOrd = autDicn.Order;
        var holOrd = dicn.Order * autOrd;

        if (n % 2 == 1)
        {
            var (b1, c1) = PrepareDihedral(n);
            var candidat = false;
            var (b2, c2) = OpsCnByUm(4, b1, c1);
            var dicn2 = Group.Generate(dicn.Name, b2.Sn, b2, c2);
            var autDicn2 = Group.AutomorphismGroup(dicn2);
            var aut1 = autDicn2.OrderByDescending(e => autDicn2.ElementsOrders[e])
                .First(aut => aut[b2].Equals(b2) && !aut[c2].Equals(c2));
            var aut2 = autDicn2.OrderByDescending(e => autDicn2.ElementsOrders[e])
                .First(aut => aut[b2].Equals(b2) && aut[c2].Equals(c2));
            var gensAut = GroupCraft.RecreateGenerators(autDicn2,
                autDicn2.OrderByDescending(e =>
                    e.Equals(aut1) ? 1000 :
                    e.Equals(aut2) ? 500 : autDicn2.ElementsOrders[e]).ToArray());
            var holGens = gensAut.ToDictionary(aut => aut, aut => UGCraft.InnerAutMatrix(aut).ToList());
            foreach (var gens in holGens.Values.MultiLoop().Select(l => l.ToArray()))
            {
                var eltsAut = GroupCraft.GenerateElementsLimited(b2.Sn, gens, autOrd);
                var listOrdAut = eltsAut.Select(e => e.Order).Order().ToArray();
                if (eltsAut.Count != autOrd || !listOrdAut.SequenceEqual(autDicn.ElementsOrdersList()))
                    continue;

                var gensF = gens.Concat(dicn2.GetGenerators()).ToArray();
                var eltsHol = GroupCraft.GenerateElementsLimited(b2.Sn, gensF, holOrd);
                var listOrdHol = eltsHol.Select(e => e.Order).Order().ToArray();
                if (eltsHol.Count == holOrd && listOrdHol.SequenceEqual(hol0.ElementsOrdersList()))
                {
                    var autDicn3 = Group.Generate($"Aut[{dicn}]", b2.Sn, gens);
                    DisplayGroup.HeadOrdersGenerators(autDicn3);
                    Console.WriteLine(
                        $"Idx gen:[{holGens.Values.Zip(gens).Select(e => e.First.FindIndex(f => f.Equals(e.Second))).Glue(", ")}]");
                    DisplayGroup.AreIsomorphics(autDicn3, autReg);
                    Console.WriteLine();
                    var hol = Group.Generate($"Hol[{dicn}](candidat)", b2.Sn, gensF);
                    DisplayGroup.HeadOrdersGenerators(hol);
                    candidat = true;
                    break;
                }
            }

            if (!candidat)
            {
                ++errors;
                Console.WriteLine("Problem Even Dicyclic group");
            }
        }
        else
        {
            var k = IntExt.PrimesDec(n)[2];
            var m = n / 2.Pow(k);
            var gensQ = QuaternionGens(1 << (k + 2));
            var (b1, c1) = gensQ.Deconstruct();
            var (a2, b2) = PrepareDihedral(m);
            var a3 = FG.PaddingRight(a2, b1.Sn.N);
            var b3 = FG.PaddingLeft(b1, b2.Sn.N);
            var c3 = FG.ConcatPerm(b2, c1);
            var dicn2 = Group.Generate(dicn.Name, a3.Sn, a3, b3, c3);
            var autDicn2 = Group.AutomorphismGroup(dicn2);
            var aut1 = autDicn2.OrderByDescending(e => autDicn2.ElementsOrders[e])
                .First(aut => aut[a3].Equals(a3) && !aut[b3].Equals(b3) && !aut[c3].Equals(c3));
            var aut2 = autDicn2.OrderByDescending(e => autDicn2.ElementsOrders[e])
                .First(aut => aut[a3].Equals(a3) && aut[b3].Equals(b3) && !aut[c3].Equals(c3));
            var aut3 = autDicn2.OrderByDescending(e => autDicn2.ElementsOrders[e])
                .First(aut => aut[a3].Equals(a3) && aut[b3].Equals(b3) && aut[c3].Equals(c3));
            var gensAut = GroupCraft.RecreateGenerators(autDicn2,
                autDicn2.OrderByDescending(e =>
                    e.Equals(aut1) ? 1000 :
                    e.Equals(aut2) ? 750 :
                    e.Equals(aut3) ? 500 : autDicn2.ElementsOrders[e]).ToArray());
            var holGens = gensAut.ToDictionary(aut => aut, aut => UGCraft.InnerAutMatrix(aut).ToList());
            var candidat = false;

            foreach (var gens in holGens.Values.MultiLoop().Select(l => l.ToArray()))
            {
                var eltsAut = GroupCraft.GenerateElementsLimited(a3.Sn, gens, autOrd);
                var listOrdAut = eltsAut.Select(e => e.Order).Order().ToArray();
                if (eltsAut.Count != autOrd || !listOrdAut.SequenceEqual(autDicn.ElementsOrdersList()))
                    continue;

                var gensF = gens.Concat(dicn2.GetGenerators()).ToArray();
                var elts = GroupCraft.GenerateElementsLimited(a3.Sn, gensF, holOrd);
                var listOrd = elts.Select(e => e.Order).Order().ToArray();
                if (elts.Count == holOrd && listOrd.SequenceEqual(hol0.ElementsOrdersList()))
                {
                    var autDicn3 = Group.Generate($"Aut[{dicn}]", a3.Sn, gens);
                    DisplayGroup.HeadOrdersGenerators(autDicn3);
                    Console.WriteLine(
                        $"Idx gen:[{holGens.Values.Zip(gens).Select(e => e.First.FindIndex(f => f.Equals(e.Second))).Glue(", ")}]");
                    DisplayGroup.AreIsomorphics(autDicn3, autReg);
                    Console.WriteLine();
                    var hol = Group.Generate($"Hol[{dicn}](candidat)", a3.Sn, gensF);
                    DisplayGroup.HeadOrdersGenerators(hol);
                    candidat = true;
                    break;
                }
            }

            if (!candidat)
            {
                ++errors;
                Console.WriteLine("Problem Even Dicyclic group");
            }
        }
    }

    DisplayGroup.HeadOrders(hol0);
    Console.WriteLine();
}

{
    GlobalStopWatch.Restart();
    (total, errors) = (0, 0);
    for (int n = 2; n <= 16; n++)
    {
        ++total;
        HolDic(n);
    }

    GlobalStopWatch.Show($"total:{total} errors:{errors}"); // # total:15 errors:0 Time:1m51s
    // |Aut[Dic7]| = 84
    // Type        NonAbelianGroup
    // BaseGroup   S18
    // 
    // Elements Orders : [1]:1, [2]:15, [3]:14, [6]:42, [7]:6, [14]:6
    // Generators of Aut[Dic7]
    // gen1 of order 6
    // [(1, 2, 7, 4, 3, 5), (8, 14, 9, 12, 13, 11)]
    // gen2 of order 14
    // [(1, 7, 6, 5, 4, 3, 2), (15, 16), (17, 18)]
    // 
    // 
}