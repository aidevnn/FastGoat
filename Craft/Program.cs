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
using FastGoat.UserGroup.GModuleN;
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

GroupSubset<T> AddElt<T>(ConcreteGroup<T> G, GroupSubset<T> set, T e) where T : struct, IElt<T>
{
    var tmp = Group.UniqueGenerators(G, set.Generators.Append(e).ToArray());
    return new(tmp.uniqueGenerators.ToHashSet(), tmp.elements);
}

ConcreteGroup<T>[] AnyFactors<T>(ConcreteGroup<T> G, ConcreteGroup<T> H, bool details = false) where T : struct, IElt<T>
{
    var oN = G.Order / H.Order;
    var hNorm = Group.IsNormalSubgroup(G, H);
    var _cosets = Group.Cosets(G, H, CosetType.Left).Values.Distinct()
        .Select(c => c.Where(x => oN % G.ElementsOrders[x] == 0).ToXSet()).ToHashSet();
    var _sub = new GroupSubset<T>([G.Neutral()], [G.Neutral()]);
    List<(HashSet<XSet<T>>, GroupSubset<T>)> all = new() { (_cosets, _sub) };

    do
    {
        var tmp = all.ToArray();
        all.Clear();
        foreach (var (cosets, sub) in tmp)
        {
            if (cosets.Count == 0)
            {
                if (all.All(c => !c.Item2.Equals(sub)))
                    all.Add((cosets, sub));

                continue;
            }

            foreach (var set in cosets.ToArray())
            {
                if (H.SuperSetOf(set))
                    continue;

                foreach (var a in set.Except(sub).OrderBy(a => G.ElementsOrders[a]))
                {
                    var cosets0 = cosets.Select(c => c.ToXSet()).ToHashSet();
                    var sub0 = AddElt(G, sub, a);
                    if (all.Any(c => c.Item2.Equals(sub0)))
                        continue;

                    foreach (var e in sub0)
                    {
                        var s = H.Select(x => G.Op(e, x)).ToHashSet();
                        cosets0.RemoveWhere(c => s.Overlaps(c));
                    }

                    if (oN % sub0.Count == 0 && sub0.Intersect(H).Count() == 1)
                    {
                        if (cosets0.Count != 0
                            || (sub0.Count == oN 
                                && (hNorm || !hNorm && Group.IsNormalSubgroup(G, Group.Generate(G, sub0.Generators.ToArray())))))
                            all.Add((cosets0, sub0));
                    }
                }
            }
        }

        if (details)
        {
            Console.WriteLine("####################################################################");
            foreach (var (cosets1, sub1) in all)
                cosets1.Println($"gens:{sub1.Generators.ToXSet()} set:{sub1.ToXSet()}");

            Console.WriteLine();
        }
    } while (all.Any(sc => sc.Item1.Count != 0));

    return all.Select((sc, k) => Group.Generate($"K{k + 1}", G, sc.Item2.Generators.ToArray())).ToArray();
}

void FirstExamples()
{
    {
        var G = FG.Symmetric(4);
        var N = Group.Generate("C2 x C2", G, G[(1, 2), (3, 4)], G[(1, 3), (2, 4)]);
        var H = Group.FactorGroup(G, N);
        H.Name = "S3";
        
        DisplayGroup.HeadOrdersGenerators(G);
        DisplayGroup.HeadOrdersGenerators(H);
        GroupCraft.DisplayIsNormal(G, H);
        foreach (var Ni in AnyFactors(G, H))
            DisplayGroup.HeadOrdersGenerators(Ni);

        DisplayGroup.HeadOrdersGenerators(N);
        GroupCraft.DisplayIsNormal(G, N);
        foreach (var Ki in AnyFactors(G, N))
            DisplayGroup.HeadOrdersGenerators(Ki);
    }

    {
        var G = FG.Symmetric(4);
        var N = Group.Generate("A4", G, G[(1, 2, 3)], G[(1, 3, 4)]);
        var H = Group.FactorGroup(G, N);
        H.Name = "C2";
        
        DisplayGroup.HeadOrdersGenerators(G);
        DisplayGroup.HeadOrdersGenerators(H);
        GroupCraft.DisplayIsNormal(G, H);
        foreach (var Ni in AnyFactors(G, H))
            DisplayGroup.HeadOrdersGenerators(Ni);

        DisplayGroup.HeadOrdersGenerators(N);
        GroupCraft.DisplayIsNormal(G, N);
        foreach (var Ki in AnyFactors(G, N))
            DisplayGroup.HeadOrdersGenerators(Ki);
    }

    {
        var G = FG.Symmetric(4);
        var H = Group.Generate("D8", G, G[(1, 2), (3, 4)], G[(1, 2, 3, 4)]); // not factor
        
        DisplayGroup.HeadOrdersGenerators(G);
        DisplayGroup.HeadOrdersGenerators(H);
        GroupCraft.DisplayIsNormal(G, H);
        foreach (var Ni in AnyFactors(G, H))
            DisplayGroup.HeadOrdersGenerators(Ni);
    }

    {
        var G = FG.MetaCyclicPg(5, 4, 2);
        var N = Group.Generate("C5", G, G[(1, 2, 3, 4, 5)]);
        var H = Group.FactorGroup(G, N);
        H.Name = "C4";
        
        DisplayGroup.HeadOrdersGenerators(G);
        DisplayGroup.HeadOrdersGenerators(H);
        GroupCraft.DisplayIsNormal(G, H);
        foreach (var Ni in AnyFactors(G, H))
            DisplayGroup.HeadOrdersGenerators(Ni);

        DisplayGroup.HeadOrdersGenerators(N);
        GroupCraft.DisplayIsNormal(G, N);
        foreach (var Ki in AnyFactors(G, N))
            DisplayGroup.HeadOrdersGenerators(Ki);
    }

    {
        var s8 = new Sn(8);
        var G = Group.Generate("SL(2,3)", s8, s8[(2, 5, 8), (4, 7, 6)], s8[(1, 2, 3, 4), (5, 8, 7, 6)]);
        var H = Group.Generate("C3", s8, s8[(2, 5, 8), (4, 7, 6)]); // factor
        
        DisplayGroup.HeadOrdersGenerators(G);
        DisplayGroup.HeadOrdersGenerators(H);
        GroupCraft.DisplayIsNormal(G, H);
        foreach (var Ni in AnyFactors(G, H))
            DisplayGroup.HeadOrdersGenerators(Ni);
    }

    {
        var s11 = new Sn(11);
        var G = Group.Generate("F(3x:8)2", s11, s11[(1, 2, 3)], s11[(2, 3), (4, 5, 6, 7, 8, 9, 10, 11)]);
        var H = Group.Generate("C4", s11, s11[(4, 6, 8, 10), (5, 7, 9, 11)]); // non split
        
        DisplayGroup.HeadOrdersGenerators(G);
        DisplayGroup.HeadOrdersGenerators(H);
        GroupCraft.DisplayIsNormal(G, H);
        foreach (var Ni in AnyFactors(G, H))
            DisplayGroup.HeadOrdersGenerators(Ni);
    }

    {
        var s11 = new Sn(11);
        var G = Group.Generate("F(3x:8)2", s11, s11[(1, 2, 3)], s11[(2, 3), (4, 5, 6, 7, 8, 9, 10, 11)]);
        var N = Group.Generate("C6", s11, s11[(1, 2, 3), (4, 8), (5, 9), (6, 10), (7, 11)]); // non split
        
        DisplayGroup.HeadOrdersGenerators(G);
        DisplayGroup.HeadOrdersGenerators(N);
        GroupCraft.DisplayIsNormal(G, N);
        foreach (var Hi in AnyFactors(G, N))
            DisplayGroup.HeadOrdersGenerators(Hi);
    }
}

void FactorsUpToOrder32()
{
    for(int ord = 1; ord <= 32; ++ord)
    {
        foreach (var G in FG.AllGroupsOfOrder(ord).Where(g => g.GroupType == GroupType.NonAbelianGroup))
        {
            Console.WriteLine($"############################################## {G.ShortName}");
            var allSubs = G.AllSubgroups();
            allSubs.Naming();
            var norms = allSubs.ProperNonTrivialNormalSubgroups();
            var allNormals = allSubs.DecomposeProducts(norms)
                .SelectMany(e => e.lhs.Conjugates)
                .ToHashSet(new GroupSetEquality<Word>());
            var allFactors = allSubs.DecomposeProducts(norms)
                .SelectMany(e => e.rhs.Conjugates)
                .ToHashSet(new GroupSetEquality<Word>());
            DisplayGroup.HeadNames(G);
            foreach (var H in allSubs.Select(n => n.Representative))
            {
                if (H.Order == 1 || H.Order == G.Order)
                    continue;

                Console.WriteLine("---------------------------------------");
                DisplayGroup.HeadOrdersGenerators(H);
                var Ns = AnyFactors(G, H);
                foreach (var N in Ns)
                    DisplayGroup.HeadOrdersGenerators(N);

                if (Ns.Length == 0)
                {
                    if (Group.IsNormalSubgroup(G, H))
                    {
                        Console.WriteLine("NonSplit");
                        if (allNormals.Contains(H))
                            throw new();
                    }
                    else
                    {
                        Console.WriteLine("NotFactor");
                        if (allFactors.Contains(H))
                            throw new();
                    }
                }
                else
                    Ns.Select(n => allSubs.First(s => s.Contains(n)))
                        .Distinct()
                        .Select(n => GroupCraft.DirectProduct(H, n.Representative)).Println($"{G} = ");

                Console.WriteLine();
            }
        }
    }
}

{
    FirstExamples();
    // FactorsUpToOrder32();
}