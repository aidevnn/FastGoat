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

ConcreteGroup<Perm> RegPermAutGroup(ConcreteGroup<Automorphism<Perm>> aut)
{
    var Dom = aut.Neutral().Domain;
    var sn = Dom.Neutral().Sn;
    if (sn.N != Dom.Count())
        throw new();

    var dic = Dom.Index().ToDictionary(e => e.Item, e => e.Index);
    var gens = aut.GetGenerators().Select(e => e.AutMap.ToDictionary(f => dic[f.Key], f => dic[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .ToArray();

    return Group.Generate(aut.Name, sn, gens);
}

ConcreteGroup<Perm> DummyPermForm<T>(AllSubgroups<T> gSubs) where T : struct, IElt<T>
{
    if (gSubs.Parent.GroupType == GroupType.AbelianGroup)
        return GroupPermutationForm.AbelianPerm(Group.AbelianGroupType(gSubs.Parent));

    var prods = gSubs.DecomposeProducts(gSubs.ProperNonTrivialNormalSubgroups())
        .DistinctBy(e => $"{e.lhs.Representative.Name} {e.rhs.Representative.Name}").ToList();
    if (prods.Any(e => e.isDirectProduct))
    {
        var (H, K, _) = prods.Where(e => e.isDirectProduct).MaxBy(e => e.lhs.Order);
        var h = DummyPermForm(gSubs.Restriction(H.Representative));
        var k = DummyPermForm(gSubs.Restriction(K.Representative));
        var nh = h.Neutral().Sn.N;
        var nk = k.Neutral().Sn.N;
        var hGens = h.GetGenerators().Select(e => UGCraft.PaddingRight(e, nk));
        var kGens = k.GetGenerators().Select(e => UGCraft.PaddingLeft(e, nh));
        var sn = new Sn(nh + nk);
        var HK = Group.Generate(gSubs.Parent.Name, sn, hGens.Concat(kGens).ToArray());
        if (!HK.IsIsomorphicTo(gSubs.Parent))
            throw new("###1");

        return HK;
    }
    else if (prods.Count == 0)
        return gSubs.Parent.ToPermGroup().Item1;
    else
    {
        var ordG = gSubs.Parent.Count();
        foreach (var (H, K, _) in prods.OrderBy(e => e.lhs.Order)
                     .ThenBy(e => e.lhs.GroupType == GroupType.NonAbelianGroup ? 0 : 1))
        {
            var k = DummyPermForm(gSubs.Restriction(K.Representative));
            var h = DummyPermForm(gSubs.Restriction(H.Representative)).ToPermGroup().Item1;
            var (snH, snK) = (h.Neutral().Sn, k.Neutral().Sn);
            var (nH, nK) = (snH.N, snK.N);

            var autHpg = RegPermAutGroup(Group.AutomorphismGroup(h)); // TODO: better automorphisms
            var homKAutH = Group.AllHomomorphisms(k, autHpg);
            foreach (var hom in homKAutH.OrderBy(e => e.Kernel().Count()))
            {
                // Console.WriteLine($"hom:{hom}");
                var img = hom.Image().ToArray();
                if (img.Length == k.Count() && k.IsIsomorphicTo(Group.Generate("Image", snH, img)))
                {
                    var HK0 = Group.Generate(gSubs.Parent.Name, snH, h.GetGenerators().Concat(img).ToArray());
                    if (HK0.IsIsomorphicTo(gSubs.Parent))
                        return HK0;
                }

                var hGens = h.GetGenerators().Select(f => UGCraft.PaddingRight(f, nK)).ToArray();
                var kGens = k.GetGenerators().Select(e => UGCraft.ConcatPerm(hom[e], e)).ToArray();
                var sn = new Sn(nH + nK);
                var hkGens = hGens.Concat(kGens).ToArray();
                // hkGens.Println();
                if (GroupCraft.GenerateElementsLimited(sn, hkGens, ordG).Count != ordG)
                    continue;

                var HK = Group.Generate(gSubs.Parent.Name, sn, hkGens);
                if (HK.IsIsomorphicTo(gSubs.Parent))
                    return HK;
            }

            Console.WriteLine(
                $"NOT FOUND  H:{H.Representative.ShortName} K:{K.Representative.ShortName} {autHpg.ShortName}");
            Console.WriteLine($"H in {h.Neutral().Sn} K in {k.Neutral().Sn}");
            throw new($"###2 {gSubs.Parent.ShortName}");
        }

        throw new($"###2 {gSubs.Parent.ShortName}");
    }
}

void MetaCyclicPermForm(int maxOrd)
{
    var (errors, total) = (0, 0);
    GlobalStopWatch.Restart();

    for (int ord = 1; ord <= maxOrd; ord++)
    {
        foreach (var (m, n, r) in GroupPermutationForm.MetaCyclicSdp(ord))
        {
            ++total;
            var Gpg = DummyPermForm(FG.MetaCyclicSdpWg(m, n, r).AllSubgroups());
            DisplayGroup.HeadElements(Gpg);
            var sn = Gpg.Neutral().Sn;
            var sn0 = GroupPermutationForm.MetaCyclicGens(m, n, r).sn;
            var test = sn.N == sn0.N;
            Console.WriteLine($"{Gpg.ShortName} in {sn} expected {sn0} Test:{test}");
            if (!test)
                ++errors;

            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show($"Errors:{errors}/{total} Max Order:{maxOrd}");
}

void AnyGroupPermFormUpTo(int maxOrd)
{
    var total = 0;
    GlobalStopWatch.Restart();

    for (int ord = 1; ord <= maxOrd; ord++)
    {
        foreach (var G in FG.AllGroupsOfOrder(ord))
        {
            ++total;
            var Gpg = DummyPermForm(G.AllSubgroups());
            DisplayGroup.HeadElements(Gpg);
            var sn = Gpg.Neutral().Sn;
            Console.WriteLine($"{G.ShortName} in {sn} Ratio:{G.Count() / (sn.N + 0.0):f3}");
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show($"Total:{total}  Max Order:{maxOrd}");
}

{
    // Optimal Form with smallest permutation order
    MetaCyclicPermForm(64); // # Errors:0/129 Max Order:64 Time:3.303s

    // SubOptimal Form
    AnyGroupPermFormUpTo(63); // # Total:319  Max Order:63 Time:11.523s
    
    // |(C3 x C3) x: S3| = 54
    // Type        NonAbelianGroup
    // BaseGroup   S9
    // 
    // Elements
    // ( 1)[1] = []
    // ( 2)[2] = [(4, 7), (5, 8), (6, 9)]
    // ( 3)[2] = [(4, 9), (5, 7), (6, 8)]
    // ( 4)[2] = [(4, 8), (5, 9), (6, 7)]
    // ( 5)[2] = [(1, 7), (2, 8), (3, 9)]
    // ( 6)[2] = [(1, 9), (2, 7), (3, 8)]
    // ( 7)[2] = [(1, 8), (2, 9), (3, 7)]
    // ( 8)[2] = [(1, 4), (2, 5), (3, 6)]
    // ( 9)[2] = [(1, 6), (2, 4), (3, 5)]
    // (10)[2] = [(1, 5), (2, 6), (3, 4)]
    // (11)[3] = [(4, 5, 6), (7, 9, 8)]
    // (12)[3] = [(4, 6, 5), (7, 8, 9)]
    // (13)[3] = [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
    // (14)[3] = [(1, 2, 3), (4, 6, 5)]
    // (15)[3] = [(1, 2, 3), (7, 9, 8)]
    // (16)[3] = [(1, 3, 2), (4, 6, 5), (7, 9, 8)]
    // (17)[3] = [(1, 3, 2), (7, 8, 9)]
    // (18)[3] = [(1, 3, 2), (4, 5, 6)]
    // (19)[3] = [(1, 4, 7), (2, 5, 8), (3, 6, 9)]
    // (20)[3] = [(1, 5, 7), (2, 6, 8), (3, 4, 9)]
    // (21)[3] = [(1, 6, 7), (2, 4, 8), (3, 5, 9)]
    // (22)[3] = [(1, 5, 9), (2, 6, 7), (3, 4, 8)]
    // (23)[3] = [(1, 6, 9), (2, 4, 7), (3, 5, 8)]
    // (24)[3] = [(1, 4, 9), (2, 5, 7), (3, 6, 8)]
    // (25)[3] = [(1, 6, 8), (2, 4, 9), (3, 5, 7)]
    // (26)[3] = [(1, 4, 8), (2, 5, 9), (3, 6, 7)]
    // (27)[3] = [(1, 5, 8), (2, 6, 9), (3, 4, 7)]
    // (28)[3] = [(1, 7, 4), (2, 8, 5), (3, 9, 6)]
    // (29)[3] = [(1, 9, 4), (2, 7, 5), (3, 8, 6)]
    // (30)[3] = [(1, 8, 4), (2, 9, 5), (3, 7, 6)]
    // (31)[3] = [(1, 8, 6), (2, 9, 4), (3, 7, 5)]
    // (32)[3] = [(1, 7, 6), (2, 8, 4), (3, 9, 5)]
    // (33)[3] = [(1, 9, 6), (2, 7, 4), (3, 8, 5)]
    // (34)[3] = [(1, 9, 5), (2, 7, 6), (3, 8, 4)]
    // (35)[3] = [(1, 8, 5), (2, 9, 6), (3, 7, 4)]
    // (36)[3] = [(1, 7, 5), (2, 8, 6), (3, 9, 4)]
    // (37)[6] = [(1, 2, 3), (4, 8, 6, 7, 5, 9)]
    // (38)[6] = [(1, 2, 3), (4, 7, 6, 9, 5, 8)]
    // (39)[6] = [(1, 2, 3), (4, 9, 6, 8, 5, 7)]
    // (40)[6] = [(1, 3, 2), (4, 9, 5, 7, 6, 8)]
    // (41)[6] = [(1, 3, 2), (4, 8, 5, 9, 6, 7)]
    // (42)[6] = [(1, 3, 2), (4, 7, 5, 8, 6, 9)]
    // (43)[6] = [(1, 9, 3, 8, 2, 7), (4, 5, 6)]
    // (44)[6] = [(1, 8, 2, 9, 3, 7), (4, 6, 5)]
    // (45)[6] = [(1, 8, 3, 7, 2, 9), (4, 5, 6)]
    // (46)[6] = [(1, 7, 2, 8, 3, 9), (4, 6, 5)]
    // (47)[6] = [(1, 9, 2, 7, 3, 8), (4, 6, 5)]
    // (48)[6] = [(1, 7, 3, 9, 2, 8), (4, 5, 6)]
    // (49)[6] = [(1, 5, 2, 6, 3, 4), (7, 9, 8)]
    // (50)[6] = [(1, 6, 3, 5, 2, 4), (7, 8, 9)]
    // (51)[6] = [(1, 5, 3, 4, 2, 6), (7, 8, 9)]
    // (52)[6] = [(1, 4, 2, 5, 3, 6), (7, 9, 8)]
    // (53)[6] = [(1, 6, 2, 4, 3, 5), (7, 9, 8)]
    // (54)[6] = [(1, 4, 3, 6, 2, 5), (7, 8, 9)]
    // 
    // |(C3 x C3) x: S3| = 54 in S9 Ratio:6.000
    // 
    // 
}