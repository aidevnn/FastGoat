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
                     .ThenBy(e => e.lhs.Representative.GroupType == GroupType.NonAbelianGroup ? 0 : 1))
        {
            var k = DummyPermForm(gSubs.Restriction(K.Representative));
            var h = DummyPermForm(gSubs.Restriction(H.Representative)).ToPermGroup().Item1; 
            var (snH, snK) = (h.Neutral().Sn, k.Neutral().Sn);
            var (nH, nK) = (snH.N, snK.N);

            var autHpg = RegPermAutGroup(Group.AutomorphismGroup(h)); // TODO: better automorphisms
            var homKAutH = Group.AllHomomorphisms(k, autHpg);
            foreach (var hom in homKAutH)
            {
                // Console.WriteLine($"hom:{hom}");
                var hkGens = h.GetGenerators().Concat(hom.Image()).ToArray();
                var _HK = GroupCraft.GenerateElementsLimited(snH, hkGens, ordG);
                if (_HK.Count == ordG)
                {
                    var HK0 = Group.Generate(gSubs.Parent.Name, snH, hkGens);
                    if (HK0.IsIsomorphicTo(gSubs.Parent))
                        return HK0;
                }

                var hGens = h.GetGenerators().Select(f => UGCraft.PaddingRight(f, nK)).ToArray();
                var kGens = k.GetGenerators().Select(e => UGCraft.ConcatPerm(hom[e], e)).ToArray();
                var sn = new Sn(nH + nK);
                hkGens = hGens.Concat(kGens).ToArray();
                // hkGens.Println();
                _HK = GroupCraft.GenerateElementsLimited(sn, hkGens, ordG);
                if (_HK.Count == 0)
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
            var sn0 = GroupPermutationForm.MetaCyclicGens(m,n,r).sn;
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
    AnyGroupPermFormUpTo(32); // # Total:144  Max Order:32 Time:4.490s
}