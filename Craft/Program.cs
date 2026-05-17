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

Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>> FactorGroup<T>(ConcreteGroup<T> G, ConcreteGroup<T> H)
    where T : struct, IElt<T>
{
    Console.WriteLine($"G:{G.ShortName} H:{H.ShortName}");
    var allFactors = GroupCraft.AllFactors(G, H);
    var HIsNormal = Group.IsNormalSubgroup(G, H);
    var facts = allFactors.Where(e => HIsNormal || Group.IsNormalSubgroup(G, e.Key)).ToDictionary();
    Console.WriteLine("############ START");
    foreach (var (sol, _) in facts)
    {
        var K = Group.Generate("K", G, sol.GetGenerators().ToArray());
        Console.WriteLine($"G:{G.ShortName} H:{H.ShortName} K:{K.ShortName}");
        var KIsNormal = Group.IsNormalSubgroup(G, K);
        Console.WriteLine($"H Is NormalSubgroup:{HIsNormal}");
        Console.WriteLine($"K Is NormalSubgroup:{KIsNormal}");
        var prodType = (!HIsNormal && !KIsNormal) || K.Intersect(H).Count() != 1
            ? "Not a Product"
            : (HIsNormal && !KIsNormal) || (!HIsNormal && KIsNormal)
                ? "Semi Direct Product"
                : "Direct Product";
        Console.WriteLine($"Is {prodType}");
        Console.WriteLine();
    }

    if (facts.Count == 0)
    {
        Console.WriteLine("No Factor Direct");
        Console.WriteLine();
    }

    Console.WriteLine("############ END");
    Console.WriteLine();

    return facts;
}

void TestFactorGroup()
{
    // < a,b | b3, a10, baba-1 >
    // Generators of F(3x:10)2 in GL(2,31)
    // var G = FG.WordGroup("F(3x:10)2", "b3, a10, baba-1");
    // < a,b | a8, b2, a2ba2b, aba-1baba-1b > 
    // Generators of MM16 x: C2 in GL(4,3)
    // var G = FG.WordGroup("MM16 x: C2", "a8, b2, a2ba2b, aba-1baba-1b");
    // < a,b,c | b2, c2, abab, a2bcbc, a3ca-1c >
    // Generators of D16 x: C2 in GL(4,5)
    var G = FG.WordGroup("D16 x: C2", "b2, c2, abab, a2bcbc, a3ca-1c");
    // < a,b,c | a4, b4, c2, abab-1, caca-1, cbcb-1 >
    // Generators of C2 x M(4x:4)3 in GL(4,5)
    // var G = FG.WordGroup("C2 x M(4x:4)3", "a4, b4, c2, abab-1, caca-1, cbcb-1");

    DisplayGroup.HeadElements(G);
    var subs = G.AllSubgroups();
    subs.Naming();
    foreach (var e in subs.Where(e => !e.IsTrivial && e.Order != G.Count()))
    {
        var H = e.Representative;
        Console.WriteLine($"############### H:{H.ShortName,-40} gensH:{H.GetGenerators().ToXSet()}");
        DisplayGroup.HeadElements(H);
        FactorGroup(G, H);
        Console.WriteLine("############### END");
        Console.WriteLine();
    }
}

{
    for (int ord = 1; ord <= 64; ord++)
    {
        foreach (var G in FG.AllGroupsOfOrder(ord).Where(e => e.GroupType == GroupType.NonAbelianGroup))
        {
            var subs = G.AllSubgroups();
            var prods = subs.DecomposeProducts(subs.ProperNonTrivialNormalSubgroups());
            prods = prods
                .Concat(prods.Where(e => e.isDirectProduct).Select(e => (lhs: e.rhs, rhs: e.lhs, e.isDirectProduct)))
                .ToList();
            foreach (var (lhs, seq) in prods.GroupBy(e => e.lhs).ToDictionary(e => e.Key, e => e.ToArray()))
            {
                var factorsExpected = seq.SelectMany(e => e.rhs.Conjugates.Select(c => c.ToSet())).ToXSet();
                var factorsFounds = FactorGroup(G, lhs.Representative).SelectMany(e => e.Value).ToXSet();
                Console.WriteLine($"facts:{factorsFounds.Count} / {factorsExpected.Count}");
                Console.WriteLine();
                if (factorsFounds.Count != factorsExpected.Count || !factorsFounds.Equals(factorsExpected))
                    throw new();
            }
        }
    }
}