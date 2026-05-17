using System.Security.Cryptography;
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

void FactorGroup<T>(ConcreteGroup<T> G, ConcreteGroup<T> H) where T : struct, IElt<T>
{
    Console.WriteLine($"G:{G.ShortName} H:{H.ShortName}");
    var allFactors = GroupCraft.AllFactors(G, H);
    var HIsNormal = Group.IsNormalSubgroup(G, H);
    var facts = allFactors.Keys.Where(e => HIsNormal || Group.IsNormalSubgroup(G, e)).Select(g => g.ToSet()).ToList();
    foreach (var sol in facts.Take(1))
    {
        var K = Group.Generate("K", G, sol.Generators.ToArray());
        DisplayGroup.HeadElements(K);
        Console.WriteLine($"G:{G.ShortName} H:{H.ShortName} K:{K.ShortName} allFactors:{facts.Count}");
        var KIsNormal = Group.IsNormalSubgroup(G, K);
        Console.WriteLine($"H Is NormalSubgroup:{HIsNormal}");
        Console.WriteLine($"K Is NormalSubgroup:{KIsNormal}");
        Console.WriteLine($"Is DirectProduct   :{(HIsNormal || KIsNormal) && K.Intersect(H).Count() == 1}");
        Console.WriteLine();
    }

    if (facts.Count == 0)
    {
        Console.WriteLine("No Factor Direct");
        Console.WriteLine();
    }
}

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