using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.TC;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void Cn(int n, bool details = false)
{
    Graph.RunTC($"a{n}", details);
}

void Abelian(int[] seq, bool details = false)
{
    var n = seq.Length.Range();
    var gens = n.Select(i => ((char)('a' + i)).ToString()).ToArray();
    var allCombs = n.SelectMany(i => n.Where(j => j > i).Select(j => (gens[i], gens[j]))).ToArray();
    var relators = gens.Select((g, i) => $"{g}{seq[i]}").Concat(allCombs.Select(e => $"{e.Item1}{e.Item2}={e.Item2}{e.Item1}"))
        .Glue(", ");

    Console.WriteLine("{0}:{1}", seq.Glue(" x ", "C{0}"), relators);
    Graph.RunTC(relators, details);
}

void Dihedral(int n, bool details = false)
{
    Console.WriteLine($"D{2 * n}");
    Graph.RunTC($"a{n},b2,abab", details);
}

void DiCyclic(int n, bool details = false)
{
    Console.WriteLine($"Dic{4 * n}");
    Graph.RunTC($"a{n} = b2, b2 = abab", details);
}

void MetaCyclicSdp(int order, bool details = false)
{
    var ms = Dividors(order).Where(d => d > 1).ToArray();
    foreach (var m in ms)
    {
        var n = order / m;
        var rs = SolveAll_k_pow_m_equal_one_mod_n(m, n);
        foreach (var r in rs)
        {
            Console.WriteLine($"MtCyc({m},{n},{r})");
            Graph.RunTC($"a{m}, b{n}, b-1ab = a{r}", details);
        }
    }
}

void Coincidences(bool details = false)
{
    Graph.RunTC("aba-1 = b2, bab-1 = a2", details);
    Graph.RunTC("a3,b3,aba2b", details);
    Graph.RunTC("a", "a3,b3,aba2b", details);
}

void PermGroup45(bool details = false)
{
    Graph.RunTC("a3,b3,abab", details); // Alt4
    Graph.RunTC("a4,b3,abab", details); // Symm4
    Graph.RunTC("a3,b3,c3,abab,bcbc,acac", details); // Alt5
    Graph.RunTC("a5, b2, abababab, a2ba-2ba2ba-2b", details); // Symm5
}

void PermGroup6()
{
    Graph.RunTC("a3,b3,c3,d3,abab,acac,adad,bcbc,bdbd,cdcd", details: false); // Alt6
    // Step:361 NbClasses:360
    // Time:371ms 

    Graph.RunTC("a6, b2, ababababab, a2ba-2ba2ba-2b, baba-1baba-1baba-1", details:false); // Symm6
    // Step:721 NbClasses:720
    // Time:979ms
}

void PermGroup7()
{
    Graph.RunTC("a3, b3, c3, d3, e3, abab, acac, adad, aeae, bcbc, bdbd, bebe, cdcd, cece, dede", details:false); // Alt7 dede
    // Step:2521 NbClasses:2520
    // Time:20.648s

    Graph.RunTC("a7, b2, abababababab, baba-1baba-1baba-1, a2ba-2ba2ba-2b", details: false); // Symm7
    // Step:5787 NbClasses:5040
    // Time:44.875s
}

{
    // Cyclics
    // for (int n = 1; n <= 8; n++) Cn(n);
    
    // Abelians
    // Abelian([2, 2]);
    // Abelian([2, 3]);
    // Abelian([2, 4]);
    // Abelian([2, 4, 6]);
    
    // Dihedrals
    // for (int n = 2; n <= 8; n++) Dihedral(n);
    
    // DiCyclics
    // for (int n = 2; n <= 8; n++) DiCyclic(n);
    
    // MetaCyclicSdp
    // for (int o = 4; o < 30; o += 2) MetaCyclicSdp(o);
    
    // PermGroup45();
    // PermGroup6();
    // PermGroup7();
}

{
    Graph.ClassesFromGroup(FG.Abelian(2));
    Graph.RunTCtable("a2");
    Graph.ClassesFromGroup(FG.Abelian(5));
    Graph.RunTCtable("a5");
    Graph.ClassesFromGroup(FG.Abelian(2, 2));
    Graph.RunTCtable("a2,b2,ab=ba");
    Graph.ClassesFromGroup(FG.Dihedral(4));
    Graph.RunTCtable("a4,b2,abab");
    
    for (int k = 2; k <= 8; ++k)
    {
        Graph.ClassesFromGroup(FG.DiCyclic(k));
        Graph.ClassesFromGroup(FG.DiCyclicSdp(k));
    }

    for (int n = 4; n <= 7; n++)
    {
        Graph.ClassesFromGroup(FG.Alternate(n));
        Graph.ClassesFromGroup(FG.Symmetric(n));
    }
}
