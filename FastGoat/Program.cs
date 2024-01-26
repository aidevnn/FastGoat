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
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void Cn(int n, bool details = false)
{
    Graph.RunToddCoxeterAlgo($"a{n}", details);
}

void Abelian(int[] seq, bool details = false)
{
    var n = seq.Length.Range();
    var gens = n.Select(i => ((char)('a' + i)).ToString()).ToArray();
    var allCombs = n.SelectMany(i => n.Where(j => j > i).Select(j => (gens[i], gens[j]))).ToArray();
    var relators = gens.Select((g, i) => $"{g}{seq[i]}").Concat(allCombs.Select(e => $"{e.Item1}{e.Item2}={e.Item2}{e.Item1}"))
        .Glue(", ");

    Console.WriteLine("{0}:{1}", seq.Glue(" x ", "C{0}"), relators);
    Graph.RunToddCoxeterAlgo(relators, details);
}

void Dihedral(int n, bool details = false)
{
    Console.WriteLine($"D{2 * n}");
    Graph.RunToddCoxeterAlgo($"a{n},b2,abab", details);
}

void DiCyclic(int n, bool details = false)
{
    Console.WriteLine($"Dic{4 * n}");
    Graph.RunToddCoxeterAlgo($"a{n} = b2, b2 = abab", details);
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
            Graph.RunToddCoxeterAlgo($"a{m}, b{n}, b-1ab = a{r}", details);
        }
    }
}

void Coincidences(bool details = false)
{
    Graph.RunToddCoxeterAlgo("aba-1 = b2, bab-1 = a2", details);
    Graph.RunToddCoxeterAlgo("a3,b3,aba2b", details);
    Graph.RunToddCoxeterAlgo("a", "a3,b3,aba2b", details);
}

void PermGroup45(bool details = false)
{
    Graph.RunToddCoxeterAlgo("a3,b3,abab", details); // Alt4
    Graph.RunToddCoxeterAlgo("a4,b3,abab", details); // Symm4
    Graph.RunToddCoxeterAlgo("a3,b3,c3,abab,bcbc,acac", details); // Alt5
    Graph.RunToddCoxeterAlgo("a5, b2, abababab, a2ba-2ba2ba-2b", details); // Symm5
}

void PermGroup6()
{
    Graph.RunToddCoxeterAlgo("a3,b3,c3,d3,abab,acac,adad,bcbc,bdbd,cdcd", details: false); // Alt6
    // Step:361 NbClasses:360
    // Time:371ms 

    Graph.RunToddCoxeterAlgo("a6, b2, ababababab, a2ba-2ba2ba-2b, baba-1baba-1baba-1", details:false); // Symm6
    // Step:721 NbClasses:720
    // Time:979ms
}

void PermGroup7()
{
    Graph.RunToddCoxeterAlgo("a3, b3, c3, d3, e3, abab, acac, adad, aeae, bcbc, bdbd, bebe, cdcd, cece, dede", details:false); // Alt7 dede
    // Step:2521 NbClasses:2520
    // Time:20.648s

    Graph.RunToddCoxeterAlgo("a7, b2, abababababab, baba-1baba-1baba-1, a2ba-2ba2ba-2b", details: false); // Symm7
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

void Tests()
{
    Graph.DefiningRelatorsOfGroup(FG.Abelian(2));
    Graph.DefiningRelatorsOfGroup(FG.Abelian(5));
    Graph.DefiningRelatorsOfGroup(FG.Abelian(2, 2));
    Graph.DefiningRelatorsOfGroup(FG.Dihedral(4));
    
    for (int k = 2; k <= 8; ++k)
    {
        Graph.DefiningRelatorsOfGroup(FG.DiCyclic(k));
        Graph.DefiningRelatorsOfGroup(FG.DiCyclicSdp(k));
    }
    
    for (int n = 4; n <= 7; n++)
    {
        Graph.DefiningRelatorsOfGroup(FG.Alternate(n));
        Graph.DefiningRelatorsOfGroup(FG.Symmetric(n));
    }
}

void L2p()
{
    Graph.DefiningRelatorsOfGroup(FG.L2p(11));
    /*
       |L2(11)| = 660

       #  Time:518ms
       All Relators
           a3
           b2
           ababababababababababab
           baba-1baba-1baba-1baba-1baba-1
           abababa-1ba-1ba-1babababa-1ba-1ba-1b
     */
    
    Graph.RunToddCoxeterAlgo("a3,b2,ababababababababababab,baba-1baba-1baba-1baba-1baba-1,abababa-1ba-1ba-1babababa-1ba-1ba-1b", details: false);
    // Step:683 NbClasses:660
    // #  Time:2.355s
    
    Graph.DefiningRelatorsOfGroup(FG.L2p(13));
    /* |L2(13)| = 1092

       #  Time:210ms
       All Relators
           a3
           b2
           ababababababababababababab
           abababa-1baba-1babababa-1baba-1b
           ababababa-1bababababa-1ba-1bababa-1ba-1b
     */
    
    Graph.RunToddCoxeterAlgo("a3,b2,ababababababababababababab,abababa-1baba-1babababa-1baba-1b,ababababa-1bababababa-1ba-1bababa-1ba-1b", details: false);
    // Step:1099 NbClasses:1092
    // #  Time:3.939s
    
    Graph.DefiningRelatorsOfGroup(FG.L2p(17));
    /* |L2(17)| = 2448
       
       #  Time:547ms
       All Relators
           a3
           b2
           ababababa-1ba-1bababa-1baba-1bababa-1ba-1b
           abababababa-1ba-1ba-1babababababa-1ba-1ba-1b
     */
    
    Graph.RunToddCoxeterAlgo("a3,b2,ababababa-1ba-1bababa-1baba-1bababa-1ba-1b,abababababa-1ba-1ba-1babababababa-1ba-1ba-1b", details: false);
    // Step:2499 NbClasses:2448
    // #  Time:14.651s
}

void u33()
{
    var s28 = new Sn(28);
    var a2 = s28[(1, 5, 7, 3, 12, 24, 11), (2, 23, 4, 27, 13, 14, 26), (6, 20, 18, 8, 25, 21, 28), (9, 10, 17, 15, 22, 16, 19)];
    var b2 = s28[(3, 4), (5, 17, 7, 16, 8, 20, 6, 13), (9, 19, 11, 14, 12, 18, 10, 15), (21, 23, 26, 28, 24, 22, 27, 25)];
    var u3_3pg = Group.Generate("U3(3)pg", s28, a2, b2);
    Graph.DefiningRelatorsOfGroup(u3_3pg);
    /* |U3(3)pg| = 6048
       
       #  Time:2.945s
       All Relators
           a8
           b7
           ba-1ba-1ba-1
           a4ba3b-1a-1b2
           a2b-1ab2abab-1a-1b
     */
    
    Graph.RunToddCoxeterAlgo("a8,b7,ba-1ba-1ba-1,a4ba3b-1a-1b2,a2b-1ab2abab-1a-1b", details: false);
    // Step:8541 NbClasses:6048
    // #  Time:1m46s
}

{
    // Graph.DefiningRelatorsOfGroup(FG.GLnp(2, 3)); // random
    // Graph.RunTC("a8,b2,a2ba-1ba-1ba", details: false);
    // L2p();
    // DisplayGroup.HeadElements(FG.AbelianWg(2, 3));
    // DisplayGroup.HeadElements(FG.DihedralWg(5));

    Graph.RunToddCoxeterAlgo("a3,b2,c2,abababab,acac,bcbcbc", details: false);
}