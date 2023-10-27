using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
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

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void SolveTwoCohom<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, bool details = false)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var autN = Group.AutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, autN);
    foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)))
    {
        var L = op.ToMapElt(autN);
        var lbl = $"Lbl{i}/{ops.Count}";
        var prod = LSolveTwoCohom(N, G, L, r: 1, lbl, details);
        LSolveTwoCohom(N, G, L, r: 2, lbl, details);
        
        var max = BigInteger.Pow(N.Count(), G.Count() - 1);
        if (details && max < 7000)
        {
            var all = CocyclesDFS.TwoCocycles(N, G, L, lbl);
            all.AllTwoCocycles();
            if (all.AllCoboundaries.Count != prod)
                throw new("###############");
            
            // foreach (var c in all.AllCoboundaries)
            //     c.Map.OrderKeys(G).Println("Coboundary");
            
            // foreach (var c in all.AllCosets.Values.SelectMany(vs => vs))
            //     c.Map.OrderKeys(G).Println("Cocycle");
        }
    }
}

int LSolveTwoCohom<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L, int r, string lbl = "test",
    bool details = false)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var (cr, cnext) = ZNSolver.LRCochain(N, G, L, r);
    // var Cnextnext = Dr(Cnext); // always zero

    Console.WriteLine($"############# {lbl,-12} #############");
    cr.OrderKeys(G).Println($"C{r}");
    cnext.OrderKeys(G).Println($"C{r + 1}");
    // Cnextnext.OrderKeys(G).Println($"C{r + 2}");

    if (r == 1)
    {
        var (prod, map) = ZNSolver.Reduce2Coboundaries(N, G, L);
        map.OrderKeys(G).Println("All 2Coboundaries");
        return prod;
    }

    if (r == 2)
    {
        ZNSolver.Reduce2Cocycles(N, G, L);
        return -1;
    }

    return -1;
}

{
    var (N, G) = (FG.Abelian(2, 2), FG.Abelian(2, 2));
    // RCochain(N, G, details: true);
    // RCochain(N, G, details: true);
    SolveTwoCohom(N, G, details: true);
}

// {
//     var N = FG.Abelian(4, 6);
//     var Nab = new AbelianDirectSum<Ep<ZnInt>>(N);
//     DisplayGroup.HeadElements(N);
//     DisplayGroup.HeadElements(Nab.AbCanonic);
//     DisplayGroup.HeadElements(Nab.AbElementaries);
//     
//     foreach (var n in N.Order())
//     {
//         Console.WriteLine($"n = {n} --Can--> {Nab.GEltToCan(n)} --Elem--> {Nab.GEltToElem(n)}");
//     }
// }