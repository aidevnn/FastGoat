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
        Console.WriteLine($"############# {lbl,-12} #############");
        var nbCobs = LSolveTwoCohom(N, G, L, r: 1);
        var nbCocs = LSolveTwoCohom(N, G, L, r: 2);
        Console.WriteLine($"#### {lbl} |B2|:{nbCobs} |Z2|:{nbCocs} |H2|:{nbCocs / nbCobs}");
        Console.WriteLine();
        
        var max = BigInteger.Pow(N.Count(), G.Count() - 1);
        if (details && max < 7000)
        {
            var all = CocyclesDFS.TwoCocycles(N, G, L, lbl);
            all.AllTwoCocycles();
            if (all.AllCoboundaries.Count != nbCobs || all.AllCosets.Sum(c => c.Value.Count) != nbCocs)
                throw new("###############");

            // foreach (var c in all.AllCoboundaries)
            //     c.Map.OrderKeys(G).Println("Coboundary");

            // foreach (var c in all.AllCosets.Values.SelectMany(vs => vs))
            //     c.Map.OrderKeys(G).Println("Cocycle");
        }
    }
}

int LSolveTwoCohom<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L, int r)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var (cr, cnext) = ZNSolver.LRCochain(N, G, L, r);
    // var Cnextnext = Dr(Cnext); // always zero

    cr.OrderKeys(G).Println($"C{r}");
    cnext.OrderKeys(G).Println($"C{r + 1}");
    // Cnextnext.OrderKeys(G).Println($"C{r + 2}"); 

    if (r == 1)
    {
        var (nb2Cobs, _) = ZNSolver.Reduce2Coboundaries(N, G, L);
        return nb2Cobs;
    }

    if (r == 2)
    {
        var (nb2Cocs, _) = ZNSolver.Reduce2Cocycles(N, G, L);
        return nb2Cocs;
    }

    return -1;
}

{
    var (N, G) = (FG.Abelian(4, 4), FG.Abelian(2));
    SolveTwoCohom(N, G, details: true);
}
