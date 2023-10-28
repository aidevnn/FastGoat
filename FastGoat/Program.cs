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
        var (nb2Cobs, map2Cobs) = ZNSolver.Reduce2Coboundaries(N, G, L);
        var (nb2Cocs, map2Cocs) = ZNSolver.Reduce2Cocycles(N, G, L);
        var nbCohs = nb2Cocs / nb2Cobs;
        Console.WriteLine($"#### {lbl} |B2|:{nb2Cobs} |Z2|:{nb2Cocs} |H2|:{nbCohs}");
        Console.WriteLine();

        var max = BigInteger.Pow(N.Count(), G.Count() - 1);
        if (details && max < 7000)
        {
            var all = CocyclesDFS.TwoCocycles(N, G, L, lbl);
            all.AllTwoCocycles();
            if (all.AllCoboundaries.Count != nb2Cobs || all.AllCosets.Count != nbCohs)
                throw new("############### Error in order H2(G, N)"); 
        }
    }
}

{
    var (N, G) = (FG.Abelian(4), FG.Abelian(2));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(2, 2), FG.Abelian(2));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(2, 4), FG.Abelian(2));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(2, 2, 2), FG.Abelian(2));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(2, 2, 4), FG.Abelian(2));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(4, 4), FG.Abelian(2));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(4), FG.Abelian(2, 2));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(2, 2), FG.Abelian(2, 2));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(3), FG.Abelian(3));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(6), FG.Abelian(3));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(9), FG.Abelian(3));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(3, 3), FG.Abelian(3));
    SolveTwoCohom(N, G, details: true);
}

{
    var (N, G) = (FG.Abelian(2, 4), FG.Abelian(2, 2));
    SolveTwoCohom(N, G, details: true);
}