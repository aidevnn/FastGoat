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

void SolveTwoCohom<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, int start = 1, bool details = false)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var autN = Group.AutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, autN);
    foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)).Skip(start - 1))
    {
        var L = op.ToMapElt(autN);
        var lbl = $"Lbl{i,2}/{ops.Count}";
        if (details)
            ExplicitTwoCohom(N, G, L);
        
        ZNSolver.Reduce2Cohomologies(N, G, L, lbl, details);
    }
}

void ExplicitTwoCohom<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var (nb2Cobs, _, maps2Cobs) = ZNSolver.Reduce2Coboundaries(N, G, L);
    var (nb2Cocs, _, maps2Cocs) = ZNSolver.Reduce2Cocycles(N, G, L);
    if (nb2Cocs > 1000)
        return;

    var g2Cobs = maps2Cobs.ToGroupMapElt("B2");
    var g2Cocs = maps2Cocs.ToGroupMapElt("Z2");
    var g2Cohs = g2Cocs.Over(g2Cobs, "H2");
    ZNSolver.DisplayCrMap(maps2Cobs);
    DisplayGroup.Head(g2Cohs);
    CocyclesDFS.DisplayMapElt("H2", g2Cohs.Select(e => e.X).ToArray());
}

// H2(D8, C4)
{
    var (N, G) = (new Cn(4), FG.Dihedral(4));
    SolveTwoCohom(N, G);
}

//############# Lbl 1/4      #############
// H2(G, N) with N:|C4| = 4 and G:|D8| = 8
// #### Lbl 1/4 |B2|:4096 |Z2|:32768 |H2|:8
// 
// B2(G,N):4096
// Z2(G,N):32768
// *H2(G,N):8 Expected:8
// ############# Lbl 2/4      #############
// H2(G, N) with N:|C4| = 4 and G:|D8| = 8
// #### Lbl 2/4 |B2|:2048 |Z2|:16384 |H2|:8
// 
// B2(G,N):2048
// Z2(G,N):16384
// %H2(G,N):32 Expected:8
// ############# Lbl 3/4      #############
// H2(G, N) with N:|C4| = 4 and G:|D8| = 8
// #### Lbl 3/4 |B2|:1024 |Z2|:16384 |H2|:16
// 
// B2(G,N):1024
// Z2(G,N):16384
// *H2(G,N):16 Expected:16
// ############# Lbl 4/4      #############
// H2(G, N) with N:|C4| = 4 and G:|D8| = 8
// #### Lbl 4/4 |B2|:2048 |Z2|:16384 |H2|:8
// 
// B2(G,N):2048
// Z2(G,N):16384
// *H2(G,N):8 Expected:8
// 

// BATCH Test
{
    var (N, G) = (FG.Abelian(4), FG.Abelian(2));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(8), FG.Abelian(2));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(16), FG.Abelian(2));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(4), FG.Abelian(2, 2));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(8), FG.Abelian(2, 2));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(2, 4), FG.Abelian(2));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(4, 4), FG.Abelian(2));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(2, 8), FG.Abelian(2));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(4, 8), FG.Abelian(2));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(4), FG.Abelian(2, 2));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(4), FG.Abelian(4));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(8), FG.Abelian(4));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(9), FG.Abelian(3));
    SolveTwoCohom(N, G);
}

{
    var (N, G) = (FG.Abelian(27), FG.Abelian(3));
    SolveTwoCohom(N, G);
}

// {
//     var (N, G) = (FG.Abelian(2, 4), FG.Abelian(4));
//     SolveTwoCohom2(N, G, details:true); // BUG
// }
//
// {
//     var (N, G) = (FG.Abelian(3, 9), FG.Abelian(3));
//     SolveTwoCohom2(N, G, details:true); // BUG
// }
//
// {
//     var (N, G) = (FG.Abelian(2, 4), FG.Abelian(2, 4));
//     SolveTwoCohom(N, G); // TODO
// }
