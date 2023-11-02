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
        ZNSolver.Reduce2Cohomologies(N, G, L, lbl);
        break;
    }
}

{
    var (N, G) = (FG.Abelian(2), FG.Abelian(2, 2, 2));
    SolveTwoCohom(N, G);
    // H2(G, N) with N:|C2| = 2 and G:|C2 x C2 x C2| = 8
    // #### Lbl1/1 |B2|:16 |Z2|:1024 |H2|:64
}

{
    var (N, G) = (FG.Abelian(2, 2), FG.Abelian(2, 2, 2));
    SolveTwoCohom(N, G);
    // H2(G, N) with N:|C2 x C2| = 4 and G:|C2 x C2 x C2| = 8
    // #### Lbl1/22 |B2|:256 |Z2|:1048576 |H2|:4096
}

{
    var (N, G) = (FG.Abelian(3), FG.Abelian(3, 3));
    SolveTwoCohom(N, G);
    // H2(G, N) with N:|C3| = 3 and G:|C3 x C3| = 9
    // #### Lbl1/1 |B2|:729 |Z2|:19683 |H2|:27
}

{
    var (N, G) = (FG.Abelian(3, 3), FG.Abelian(3, 3));
    SolveTwoCohom(N, G);
    // H2(G, N) with N:|C3 x C3| = 9 and G:|C3 x C3| = 9
    // #### Lbl1/33 |B2|:531441 |Z2|:387420489 |H2|:729
}
