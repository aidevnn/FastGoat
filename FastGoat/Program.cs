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

void SolveRCohomologies<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, int r = 2)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var CN = Group.Zentrum(N);
    var (autN, innN, outN) = Group.OuterAutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, outN);
    foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)))
    {
        var L = op.ToMapElt(outN);
        var lbl = $"H{r}(G,N) Lbl{i}/{ops.Count}";
        var (hr, (_, _, br), (_, _, zr)) = ZNSolver.ReduceCohomologies(CN, G, L, r, lbl);
        // ZNSolver.DisplayCrMap($"B{r}, Z{r}, and H{r} reprs", new[] { br, zr }.Concat(hr).ToArray());
    }
}

{
    var (N, G) = (FG.Abelian(2), FG.Abelian(2));
    SolveRCohomologies(N, G, r: 1);
    SolveRCohomologies(N, G, r: 2);
    SolveRCohomologies(N, G, r: 3);
    SolveRCohomologies(N, G, r: 4);
}

{
    var (N, G) = (FG.Abelian(4), FG.Abelian(2));
    SolveRCohomologies(N, G, r: 1);
    SolveRCohomologies(N, G, r: 2);
    SolveRCohomologies(N, G, r: 3);
    SolveRCohomologies(N, G, r: 4);
}

{
    var (N, G) = (FG.Abelian(2, 2), FG.Abelian(2));
    SolveRCohomologies(N, G, r: 1);
    SolveRCohomologies(N, G, r: 2);
    SolveRCohomologies(N, G, r: 3);
    SolveRCohomologies(N, G, r: 4);
}

{
    var (N, G) = (FG.Abelian(2), FG.Abelian(2, 2));
    SolveRCohomologies(N, G, r: 1);
    SolveRCohomologies(N, G, r: 2);
    SolveRCohomologies(N, G, r: 3);
}

{
    var (N, G) = (FG.Abelian(3), FG.Abelian(3));
    SolveRCohomologies(N, G, r: 1);
    SolveRCohomologies(N, G, r: 2);
    SolveRCohomologies(N, G, r: 3);
    SolveRCohomologies(N, G, r: 4);
}

// {
//     var (N, G) = (FG.Abelian(4), FG.Abelian(2, 2));
//     SolveRCohomologies(N, G, r: 1);
//     SolveRCohomologies(N, G, r: 2);
//     SolveRCohomologies(N, G, r: 3);
// }
//
// {
//     var (N, G) = (FG.Abelian(9), FG.Abelian(3));
//     SolveRCohomologies(N, G, r: 1);
//     SolveRCohomologies(N, G, r: 2);
//     SolveRCohomologies(N, G, r: 3);
//     SolveRCohomologies(N, G, r: 4);
// }
//
// {
//     var (N, G) = (FG.Abelian(3, 3), FG.Abelian(3));
//     SolveRCohomologies(N, G, r: 1);
//     SolveRCohomologies(N, G, r: 2);
//     SolveRCohomologies(N, G, r: 3);
//     SolveRCohomologies(N, G, r: 4);
// }
