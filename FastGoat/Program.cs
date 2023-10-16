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
using FastGoat.UserGroup.GModuleM;
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

void TwoCocycles<Tn, Tg>(ConcreteGroup<Tn> N1, ConcreteGroup<Tg> G)
    where Tn : struct, IElt<Tn>
    where Tg : struct, IElt<Tg>
{
    var N = Group.Zentrum(N1);
    var map = Solver.MapCocycles(N, G);
    var sys = Solver.TwoCocycleCondition(map);

    var autN = Group.AutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, autN);
    var hom = ops.First();
    var L = new MapElt<Tg, Automorphism<Tn>>(G, autN, new(hom.HomMap));
    
    var eq = sys.Order().First();
    var sols = Solver.SolveEq2Cocycle(eq, L);
    hom.HomMap.Println("Solve for L");
    foreach (var sol in sols)
    {
        Console.WriteLine($"Sol:{sol.Glue(",")}");
        var newSys = sys.Select(gz0 => gz0.Substitute(sol).Act(L)).Where(gz0 => !gz0.IsZero()).Distinct().ToHashSet();
        newSys.Order().Println($"New Sys:{newSys.Count()}");
    }

    Console.WriteLine();
}

void TwoCoboundaries<Tn, Tg>(ConcreteGroup<Tn> N1, ConcreteGroup<Tg> G)
    where Tn : struct, IElt<Tn>
    where Tg : struct, IElt<Tg>
{
    var N = Group.Zentrum(N1);
    var map = Solver.MapCoboundaries(N, G);
    Solver.TwoCocycleCondition(map); // always empty system
}

{
    var (N, G) = (new Cn(2), new Cn(2));
    TwoCoboundaries(N, G);
    TwoCocycles(N, G);
}

{
    var (N, G) = (new Cn(4), new Cn(2));
    TwoCoboundaries(N, G);
    TwoCocycles(N, G);
}

{
    var (N, G) = (FG.Abelian(2, 2), new Cn(2));
    TwoCoboundaries(N, G);
    TwoCocycles(N, G);
}

{
    var (N, G) = (FG.Abelian(2, 2, 2), new Cn(2));
    TwoCoboundaries(N, G);
    TwoCocycles(N, G);
}

{
    var (N, G) = (new Cn(2), FG.Abelian(2, 2));
    TwoCoboundaries(N, G);
    TwoCocycles(N, G);
}

{
    var (N, G) = (FG.Abelian(4, 4), new Cn(2));
    TwoCoboundaries(N, G);
    TwoCocycles(N, G);
}

{
    Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
    var (N, G) = (FG.Quaternion(16), new Cn(2));
    TwoCoboundaries(N, G);
    TwoCocycles(N, G);
}