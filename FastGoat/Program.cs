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

void TwoCocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
    where Tn : struct, IElt<Tn>
    where Tg : struct, IElt<Tg>
{
    var CN = Group.Zentrum(N);
    var autN = Group.AutomorphismGroup(CN);
    var allOps = Group.AllHomomorphisms(G, autN);
    int lbl = 1;
    foreach (var hom in allOps)
    {
        var L = new MapElt<Tg, Automorphism<Tn>>(G, autN, new(hom.HomMap));
        L.map.Println();
        Console.WriteLine();
        var sols = Solver.SolveEq2Cocycles(CN, G, L);
        Console.WriteLine($"New meth Sol{lbl} total {sols.ToList().Count}");
        var homol = new CocyclesDFS.TwoCocyclesDFS<Tn, Tg>(N, G, L, $"Lbl{++lbl}/{allOps.Count}");
        var all2Cocycles = homol.AllTwoCocycles();
        Console.WriteLine();
    }
    
    Console.WriteLine();
    CocyclesDFS.All_2_Cocycles_N_G(CN, G);
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
    TwoCocycles(N, G);
}

{
    var (N, G) = (new Cn(4), new Cn(2));
    TwoCocycles(N, G);
}

{
    var (N, G) = (FG.Abelian(2, 2), new Cn(2));
    TwoCocycles(N, G);
}

{
    var (N, G) = (new Cn(2), FG.Abelian(2, 2));
    TwoCocycles(N, G);
}

{
    var (N, G) = (FG.Abelian(4, 4), new Cn(2));
    TwoCocycles(N, G);
}

{
    var (N, G) = (FG.Abelian(4, 2), new Cn(4));
    TwoCocycles(N, G);
}

// {
//     Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
//     var (N, G) = (new Cn(4), FG.Dihedral(4));
//     TwoCocycles(N, G);
// }

// {
//     Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
//     var (N, G) = (FG.Quaternion(16), new Cn(2));
//     TwoCoboundaries(N, G);
//     TwoCocycles(N, G);
// }