using System.Collections;
using System.ComponentModel;
using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
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
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void TwoCohomologyCyclicgroup(int p)
{
    if (!Primes10000.Contains(p))
        throw new();

    var cp = new Cn(p);
    var autCp = Group.AutomorphismGroup(cp);
    var L = new MapElt<ZnInt, Automorphism<ZnInt>>(cp, autCp); // trivial

    var (cohs1, _, _) = ZNSolver.ReduceCohomologies(cp, cp, L, 1);
    var (cohs2, sysCobs2, _) = ZNSolver.ReduceCohomologies(cp, cp, L, 2);
    ZNSolver.DisplayCrMap("H2", cohs2);

    var cobs2 = sysCobs2.allMaps;
    var ind = cobs2.ZZero.Indeterminates;
    var set1 = cohs1.Grid2D().Select(e => e.t1.Mul(e.t2).Recreate(ind)).ToHashSet();
    ZNSolver.DisplayCrMap(set1.Prepend(cobs2).ToArray());
    var set2 = new HashSet<CrMap<ZnInt, ZnInt>>();
    foreach (var m1 in set1)
    {
        var m2 = ZNSolver.SysRepresentative(cobs2.Add(m1));
        set2.Add(m2);
    }
    
    ZNSolver.DisplayCrMap(set2.ToArray());
    Console.WriteLine();
}

void TestSylowsAndCohomology<Tn, Tg>(ConcreteGroup<Tg> G, ConcreteGroup<Tn> N)
    where Tn : struct, IElt<Tn>
    where Tg : struct, IElt<Tg>
{
    var autG = Group.AutomorphismGroup(G);
    var autN = Group.AutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, autN).ToHashSet(FG.EqOpByAut(G, autG, autN));
    var allSubgs = G.AllSubgroups();
    allSubgs.Naming();
    var pSylows = allSubgs.AllSylows().Values.SelectMany(cj => cj).Select(sg => sg.Representative).ToArray();
    foreach (var L in ops)
    {
        var L0 = L.ToMapElt(autN);
        pSylows.Select(sg => ZNSolver.ReduceCohomologies(N, sg, L0)).ToArray();
        ZNSolver.ReduceCohomologies(N, G, L0);
        Console.WriteLine();
    }

    Console.WriteLine();
}

{
    TwoCohomologyCyclicgroup(2);
    TwoCohomologyCyclicgroup(3);
    
    // TestSylowsAndCohomology(FG.Symmetric(3), FG.Abelian(2));
    // TestSylowsAndCohomology(Group.SemiDirectProd(new Cn(3), new Cn(4)), FG.Abelian(2));
    // TestSylowsAndCohomology(FG.Alternate(4), FG.Abelian(2));
}
