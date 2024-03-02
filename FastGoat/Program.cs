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

void RCohomologyCyclicgroup(int p, int r)
{
    if (r < 1 || !Primes10000.Contains(p))
        throw new();

    var cp = new Cn(p);
    var autCp = Group.AutomorphismGroup(cp);
    var L = new MapElt<ZnInt, Automorphism<ZnInt>>(cp, autCp); // trivial

    Console.WriteLine($"H{r}({cp}, {cp})");
    var solsCobs = ZNSolver.ReduceCoboundaries(cp, cp, L, r - 1);
    var solsCocs = ZNSolver.ReduceCocycles(cp, cp, L, r);
    var (nbCobs, sredCobs, mapCobs) = solsCobs;
    var (nbCocs, sredCocs, mapCocs) = solsCocs;
    var nb2Cohs = nbCocs / nbCobs;
    Console.WriteLine($"|B{r}|:{nbCobs} |Z{r}|:{nbCocs} |H{r}|:{nb2Cohs}");

    if (nb2Cohs != p)
        throw new();
    
    Console.WriteLine();
}

{
    foreach (var p in new[] { 2, 3, 5, 7 })
    {
        for (int r = 1; p.Pow(r + 1) < 1000; ++r)
            RCohomologyCyclicgroup(p, r); // H^r(Cp, Cp) = Cp
    }
}