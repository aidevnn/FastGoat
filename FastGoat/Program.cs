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

void DiagAction<Tg, Tn>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
    where Tn : struct, IElt<Tn>
    where Tg : struct, IElt<Tg>
{
    var nab = N.AbelianDirectSum();
    var subgs = nab.ElementarySubgroups();

    var autG = Group.AutomorphismGroup(G);
    var autN = Group.AutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, autN).ToHashSet(FG.EqOpByAut(G, autG, autN));
    var diagOps = ops.Where(op => G.All(g => subgs.All(sg => sg.SetEquals(sg.Select(e => op[g][e]))))).ToArray();
    foreach (var (L, k) in diagOps.Select((op, k) => (op.ToMapElt(autN), k + 1)))
    {
        var lbl = $"Lbl{k}/{diagOps.Length}/{ops.Count}";
        var decomp = subgs.Select(sg => ZNSolver.ReduceCohomologies(sg, G, L, r: 2, lbl)).ToArray();
        var sol = ZNSolver.ReduceCohomologies(N, G, L, r: 2, lbl);
        
        var nbH2 = sol.solsCohs.Length;
        var prod = decomp.Select(s => s.solsCohs.Length).Aggregate((a0, a1) => a0 * a1);
        if (nbH2 == prod)
        {
            Console.WriteLine("PASS");
            var cobsExpected = sol.solsCobs.allMaps.ToGroupMapElt();
            var cocsExpected = sol.solsCocs.allMaps.ToGroupMapElt();
            var cohsExpected = sol.solsCohs.Select(c => c.ToMapElt).ToHashSet();

            var gb = (MapGroupBase<Ep<Tg>, Tn>)cobsExpected.BaseGroup;

            var cobsSeq = decomp.Select(e => e.solsCobs.allMaps.ToGroupMapElt()).ToArray();
            var cobsActual = cobsSeq.MultiLoop().Select(l => gb.OpSeq(l)).ToHashSet();
            if (!cobsExpected.SetEquals(cobsActual))
                throw new("#Coboundaries");
            
            var cocsSeq = decomp.Select(e => e.solsCocs.allMaps.ToGroupMapElt()).ToArray();
            var cocsActual = cocsSeq.MultiLoop().Select(l => gb.OpSeq(l)).ToHashSet();
            if (!cocsExpected.SetEquals(cocsActual))
                throw new("#Cocycles");

            var cohsSeq = decomp.Select(e => e.solsCohs.Select(c => c.ToMapElt)).ToArray();
            var cohsActual = cohsSeq.MultiLoop().Select(l => gb.OpSeq(l)).ToHashSet();
            if (!cohsExpected.SetEquals(cohsActual))
                throw new("#Cohomologies");
        }
        else
            throw new("FAIL");
        
        Console.WriteLine();
    }
}

{
    var (G, N) = (FG.Abelian(2), FG.Abelian(2, 4));
    DiagAction(N, G);
}

{
    var (G, N) = (FG.Abelian(2), FG.Abelian(4, 4));
    DiagAction(N, G);
}

{
    var (G, N) = (FG.Abelian(4), FG.Abelian(2, 2));
    DiagAction(N, G);
}

{
    var (G, N) = (FG.Abelian(4), FG.Abelian(2, 4));
    DiagAction(N, G);
}

{
    var (G, N) = (FG.Abelian(3), FG.Abelian(3, 3));
    DiagAction(N, G);
}

{
    var (G, N) = (FG.Abelian(3), FG.Abelian(3, 9));
    DiagAction(N, G);
}