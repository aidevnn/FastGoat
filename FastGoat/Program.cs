using System.Globalization;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
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

{
    // FG.CharactersTable(Group.SemiDirectProd(FG.Abelian(3, 3), new Cn(4))).DisplayCells();
    // FG.CharactersTable(FG.DiCyclic(3)).DisplayCells();
    // a3 = b3 = c4 = 1, ab = ba, cac−1 = b, cbc−1 = a2

    // FG.CharactersTable(FG.Abelian(2)).DisplayCells();
    // FG.CharactersTable(FG.Abelian(3)).DisplayCells();
    // FG.CharactersTable(FG.Symmetric(3)).DisplayCells();
    // var table = FG.CharactersTable(Group.SemiDirectProd(new Cn(7), new Cn(3)));
    // table.DisplayCells();
    // var e = table.CnfCells[3, 3].E;
    // Console.WriteLine(e);
    // Console.WriteLine(e.Module2);
    // Console.WriteLine(e.Module);
    // Console.WriteLine(e * e.Conj);
    // Console.WriteLine(e.Re);
    // Console.WriteLine(e.Re.Module2);
    // Console.WriteLine(e.Im);
    // Console.WriteLine(e.Im.Module2);

    // var gr = Group.SemiDirectProd(new Cn(7), new Cn(3));
    // var gr = FG.Dihedral(4);
    // var gr = FG.Symmetric(3);
    var gr = FG.Alternate(5);
    var t = FG.CharactersTable(gr);
    t.DisplayCells();

    // var k = Group.Derived(gr);
    // DisplayGroup.HeadElements(k);
    // var t0 = FG.CharactersTable(k);
    // t0.DisplayCells();
    //
    // var inducedH = new Dictionary<int, Dictionary< Ep2<ZnInt,ZnInt>, Cnf>>();
    // foreach (var i in t0.IndexesTint.Values)
    // {
    //     var indGH = t.IndexesTint.ToDictionary(e => e.Key, _ => Cnf.CnfZero);
    //     foreach (var g in t.CClasses.GetRepresentatives())
    //     {
    //         foreach (var y in gr)
    //         {
    //             var g0 = gr.Op(gr.Invert(y), gr.Op(g, y));
    //             if (!k.Contains(g0))
    //                 continue;
    //
    //             var j = t0.IndexesTint[t0.CClasses.GetRepresentative(g0)];
    //             indGH[g] += t0.CnfCells[i, j].E;
    //         }
    //
    //         indGH[g] /= k.Count();
    //     }
    //
    //     inducedH[i] = indGH;
    // }
    //
    // var restrictedH = t.IndexesTint.Values.ToDictionary(
    //     i => i,
    //     i => t0.IndexesTint.ToDictionary(e => e.Key, e => t.CnfCells[i, t.IndexesTint[t.CClasses.GetRepresentative(e.Key)]].E)
    // );
    //
    // foreach (var (a, indGH) in inducedH)
    // {
    //     indGH.ToDictionary(e => t.CClasses.GetClassName(e.Key), e => e.Value).Println($"{a}");
    // }
    //
    // foreach (var (a, resGH) in restrictedH)
    // {
    //     resGH.ToDictionary(e => t0.CClasses.GetClassName(e.Key), e => e.Value).Println($"{a}");
    // }
    //
    // foreach (var (j, indGH) in inducedH)
    // {
    //     foreach (var (i, resGH) in restrictedH)
    //     {
    //         var sumG = Cnf.CnfZero;
    //         foreach (var g in gr)
    //         {
    //             var idxg = t.IndexesTint[t.CClasses.GetRepresentative(g)];
    //             sumG += indGH[t.CClasses.GetRepresentative(g)] * t.CnfCells[i, idxg].E.Conj;
    //         }
    //         sumG /= gr.Count();
    //
    //         var sumH = Cnf.CnfZero;
    //         foreach (var g in k)
    //         {
    //             var idxg = t0.IndexesTint[t0.CClasses.GetRepresentative(g)];
    //             sumH += t0.CnfCells[j, idxg].E * resGH[t0.CClasses.GetRepresentative(g)].Conj;
    //         }
    //
    //         sumH /= k.Count();
    //         Console.WriteLine(new { sumG, sumH });
    //     }
    // }
}