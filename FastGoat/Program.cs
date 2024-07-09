using System.CodeDom;
using System.Collections;
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
using FastGoat.UserGroup.Floats;
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
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void TensorTable<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    var gSubgrs = g.AllSubgroups();
    gSubgrs.Naming();
    var dirProd = gSubgrs.DecomposeProducts(gSubgrs.ProperNonTrivialNormalSubgroups())
        .Where(e => e.Item3).Take(1).ToArray();
    if (dirProd.Length == 0)
        return;

    var (H, K) = (dirProd[0].Item1.Representative, dirProd[0].Item2.Representative);
    var ctH = FG.CharacterTable(H);
    var ctK = FG.CharacterTable(K);
    ctH.DisplayCells();
    ctK.DisplayCells();

    var ctG = FG.CharacterTableEmpty(g);
    ctG.DerivedSubGroupLift();
    // ctG.InductionFromStabilizers();
    ctG.DisplayCells();

    foreach (var (chi1, chi2) in ctH.AllCharacters.Grid2D(ctK.AllCharacters))
    {
        var map = chi1.Gr.Grid2D(chi2.Gr).Select(e => (e, ctG.Classes.GetRepresentative(g.Op(e.t1, e.t2))))
            .DistinctBy(e => e.Item2)
            .ToDictionary(e => e.Item2, e => (Cnf?)(chi1[e.e.t1]!.Value * chi2[e.e.t2]!.Value));

        var state = ctG.AddCharacter(new(ctG.Classes, map));
        if (state == AddCharacterState.TableFull)
            break;
    }

    ctG.DisplayCells();

    Console.WriteLine("----------------------------------------------------------------------------------------------");
    
    var ctG2 = FG.CharacterTableEmpty(g);
    ctG2.DerivedSubGroupLift();
    ctG2.InductionFromStabilizers();
    ctG2.InductionFromSubGroups(gSubgrs);
    ctG2.DisplayCells();

    Console.WriteLine($"END {g.ShortName}");
    Console.WriteLine();
}

{
    TensorTable(FG.Abelian(2, 3));
    TensorTable(Product.Generate(FG.AbelianPerm(3), FG.Dihedral(4)));
}