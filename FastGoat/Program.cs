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

void TestOrth<T>(ConcreteGroup<T> g, bool indStab = true, params (int dim, int[] linIdx)[] infos)
    where T : struct, IElt<T>
{
    var ct = FG.CharacterTableEmpty(g);
    ct.DerivedSubGroupLift();
    if (indStab)
        ct.InductionFromStabilizers();

    Logger.SetOff();
    DisplayGroup.HeadOrders(g);
    ct.DisplayCells(tableOnly: true);
    if (ct.TableComplete)
        return;

    ct.SolveSumSquare();
    ct.SolveOrthogonality(infos);
    ct.DisplayCells();
    Console.Beep();
}

void LinChisMulGroup<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    Console.WriteLine(g.ShortName);
    var ct = FG.CharacterTable(g);
    if (!ct.TableComplete)
        ct.SolveOrthogonality();

    var chG = Group.MulGroup($"Ch({g})", ct.DoneChis.Where(chi => chi.IsLinear).ToArray());
    var gtype = Group.AbelianGroupType(chG);
    chG.Name = gtype.Glue(" x ", "C{0}");
    DisplayGroup.HeadElements(chG);
    DisplayGroup.Generators(chG);
}

void runOrth()
{
    GlobalStopWatch.Restart();
    TestOrth(FG.Quaternion(16), indStab: false);
    TestOrth(FG.MetaCyclicSdp(4, 4, 3), indStab: false);
    TestOrth(FG.DiCyclic(5), indStab: false, infos: [(2, [0, 2]), (2, [0, 2])]);
    TestOrth(FG.DiCyclic(6), indStab: false, infos: [(2, [0, 2]), (2, [0, 2])]);

    TestOrth(FG.SL2p(3));
    TestOrth(FG.GL2p(3));
    TestOrth(Group.SemiDirectProd(FG.SL2p(3), FG.AbelianMat(2)), infos: (2, [0, 1, 2, 3, 4, 5]));
    TestOrth(Product.Generate(FG.SL2p(3), FG.AbelianMat(2)), infos: (2, [0, 1, 2, 3, 4, 5]));

    GlobalStopWatch.Show(); // Time:54.696s
}

void KerChisGroup<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    var ct = FG.CharacterTable(g);
    if (!ct.TableComplete)
        ct.SolveOrthogonality((2, [0, 1, 2]));

    ct.DisplayCells();
    
    var listKers = ct.AllCharacters.Select((chi, k) => (chi, k))
        .Select(e => Group.Generate($"Ker(Ꭓ.{e.k + 1})", g, e.chi.Kernel().ToArray()))
        .ToList();

    var start = new HashSet<ConcreteGroup<T>>(listKers, new GroupSetEquality<T>());
    var allNormals = new HashSet<ConcreteGroup<T>>(listKers, new GroupSetEquality<T>());
    var sz = 0;
    while (sz != allNormals.Count)
    {
        sz = allNormals.Count;
        foreach (var n1 in allNormals.ToArray())
        {
            foreach (var n2 in start)
            {
                if (n2.SubSetOf(n1) || n1.SubSetOf(n2))
                    continue;

                var n3 = Group.Generate($"{n1} ∩ {n2}", g, n1.Intersect(n2).ToArray());
                allNormals.Add(n3);
            }
        }
    }

    var gAllSubgrs = g.AllSubgroups();
    Console.WriteLine(gAllSubgrs.Infos);
    allNormals.OrderBy(n => n.Count()).Println(e => e.ShortName, $"Total:{allNormals.Count}");
    Console.WriteLine();
    if (gAllSubgrs.Infos.AllNorms != allNormals.Count)
        throw new();
}

{
    // runOrth();
    
    foreach (var g in FG.AllGroupsOfOrder(24, 32))
        KerChisGroup(g);
    
}