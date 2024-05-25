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
    if (g.GroupType == GroupType.AbelianGroup)
        return;
    
    var ct = FG.CharacterTableEmpty(g);
    ct.DerivedSubGroupLift();
    ct.InductionFromStabilizers();
    if(ct.TableComplete)
        return;
    
    ct.DisplayCells();

    var gr = Group.MulGroup("LinChi", ct.DoneChis.Where(chi => chi.IsLinear).ToArray());
    DisplayGroup.HeadElementsNames(gr);
    DisplayGroup.Generators(gr);
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

{
    // runOrth();
    foreach (var g in FG.AllGroupsOfOrder(1, 63))
        LinChisMulGroup(g);
}