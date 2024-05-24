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

void TestOrth<T>(ConcreteGroup<T> g, bool indStab = true, params (int dim, int[] linIdx)[] infos) where T : struct, IElt<T>
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

{
    GlobalStopWatch.Restart();
    TestOrth(FG.SL2p(3), infos: (2, [0, 1, 2]));
    TestOrth(FG.GL2p(3), infos: (2, [0, 1]));
    TestOrth(Group.SemiDirectProd(FG.SL2p(3), FG.AbelianMat(2)), infos: (2, [0, 1, 2, 3, 4, 5]));
    TestOrth(Product.Generate(FG.SL2p(3), FG.AbelianMat(2)), infos: (2, [0, 1, 2, 3, 4, 5]));

    TestOrth(FG.Quaternion(16), indStab: false, infos: (2, [0, 1]));
    TestOrth(FG.MetaCyclicSdp(4, 4, 3), indStab: false, infos: (2, [0, 4]));
    TestOrth(FG.DiCyclic(5), indStab: false, infos: [(2, [0, 2]), (2, [0, 2])]);
    TestOrth(FG.DiCyclic(6), indStab: false, infos: [(2, [0, 2]), (2, [0, 2])]);
    GlobalStopWatch.Show(); // Time:54.696s
}