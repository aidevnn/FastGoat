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

void TestSL(ConcreteGroup<Mat> g, ConcreteGroup<Mat> super)
{
    DisplayGroup.HeadNames(g);
    
    GlobalStopWatch.AddLap();
    var ct0 = FG.CharacterTableEmpty(g);
    ct0.DerivedSubGroupLift();
    ct0.InductionFromStabilizers();
    if (g.Count() == 24)
        ct0.SolveOrthogonality((2, 3.Range()));
    else
        ct0.SolveOrthogonality((2, 6.Range()));

    ct0.DisplayCells(tableOnly: true);
    GlobalStopWatch.Show("With orthogonality");
    Console.WriteLine();

    GlobalStopWatch.AddLap();
    var ct2 = FG.CharacterTable(super);
    var ct1 = FG.CharacterTableEmpty(g);
    ct1.DerivedSubGroupLift();
    ct1.InductionFromStabilizers();
    ct1.RestrictionFromSuperGroup(ct2);
    ct1.DisplayCells(tableOnly: true);
    GlobalStopWatch.Show("With restriction from super group");
    Console.WriteLine();

    ct0.AllCharacters.Zip(ct1.AllCharacters)
        .Select(e => (e.First - e.Second).Simplify())
        .Where(chi => !chi.IsZero())
        .Println("Diff tables"); // empty
}

{
    // Console.WriteLine(48.Range(1).Sum(i => GroupExt.A000001[i])); // 250 groups up to order 48
    
    // Generators of SL(2,3) in GL(4,5)
    // gen1 of order 3
    // [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0]
    // gen2 of order 4
    // [0, 0, 4, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 4, 0, 0]
    //
    // [4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0]
    // 
    
    // Generators of SL(2,3) x: C2 in GL(4,5)
    // gen1 of order 2
    // [0, 0, 2, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 2, 0, 0]
    // gen2 of order 3
    // [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0]
    //
    // [4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0]
    // 
    
    // Generators of C2 x SL(2,3) in GL(5,5)
    // gen1 of order 2
    // [4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1]
    // gen2 of order 3
    // [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0]
    // gen3 of order 4
    // [1, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 4, 0, 0]
    //
    // [1, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]


    var gl45 = new GL(4, 5);
    var a0 = gl45[0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0];
    var a1 = gl45[0, 0, 4, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 4, 0, 0];
    var a2 = gl45[4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0];
    var a3 = gl45[0, 0, 2, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 2, 0, 0];
    var a4 = gl45[0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0];
    var sl23 = Group.Generate("SL(2,3)", gl45, a0, a1);
    var gl23 = Group.Generate("GL(2,3)", gl45, a0, a1, a2);
    var sl23byc2 = Group.Generate("SL(2,3) x: C2", gl45, a3, a4);
    var gl23byc2 = Group.Generate("GL(2,3) x: C2", gl45, a3, a0, a2);

    var gl55 = new GL(5, 5);
    var b0 = gl55[4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
    var b1 = gl55[1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0];
    var b2 = gl55[1, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 4, 0, 0];
    var b3 = gl55[1, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0];
    
    var c2sl23 = Group.Generate("C2 x SL(2,3)", gl55, b0, b1, b2);
    var c2gl23 = Group.Generate("C2 x GL(2,3)", gl55, b0, b1, b2, b3);

    GlobalStopWatch.Restart();
    Logger.Level = LogLevel.Level1;

    TestSL(sl23, gl23);
    TestSL(sl23byc2, gl23byc2);
    TestSL(c2sl23, c2gl23);
}
