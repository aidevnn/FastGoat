using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
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

void NamingExtsExamples()
{
    FG.AllExtensions((FG.Abelian(2, 2), FG.Abelian(3)), (FG.Abelian(6), FG.Abelian(2))).NamingExts().DisplayExts();
    FG.AllExtensions(
            (FG.Abelian(2, 2, 2), FG.Abelian(2)),
            (FG.Abelian(4, 2), FG.Abelian(2)),
            (FG.Abelian(8), FG.Abelian(2)))
        .NamingExts()
        .DisplayExts();

    FG.AllExtensions((FG.Abelian(5), FG.Abelian(2, 2)), (FG.Abelian(5), FG.Abelian(4))).NamingExts().DisplayExts();
    FG.AllExtensions((FG.Abelian(2, 5), FG.Abelian(2, 2)), (FG.Abelian(2, 5), FG.Abelian(4))).NamingExts().DisplayExts();
    FG.AllExtensions((FG.Abelian(4, 5), FG.Abelian(2, 2)), (FG.Abelian(4, 5), FG.Abelian(4))).NamingExts().DisplayExts();
}

// void Ord64()
{
    var nbOpsMax = 10000;
    var s32 = new Sn(32);
    var g1 = s32[(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16), (17, 18, 19, 20), (21, 22, 23, 24),
        (25, 26, 27, 28), (29, 30, 31, 32)];
    var g2 = s32[(1, 17, 15, 6), (2, 18, 16, 7), (3, 19, 13, 8), (4, 20, 14, 5), (9, 29, 24, 28), (10, 30, 21, 25),
        (11, 31, 22, 26), (12, 32, 23, 27)];
    var g3 = s32[(1, 26, 15, 31), (2, 32, 16, 27), (3, 28, 13, 29), (4, 30, 14, 25), (5, 12, 20, 23),
        (6, 24, 17, 9), (7, 10, 18, 21), (8, 22, 19, 11)];

    var sm3232 = Group.Generate("[(C4 x C4) . C2]pg", s32, g1, g2, g3).ToCGTable();

    var c2 = FG.Abelian(2);
    var sdps1 = Group.AllSemiDirectProd(FG.Abelian(8), FG.Abelian(4)).Select(e => e.ToCGTable()).Naming().ToArray();
    var sdps2 = Group.AllSemiDirectProd(FG.Abelian(4), FG.Abelian(8)).Select(e => e.ToCGTable()).Naming().ToArray();
    var sdp3 =  Product.Generate(c2, FG.ModularMaxSdp(4)).ToCGTable() ;
    var sdp4 =  Group.AllSemiDirectProd(FG.Abelian(2, 4), FG.Abelian(4)).Select(e => e.ToCGTable()).Naming().ToArray()[2];

    var e245 = (10, 30, sm3232, c2); // Group64 Id 172,245
    var tuplesC2 = sdps1.Concat([..sdps2, sdp3, sdp4]).Select(e => (nbOpsMax, 0, e, c2)).ToArray(); // Group64 Id 11,13,14,22,79,81,82,160,172,180
    var tuplesC4 = new[] { 
        (nbOpsMax, 0, FG.Abelian(8, 2).ToCGTable(), FG.Abelian(4)), 
        (nbOpsMax, 0, FG.Abelian(4, 4).ToCGTable(), FG.Abelian(4)) }; // Group64 Id 19,22,37,45
    var tuplesC2C2 = new[] { (nbOpsMax, 0, FG.Abelian(2, 8).ToCGTable(), FG.Abelian(2, 2)) }; // Group64 Id 49,45,180,160,168,43,172 // Minor Bugs
    
    GlobalStopWatch.Restart();
    FG.AllExtensions([..tuplesC2, ..tuplesC4, ..tuplesC2C2, e245]).NamingExts().DisplayExts();
    GlobalStopWatch.Show("End");
    Console.Beep();
    // Total Exts:192
    // # End Time:2770743 ms ~ 46 min
    // Need ~2048 KB for console size
}
