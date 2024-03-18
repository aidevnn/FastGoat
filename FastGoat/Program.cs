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

void TestLoggerOffExts()
{
    var exts = FG.AllExtensions(
            (FG.Abelian(1), FG.Abelian(8)),
            (FG.Abelian(8), FG.Abelian(1)),
            (FG.Abelian(2), FG.Abelian(2, 2))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .FilterIsomorphic()
        .Naming()
        .DisplayNames(rename: true);
}

void TestLoggerLvl1Exts()
{
    Logger.Level = LogLevel.Level1;
    var exts = FG.AllExtensions(
            (FG.Abelian(1), FG.Abelian(8)),
            (FG.Abelian(8), FG.Abelian(1)),
            (FG.Abelian(2), FG.Abelian(2, 2))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .FilterIsomorphic()
        .Naming()
        .DisplayNames(rename: true);
}

void TestLoggerOffWG()
{
    DisplayGroup.HeadElements(FG.DihedralWg(5));
}

void TestLoggerLvl2WG()
{
    Logger.Level = LogLevel.Level2;
    DisplayGroup.HeadElements(FG.DihedralWg(5));
}

{
    // TestLoggerOffExts();
    // TestLoggerLvl1Exts();
    // TestLoggerOffWG();
    // TestLoggerLvl2WG();
}
