using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.EllCurve;
using static FastGoat.Commons.IntExt;
using RegXGroup = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
RecomputeAllPrimesUpTo(200000);

void testEllDB()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    GlobalStopWatch.Restart();

    var ellDB = EllipticExt.LMFDB_Ell_Q().Select(e => new EllDB(e.name, e.conductor, e.rank, e.torsType, e.model))
        .ToArray();
    foreach (var e in ellDB)
    {
        Console.WriteLine(e);
        var ellAnRank = EllipticCurvesPart2.EllAnalyticRank(e.model);
        var ngl = EC.NagellLutzTorsionGroup(ellAnRank.E.ToEllGroup());
        if (ellAnRank.rank != e.rank || ellAnRank.N.Num != e.conductor || !ngl.abType.SequenceEqual(e.torsType))
            throw new($"N={ellAnRank.N} rank={ellAnRank.rank} torsType=[{ngl.abType.Glue(", ")}]");
    }

    GlobalStopWatch.Show($"EllDB {ellDB.Length} curves");
    Console.WriteLine();
}