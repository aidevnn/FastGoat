using System.Numerics;
using System.Reflection;
using System.Text;
using Craft;
using Craft.Craft;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.Tools;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Perm.Style = DisplayPerm.CyclesComma;

{
    var seq1 = UGCraft.HolCnMatrix(FG.ConcatPerm(FG.Cycles(2), FG.Cycles(2), FG.Cycles(4))).Order().ToArray();
    var seq2 = UGCraft.HolCn(FG.ConcatPerm(FG.Cycles(2), FG.Cycles(2), FG.Cycles(4))).Order().ToArray();
    // seq1.Println();
    Console.WriteLine($"seq:{seq1.Length} {seq1.SequenceEqual(seq2)}");
    Console.WriteLine();
}

{
    var seq1 = UGCraft.HolCnMatrix(FG.ConcatPerm(FG.Cycles(2), FG.Cycles(3), FG.Cycles(2), FG.Cycles(3), FG.Cycles(2))).Order().ToArray();
    var seq2 = UGCraft.HolCn(FG.ConcatPerm(FG.Cycles(2), FG.Cycles(3), FG.Cycles(2), FG.Cycles(3), FG.Cycles(2))).Order().ToArray();
    // seq1.Println();
    Console.WriteLine($"seq:{seq1.Length} {seq1.SequenceEqual(seq2)}");
    Console.WriteLine();
}

{
    var mt = FG.MetaCyclicPg(4, 4, 3);
    var autM = Group.AutomorphismGroup(mt);
    Console.WriteLine(mt.ShortName);
    Console.WriteLine(autM.ShortName);
    autM.GetGenerators().SelectMany(aut => UGCraft.InnerAutMatrix(aut)).ToArray().Println();
    Console.WriteLine();
}

{
    var mt = FG.MetaCyclicPg(8, 2, 3);
    var autM = Group.AutomorphismGroup(mt);
    Console.WriteLine(mt.ShortName);
    Console.WriteLine(autM.ShortName);
    autM.GetGenerators().SelectMany(aut => UGCraft.InnerAutMatrix(aut).ToArray()).Println();
    Console.WriteLine();
}

{
    var mt = FG.MetaCyclicPg(8, 2, 5);
    var autM = Group.AutomorphismGroup(mt);
    Console.WriteLine(mt.ShortName);
    Console.WriteLine(autM.ShortName);
    autM.GetGenerators().SelectMany(aut => UGCraft.InnerAutMatrix(aut)).ToArray().Println();
    Console.WriteLine();
}

{
    var mt = FG.MetaCyclicPg(7, 3, 2);
    var autM = Group.AutomorphismGroup(mt);
    Console.WriteLine(mt.ShortName);
    Console.WriteLine(autM.ShortName);
    autM.GetGenerators().SelectMany(aut => UGCraft.InnerAutMatrix(aut)).ToArray().Println();
    Console.WriteLine();
}

{
    var mt = FG.MetaCyclicPg(3, 8, 2);
    var autM = Group.AutomorphismGroup(mt);
    Console.WriteLine(mt.ShortName);
    Console.WriteLine(autM.ShortName);
    autM.GetGenerators().SelectMany(aut => UGCraft.InnerAutMatrix(aut)).ToArray().Println();
    Console.WriteLine();
}

{
    var mt = FG.MetaCyclicPg(6, 4, 5);
    var autM = Group.AutomorphismGroup(mt);
    Console.WriteLine(mt.ShortName);
    Console.WriteLine(autM.ShortName);
    autM.GetGenerators().SelectMany(aut => UGCraft.InnerAutMatrix(aut)).ToArray().Println();
    Console.WriteLine();
}

{
    var mt = FG.MetaCyclicPg(12, 2, 11);
    var autM = Group.AutomorphismGroup(mt);
    Console.WriteLine(mt.ShortName);
    Console.WriteLine(autM.ShortName);
    autM.GetGenerators().SelectMany(aut => UGCraft.InnerAutMatrix(aut)).ToArray().Println();
    Console.WriteLine();
}