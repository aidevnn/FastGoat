using System.Collections;
using FastGoat;
using FastGoat.Examples;
using FastGoat.Gp;
using FastGoat.ToddCoxeter;
using FastGoat.UserGroup;
using static FastGoat.IntExt;
using static FastGoat.EnumerableExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    // Algebre Tome 1, Daniel Guin – Thomas Hausberger
    ToddCoxeterAlgo.Run("a", "a4, b3, abab", details: true); // p101 step by step algorithm with subgroup H=<a>
    ToddCoxeterAlgo.Run("a", "a3, b3, abab", details: true); // p106 step by step algorithm with subgroup H=<a>
    ToddCoxeterAlgo.Run("a", "a3, b3, aba2b", details: true); // p107 step by step algorithm with subgroup H=<a>
    
    ToddCoxeterAlgo.Run("b", "a7, b3, a2=bab-1", details: true); // C7 : C3 step by step algorithm with subgroup H=<b>
    ToddCoxeterAlgo.Run("a7, b3, a2=bab-1", details: true); // C7 : C3 step by step algorithm with subgroup H=<Id>
    
    ToddCoxeterAlgo.Run("a2, b4, ab=ba", details: true); // Dihedral 8 step by step algorithm
    ToddCoxeterAlgo.Run("a4, a2=b2, b=aba").DisplayOps(); // Quartenion Table
    ToddCoxeterAlgo.Run("a4, b4, c2=b2, ab=ba, cac-1=ab2, cbc-1=a2b").DisplayOps(); // (C4 x C4) . C2 Table
    ToddCoxeterAlgo.Run("a, b", "a4, b4, c2=b2, ab=ba, cac-1=ab2, cbc-1=a2b", details: true); // C2 Table with subgroup H = C4 x C4=<a,b>
}