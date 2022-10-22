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
    var gr = new WordGroup("Q32", "a16, b2=a8, bab-1=a-1");
    var chain = new ZentrumChain<Word>(gr);
}