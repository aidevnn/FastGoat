using System.Collections;
using System.Diagnostics;
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
    var gl = new GL(2, 3);
    var a = gl[2, 2, 2, 1];
    var b = gl[0, 2, 1, 0];
    var q8 = Group.Generate("Q8", gl, a, b);
    DisplayGroup.HeadElements(q8);

    var wgQ8 = new WordGroup("wgQ8", "a2=b2, a-1=bab-1");
    Console.WriteLine(wgQ8.IsIsomorphicTo(q8));
}