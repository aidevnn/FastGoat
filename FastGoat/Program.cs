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
    var a = gl[2, 0, 0, 1];
    var b = gl[2, 1, 2, 0];

    var g1 = Group.Generate(gl, a, b);
    DisplayGroup.HeadElements(g1); // 48 elements
}

{
    var gl = new GL(2, 3);
    var a = gl[1, 1, 0, 1];
    var b = gl[0, 1, 2, 0];
    
    var g1 = Group.Generate(gl, a, b);
    DisplayGroup.Head(g1); // 24 elements
}

{
    var gl = new GL(3, 2);
    var a = gl.At((0, 1, 4, 8), (1, 1, 1, 1));
    var b = gl.At((2, 3, 7), (1, 1, 1));
    var g1 = Group.Generate(gl, a, b);
    DisplayGroup.Head(g1); // 168 elements
}