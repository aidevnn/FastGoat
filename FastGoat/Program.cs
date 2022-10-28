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
    var cn = new Cn(40);
    var h = Group.Generate("H", cn, cn[5]);
    DisplayGroup.HeadCosets(cn.Over(h), details: true);
}