using System.ComponentModel.DataAnnotations;
using FastGoat;
using FastGoat.Examples;
using FastGoat.Gp;
using FastGoat.UserGroup;
using static FastGoat.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    var c2 = new Cn(2);
    var c3 = new Cn(3);
    var c4 = new Cn(4);
    var c6 = new Cn(6);
    var s4 = new Symm(4);
    var s5 = new Symm(5);
    var a4 = Group.Generate("A4", s4[(1, 2, 3)], s4[(1, 2), (3, 4)]);
    var c3c4 = Group.SemiDirectProd(c3, c4);
    var d12 = Group.SemiDirectProd(c6, c2);
    var d12b = Group.Generate(s5, s5[(1, 4, 5)], s5[(2, 3)], s5[(4, 5)]);
    
    DisplayGroup.HeadElements(a4);
    DisplayGroup.HeadElements(d12b);
    DisplayGroup.HeadElementsSdp(c3c4);
    DisplayGroup.HeadElementsSdp(d12);

    Console.WriteLine(a4.IsIsomorphicTo(c3c4));
    Console.WriteLine(a4.IsIsomorphicTo(d12));
    Console.WriteLine(c3c4.IsIsomorphicTo(d12b));
    Console.WriteLine(d12.IsIsomorphicTo(d12b));
    
}