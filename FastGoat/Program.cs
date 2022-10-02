using FastGoat;
using FastGoat.Gp;
using FastGoat.UserGroup;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
{
    var s4 = new Sn(4);
    var S4 = Group.Generate("S4", s4[(1, 2)], s4[(1, 2, 3, 4)]);
    var V = Group.Generate("V", S4, s4[(1, 2), (3, 4)], s4[(1, 3), (2, 4)]);
    var Q = S4.Over(V);
    DisplayGroup.HeadClassesTable(Q);
}