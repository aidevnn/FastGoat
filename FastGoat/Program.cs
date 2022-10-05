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
    var g = Group.Create(Product.Group(new Cn(2), new Cn(2)));
    Console.WriteLine(g.LongestCycles.Keys.Glue(", "));
    var gAut = Group.Aut(g);
    var a0 = gAut[(g[1, 0], g[1, 0]), (g[0, 1], g[1, 1])];
    var a1 = gAut[(g[1, 0], g[0, 1]), (g[0, 1], g[1, 0])];
    var a2 = gAut[(g[1, 0], g[1, 1]), (g[1, 1], g[1, 0])];
    Console.WriteLine(gAut);
    Console.WriteLine(gAut.Neutral());
    Console.WriteLine();
    
    Console.WriteLine(a0);
    Console.WriteLine(a1);
    Console.WriteLine();
    
    Console.WriteLine(gAut.Invert(a0));
    Console.WriteLine(gAut.Invert(a1));
    Console.WriteLine();

    Console.WriteLine(gAut.Op(a0, a1));

    var ga0 = Group.Generate(a0);
    DisplayGroup.HeadElements(ga0);
    var ga1 = Group.Generate(a1);
    DisplayGroup.HeadElements(ga1);
    var ga01 = Group.Generate(a0, a1);
    DisplayGroup.HeadElements(ga01);
    Console.WriteLine(ga01.Contains(a2));
    var a3 = ga01[(g[1, 1], g[0, 1]), (g[0, 1], g[1, 1])];
    Console.WriteLine(a3.Equals(a0));
    Console.WriteLine(a3.Equals(a1));
}