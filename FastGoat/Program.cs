using FastGoat;
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
    var g1 = Group.SemiDirectProd(Product.Generate(new Cn(3), new Cn(3)), new Cn(2));
    DisplayGroup.HeadElementsSdp(g1);

    var g2 = Group.SemiDirectProd(new Cn(9), new Cn(2));
    DisplayGroup.HeadElementsSdp(g2);

    var g3 = Group.SemiDirectProd(new Cn(8), Product.Generate(new Cn(2), new Cn(2)));
    DisplayGroup.HeadElementsSdp(g3);

    var g4 = Group.SemiDirectProd(Product.Generate(new Cn(2), new Cn(2), new Cn(2)), Product.Generate(new Cn(2), new Cn(2)));
    DisplayGroup.HeadElementsSdp(g4);

    var g5 = Group.SemiDirectProd(new Cn(3), new Cn(8));
    var g6 = Group.SemiDirectProd(g5, new Cn(2));
    DisplayGroup.HeadSdp(g6);
}