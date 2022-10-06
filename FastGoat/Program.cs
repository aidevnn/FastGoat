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
    var sn = new Sn(3);
    var autSn = Group.Aut(sn, sn.GetGenerators().ToArray());
    DisplayGroup.HeadElements(autSn);

    var v = Product.Group(new Zn(2), new Zn(2));
    var autV = Group.Aut("V", v, v.GetGenerators().ToArray());
    DisplayGroup.HeadElements(autV);

    var k8 = Product.Group(new Zn(2), new Zn(2), new Zn(2));
    var autK8 = Group.Aut("K8", k8, k8.GetGenerators().ToArray());
    DisplayGroup.Head(autK8);
    
    var gZ3Z3=Product.Group(new Zn(3), new Zn(3));
    var autZ3Z3 = Group.Aut(gZ3Z3, gZ3Z3.GetGenerators().ToArray());
    DisplayGroup.Head(autZ3Z3);
}