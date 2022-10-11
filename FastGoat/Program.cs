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
    var v = Product.Generate(c2, c2);
    var c6 = new Cn(6);
    var s4 = new Symm(4);
    var s5 = new Symm(5);
    var s6 = new Symm(6);
    var a4 = Group.Generate("A4", s4[(1, 2, 3)], s4[(1, 2), (3, 4)]);
    var c3c4 = Group.SemiDirectProd(c3, c4);
    var c3c4b = Group.SemiDirectProd(c3, v);
    var d12 = Group.SemiDirectProd(c6, c2);
    var d12b = Group.Generate(s5, s5[(1, 4, 5)], s5[(2, 3)], s5[(4, 5)]);
    var d12c = Group.Generate(s6, s6[(1, 2, 3, 4, 5, 6)], s6[(1, 6), (2, 5), (3, 4)]);

    DisplayGroup.HeadElements(a4);
    DisplayGroup.HeadElementsSdp(c3c4);
    DisplayGroup.HeadElementsSdp(c3c4b); // No sense C3 : (C2 x C2)
    DisplayGroup.HeadElementsSdp(d12);
    DisplayGroup.HeadElements(d12b);
    DisplayGroup.HeadElements(d12c);

    Console.WriteLine(a4.IsIsomorphicTo(c3c4));
    Console.WriteLine(a4.IsIsomorphicTo(d12));
    Console.WriteLine(c3c4.IsIsomorphicTo(d12b));
    Console.WriteLine(c3c4b.IsIsomorphicTo(d12b)); // No sense C3 : (C2 x C2)
    Console.WriteLine(d12.IsIsomorphicTo(d12b));
    Console.WriteLine(d12c.IsIsomorphicTo(d12b));

    var allIsos = Group.AllIsomorphisms(d12b, d12c);
    Console.WriteLine();
    Console.WriteLine("Isomorsphims between D12 from Symm5 and Symm6. Total : {0}", allIsos.Count);
    foreach (var iso in allIsos)
    {
        foreach (var kp in iso.OrderBy(kp => d12b.ElementsOrders[kp.Key]))
        {
            Console.WriteLine($"{kp.Key,20} => {kp.Value}");
        }

        Console.WriteLine();
    }
}