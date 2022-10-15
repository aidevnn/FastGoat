using System.CodeDom;
using System.Threading.Channels;
using FastGoat;
using FastGoat.Examples;
using FastGoat.Gp;
using FastGoat.UserGroup;
using static FastGoat.IntExt;
using static FastGoat.EnumerableExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
//
// {
//     GlobalStopWatch.Restart();
//     var n = 5;
//     var sn = new Symm(n);
//     var Xn = new Cn(n);
//
//     for (int k = 1; k <= n; k++)
//     {
//         var Xnk = Product.Gp(Enumerable.Repeat(Xn, k).ToArray());
//         var setXnk = Xnk.Where(xnk => xnk.Ei.Distinct().Count() == Xnk.Gi.Length).ToArray();
//         Console.WriteLine(setXnk.Length);
//         Ep<ZnInt> Image(Perm g, Ep<ZnInt> x) => Product.Ep(x.Ei.Select(e => Xn[g.Table[e.K]]).ToArray());
//         foreach (var gens in sn.SubGroupsGenerators())
//         {
//             var g1 = Group.Generate(sn, gens.ToArray());
//             var allOrbits = Group.AllOrbits(g1, setXnk, Image);
//             if (allOrbits.Count == 1)
//             {
//                 DisplayGroup.Head(g1);
//                 Group.DisplayOrbx(allOrbits);
//                 Console.WriteLine();
//             }
//         }
//
//         GlobalStopWatch.Show($"X{n}{k}");
//     }
//
//     GlobalStopWatch.Show($"{sn}");
// }
//
// {
// }

{
    var s4 = new Symm(4);
    Group.DisplayOrbx(s4, Group.ByConjugate(s4));
    Console.WriteLine();

    var c3c4 = Product.Generate(new Cn(3), new Cn(4));
    Group.DisplayOrbx(c3c4, Group.ByConjugate(c3c4));
    Console.WriteLine();

    var grG = new Symm(4);
    var grH = Group.Generate("H", grG, grG[(1, 2, 3)]);
    var cs = Group.Cosets(grG, grH, Coset.Left);
    
    // transitivity
    Group.DisplayOrbx(grG, cs.Values.Distinct().ToArray(), Group.ByLeftCoset(grG, grH)); 
}