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

{
    GlobalStopWatch.Restart();
    var n = 5;
    var sn = new Symm(n);
    var Xn = new Cn(n);

    for (int k = 1; k <= n; k++)
    {
        var Xnk = Product.Gp(Enumerable.Repeat(Xn, k).ToArray());
        var setXnk = Xnk.Where(xnk => xnk.Ei.Distinct().Count() == Xnk.Gi.Length).ToArray();
        Console.WriteLine(setXnk.Length);
        Ep<ZnInt> Image(Perm g, Ep<ZnInt> x) => Product.Ep(x.Ei.Select(e => Xn[g.Table[e.K]]).ToArray());
        foreach (var gens in sn.SubGroupsGenerators())
        {
            var g1 = Group.Generate(sn, gens.ToArray());
            // DisplayGroup.Head(g1);
            Group.AllOrbits(g1, setXnk, Image);
            // Console.WriteLine();
        }
        GlobalStopWatch.Show($"X{n}{k}");
    }
    
    GlobalStopWatch.Show($"{sn}");
}
