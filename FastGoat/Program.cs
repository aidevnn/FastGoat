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
    // gap> StructureDescription(SmallGroup(18,4));
    // "(C3 x C3) : C2"
    // gap> SortedList(List(SmallGroup(18,4),Order));
    // [ 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3 ]

    var n = Product.Generate(new Cn(3), new Cn(3));
    var g = new Cn(2);
    var thetas = Group.AllOpsByAutomorphisms(g, n);
    foreach (var theta in thetas)
    {
        var g1 = new SemiDirectProduct<Ep2<ZnInt, ZnInt>, ZnInt>("G", n, theta, g);
        DisplayGroup.HeadSdp(g1);
        Console.WriteLine("Sorted Orders ({0})", g1.ElementsOrders.Values.Ascending().Glue(", "));

        // Quotient by normal group
        DisplayGroup.Head(g1.Over(g1.Ncan)); 

        // // Comparing theta(g0)(n0) with g0.n0.g0^1
        foreach (var g0 in g)
        {
            var g01 = Product.Elt(n.Neutral(), g0);
            foreach (var n0 in n)
            {
                var g0n0 = theta[g0][n0];
                var n01 = Product.Elt(n0, g.Neutral());
                var g0n0g0i = g1.Op(g01, g1.Op(n01, g1.Invert(g01)));
                if (!g0n0.Equals(g0n0g0i.E1) || !g0n0g0i.E2.Equals(g.Neutral()))
                    throw new Exception();
            }
        }

        // Split extension 1 -> N -> N:G -> G -> 1
        // with i = N -> N:G and p = N:G -> G
        var mapI = Group.PartialMap((n[1, 0], g1[(1, 0), 0]), (n[0, 1], g1[(0, 1), 0]));
        var homI = Group.HomomorphismMap(n, g1, mapI);
        var mapP = Group.PartialMap((g1[(0, 1), 0], g[0]), (g1[(0, 1), 1], g[1]), (g1[(1, 0), 0], g[0]),
            (g1[(1, 0), 1], g[1]));
        var homP = Group.HomomorphismMap(g1, g, mapP);
        Console.WriteLine("i = N -> N:G :({0}", homI.GlueMap());
        Console.WriteLine("p = N:G -> G : ({0}", homP.GlueMap());
        var imI = homI.Values.Ascending().ToArray();
        var kerP = homP.Where(kp => kp.Value.Equals(g.Neutral())).Select(kp => kp.Key).Ascending().ToArray();
        Console.WriteLine("Im(i)  = {0}", imI.Glue(" "));
        Console.WriteLine("Ker(p) = {0}", kerP.Glue(" "));
        Console.WriteLine("Im(i)  = Ker(p) : {0}",imI.SequenceEqual(kerP));

        Console.WriteLine("#########################################################################");
    }
}