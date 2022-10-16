using FastGoat.Gp;
using FastGoat.UserGroup;

namespace FastGoat.Examples;

public static class GroupAction
{
    public static void C7_C3()
    {
        var c7 = new Cn(7);
        var c3 = new Cn(3);

        var u7 = new Un(7); // U7 = Aut(C7) 
        var y = u7[(c7[1], c7[2])];
        var pMap = Group.PartialMap((c3[1], y));
        var homoC3toU7 = Group.HomomorphismMap(c3, u7, pMap); // an element of Hom(C3, U7)

        Console.WriteLine("Action of C3 on C7");
        foreach (var kp in homoC3toU7)
            Console.WriteLine("g={0}, Y(g) = [{1}]", kp.Key, kp.Value);

        Console.WriteLine();
        Console.WriteLine("Result");
        var sdp = Group.SemiDirectProd(c7, c3);
        foreach (var n in c7)
        {
            foreach (var g in c3)
            {
                var n0 = homoC3toU7[g][n];
                var n1 = sdp.Act(n, g);
                Console.WriteLine("({0},{1}) -> ({2},{1}) ~ {3}", n, g, n0, n1);
            }
        }
    }

    public static void GroupOrder72()
    {
        var c2c4 = Product.Generate(new Cn(2), new Cn(4));
        var e1 = c2c4[1, 0];
        var e2 = c2c4[0, 1];
        var oe1 = c2c4.ElementsOrders[e1]; // Obvious, it is 2
        var oe2 = c2c4.ElementsOrders[e2];

        var c3c3 = Product.Generate(new Cn(3), new Cn(3));
        var autC3C3 = Group.AutomorphismGroup(c3c3);
        var byOrder = autC3C3.GroupBy(a => autC3C3.ElementsOrders[a]).ToDictionary(a => a.Key, b => b.ToArray());

        // Automorphism order divides generator order, it is a necessary condition to create 
        // a valid homomorphism but it is not a sufficient condition
        var orders1 = byOrder.Where(kp => kp.Key != 1 && oe1 % kp.Key == 0).SelectMany(kp => kp.Value).ToArray();
        var orders2 = byOrder.Where(kp => kp.Key != 1 && oe2 % kp.Key == 0).SelectMany(kp => kp.Value).ToArray();

        var allPossibles = orders1.SelectMany(a1 => orders2.Select(a2 => (a1, a2))).ToArray();
        var pMap = (Automorphism<Ep2<ZnInt, ZnInt>> a1, Automorphism<Ep2<ZnInt, ZnInt>> a2) =>
            Group.PartialMap((e1, a1), (e2, a2));

        var homC2C4ToAutC3C3 = allPossibles.Select(e => pMap(e.a1, e.a2))
            .Select(pm => Group.HomomorphismMap(c2c4, autC3C3, pm))
            .First(map => map.Count != 0);

        Console.WriteLine("(C3 x C3) :y (C2 x C4)");
        Console.WriteLine("y = Hom(C2 x C4, Aut(C3 x C3))");
        foreach (var kp in homC2C4ToAutC3C3)
        {
            Console.WriteLine("g={0} y(g) = [{1}]", kp.Key, kp.Value);
        }

        var y = homC2C4ToAutC3C3;
        var invert = (Ep2<Ep2<ZnInt, ZnInt>, Ep2<ZnInt, ZnInt>> x) =>
        {
            var n = x.E1;
            var g = x.E2;
            var ni = c3c3.Invert(n);
            var gi = c2c4.Invert(g);
            var gini = y[gi][ni];
            return Product.Elt(gini, gi);
        };

        var op = (Ep2<Ep2<ZnInt, ZnInt>, Ep2<ZnInt, ZnInt>> x1, Ep2<Ep2<ZnInt, ZnInt>, Ep2<ZnInt, ZnInt>> x2) =>
        {
            var n1 = x1.E1;
            var g1 = x1.E2;
            var n2 = x2.E1;
            var g2 = x2.E2;

            var yg1n2 = y[g1][n2];
            return Product.Elt(c3c3.Op(n1, yg1n2), c2c4.Op(g1, g2));
        };

        var group = Product.Group(c3c3, c2c4).ToHashSet();
        var cayleyTable = new Ep2<Ep2<ZnInt, ZnInt>, Ep2<ZnInt, ZnInt>>[group.Count, group.Count];
        bool isabelian = true;
        int i = 0;
        foreach (var a in group)
        {
            int j = 0;
            foreach (var b in group)
            {
                var ab = op(a, b);
                if (isabelian)
                {
                    var ba = op(b, a);
                    isabelian &= ab.Equals(ba);
                }

                cayleyTable[i, j] = ab;
                ++j;
            }

            ++i;
        }

        var lt = Enumerable.Range(0, group.Count).ToArray();
        var cayleyRows = lt.Select(i0 => lt.Select(j0 => cayleyTable[i0, j0]).ToHashSet()).All(group.SetEquals);
        var cayleyCols = lt.Select(i0 => lt.Select(j0 => cayleyTable[j0, i0]).ToHashSet()).All(group.SetEquals);

        if (!cayleyRows || !cayleyCols)
            throw new Exception();

        Console.WriteLine();
        Console.WriteLine("Test Group Passed with success !!!");
        Console.WriteLine("Group elements : {0}", group.Count);
        Console.WriteLine(isabelian ? GroupType.AbelianGroup : GroupType.NonAbelianGroup);
    }
}