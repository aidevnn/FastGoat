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

        var c3c3 = Product.Group(new Cn(3), new Cn(3));
        var autC3C3 = Group.Aut(c3c3[1, 0], c3c3[0, 1]);
        var byOrder = autC3C3.GroupBy(a => autC3C3.ElementsOrders[a]).ToDictionary(a => a.Key, b => b.ToArray());
        var allPossibles = byOrder[oe1].SelectMany(a1 => byOrder[oe2].Select(a2 => (a1, a2))).ToArray();
        var pMap = (Automorphism<Ep2<ZnInt, ZnInt>> a1, Automorphism<Ep2<ZnInt, ZnInt>> a2) =>
            Group.PartialMap((e1, a1), (e2, a2));

        Console.WriteLine(allPossibles.Length);
        var solution = allPossibles.Select(e => pMap(e.a1, e.a2)).Select(pm => Group.HomomorphismMap(c2c4, autC3C3, pm))
            .First(map => map.Count != 0);
        
        var homC2C4ToC3C3 = solution;
        foreach (var kp in homC2C4ToC3C3)
        {
            Console.WriteLine("g={0} y(g) = [{1}]", kp.Key,kp.Value);
        }
        
        foreach (var n in c3c3)
        {
            foreach (var g in c2c4)
            {
                var n0 = homC2C4ToC3C3[g][n];
                Console.WriteLine("({0},{1}) -> ({2},{1})", n, g, n0);
            }
        }
    }
}