using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;

namespace FastGoat.Examples;

public static class PSL2q
{
    static void L2pCompute(int p)
    {
        var gl = new GL(2, p);
        var x = (p - 2).Range(2).First(i => IntExt.PowMod(i, p - 1, p) == 1);
        var xp = IntExt.PowMod(x, p - 2, p);

        // SL(2,q) generators from H.E. Rose, page 271, Problem 12.5
        var a = gl[x, 0, 0, xp];
        var b = gl[p - 1, 1, p - 1, 0];

        var SL2q = Group.Generate($"SL(2,{p})", gl, a, b);
        var z = Group.Zentrum(SL2q);
        var L2q = SL2q.Over(z, $"L2({p})");
        DisplayGroup.Head(L2q);
    }

    static void L2qCompute(int q)
    {
        var gl = new GLnq(2, q);
        if (gl.Fq.M == 1)
            return;
        var x = gl.Fq['x'];

        // SL(2,q) generators from H.E. Rose, page 271, Problem 12.5
        var a = gl[x, 0, 0, x.Pow(q - 2)];
        var b = gl[q - 1, 1, q - 1, 0];

        var SL2q = Group.Generate($"SL(2,{q})", gl, a, b);
        var z = Group.Zentrum(SL2q);
        var L2q = SL2q.Over(z, $"L2({q})");
        DisplayGroup.Head(L2q);
    }

    public static void L2p()
    {
        foreach (var p in IntExt.Primes10000.Where(p => p is > 3 and < 54))
            L2pCompute(p);
    }

    public static void L2q()
    {
        foreach (var q in new[] { 4, 8, 9, 16, 25, 27, 32, 49 })
            L2qCompute(q);
    }

    public static void L2_4()
    {
        var gl = new GLnq(2, 4);
        var a = gl['x', 0, 0, 1];
        var b = gl[1, 1, 1, 0];

        var gl2_4 = Group.Generate(gl, a, b);
        DisplayGroup.HeadOrders(gl2_4);

        var one = gl.Fq.One;
        var detPositive = gl2_4.Where(e => gl.Determinant(e).Equals(one)).ToArray();
        Console.WriteLine("SL count :{0}", detPositive.Length);
        var a1 = gl['x', 0, 0, ('x', 2)];
        var sl2_4 = Group.Generate("SL2(4)", gl, a1, b);
        DisplayGroup.HeadOrders(sl2_4);

        var zg2_4 = Group.Zentrum(sl2_4);
        DisplayGroup.Head(zg2_4);
        var l2_4 = sl2_4.Over(zg2_4, "L2(4)");
        DisplayGroup.Head(l2_4);

        DisplayGroup.AreIsomorphics(l2_4, FG.Alternate(5));
    }

    public static void L2_8()
    {
        var gl = new GLnq(2, 8);
        var a = gl['x', 0, 0, 1];
        var b = gl[1, 1, 1, 0];

        var gl2_8 = Group.Generate(gl, a, b);
        DisplayGroup.HeadOrders(gl2_8); // |Gl(2, F8)| = 3528

        var one = gl.Fq.One;
        var detPositive = gl2_8.Where(e => gl.Determinant(e).Equals(one)).ToArray();
        Console.WriteLine("SL count :{0}", detPositive.Length);
        var a1 = gl['x', 0, 0, ('x', 6)];
        var sl2_8 = Group.Generate("SL2(8)", gl, a1, b);
        DisplayGroup.HeadOrders(sl2_8);

        var zg2_8 = Group.Zentrum(sl2_8);
        DisplayGroup.Head(zg2_8);
        var l2_8 = sl2_8.Over(zg2_8, "L2(8)");
        DisplayGroup.Head(l2_8);
    }

    public static void L2_9()
    {
        var gl = new GLnq(2, 9);
        var a = gl['x', 0, 0, 1];
        var b = gl[2, 1, 2, 0];

        var gl2_9 = Group.Generate(gl, a, b);
        DisplayGroup.HeadOrders(gl2_9); // |Gl(2, F9)| = 5760

        var a1 = gl['x', 0, 0, ('x', 7)];
        var sl2_9 = Group.Generate("SL2(9)", gl, a1, b);
        DisplayGroup.HeadOrders(sl2_9); // |SL2(9)| = 720

        var s6 = new Symm(6);
        DisplayGroup.HeadOrders(s6);
        DisplayGroup.AreIsomorphics(s6, sl2_9);

        var zg2_9 = Group.Zentrum(sl2_9);
        var l2_9 = sl2_9.Over(zg2_9, "L2(9)");
        DisplayGroup.HeadOrders(l2_9);

        DisplayGroup.AreIsomorphics(l2_9, FG.Alternate(6));
    }

    public static void L2_16()
    {
        var gl = new GLnq(2, 16);
        var a = gl['x', 0, 0, 1];
        var b = gl[1, 1, 1, 0];

        // var gl216 = Group.Generate(gl, a, b);
        // DisplayGroup.HeadOrders(gl216); // |Gl(2, F16)| = 61200

        var a1 = gl['x', 0, 0, ('x', 14)];
        var sl2_16 = Group.Generate("SL2(16)", gl, a1, b);
        DisplayGroup.HeadOrders(sl2_16); // |SL2(16)| = 4080

        var zg2_16 = Group.Zentrum(sl2_16);
        DisplayGroup.Head(zg2_16);
        var l2_16 = sl2_16.Over(zg2_16, "L2(16)");
        DisplayGroup.Head(l2_16); // |L2(16)| = 4080
    }

    public static void L2_25()
    {
        var gl = new GLnq(2, 25);

        var a = gl['x', 0, 0, ('x', 23)];
        var b = gl[4, 1, 4, 0];
        var sl2_25 = Group.Generate("SL2(25)", gl, a, b);
        DisplayGroup.HeadOrders(sl2_25); // |SL2(25)| = 15600


        var zg2_25 = Group.Zentrum(sl2_25);
        DisplayGroup.Head(zg2_25); // |Z(SL2(25))| = 2
        var l2_25 = sl2_25.Over(zg2_25, "L2(25)");
        DisplayGroup.Head(l2_25); // |L2(25)| = 7800
    }

    public static void L2_27()
    {
        var gl = new GLnq(2, 27);

        var a = gl['x', 0, 0, ('x', 25)];
        var b = gl[2, 1, 2, 0];
        var sl2_27 = Group.Generate("SL2(27)", gl, a, b);
        DisplayGroup.HeadOrders(sl2_27); // |SL2(27)| = 19656

        var zg2_27 = Group.Zentrum(sl2_27);
        DisplayGroup.Head(zg2_27); // |Z(SL2(27))| = 2
        var l2_27 = sl2_27.Over(zg2_27, "L2(27)");
        DisplayGroup.Head(l2_27); // |L2(27)| = 9828
    }

    public static void L2pWordGroup()
    {
        foreach (var p in new[] { 5, 7, 11, 13 })
        {
            if (!IntExt.Primes10000.Contains(p) || p <= 3)
                throw new();

            var s = $"a4ba{(p + 1) / 2}b";
            var rel = $"a{p}, {s}{s}, ababab=b2";
            var wg = FG.WordGroup($"L({p})", rel);

            Console.WriteLine($"{wg.ShortName}");
            Console.WriteLine($"{wg.Definition}");
            DisplayGroup.AreIsomorphics(wg, FG.L2p(p));
            Console.WriteLine();
        }
    }
}