using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;

namespace FastGoat.Examples;

// H.E.Rose, A Course on Finite Groups, page 266 
public static class UnitaryGroup
{
    public static void U3_2()
    {
        var gl34 = new GLnq(3, 4);
        var x = gl34.Fq['x'];
        // gap> GeneratorsOfGroup(SU(3,2));
        var a1 = gl34[
            1, x, x,
            0, 1, x + 1,
            0, 0, 1
        ];
        var b1 = gl34[
            x, 1, 1,
            1, 1, 0,
            1, 0, 0
        ];

        var su3_2 = Group.Generate("SU3(2)", gl34, a1, b1);
        DisplayGroup.Head(su3_2);

        var z = Group.Zentrum(su3_2);
        DisplayGroup.Head(z);
        var u3_2 = su3_2.Over(z, "U3(2)");
        DisplayGroup.Head(u3_2);

        var u3_2sdp = Group.SemiDirectProd(FG.Abelian(3, 3), FG.Quaternion(8));
        DisplayGroup.Head(u3_2sdp);
        DisplayGroup.AreIsomorphics(u3_2, u3_2sdp);
    }

    public static void U3_3()
    {
        var gl39 = new GLnq(3, 9);
        var x = gl39.Fq['x'];
        
        // H.E.Rose, A Course on Finite Groups, page 266
        var a1 = gl39[
            x + 2, 1, 1,
            1, 2, 0,
            1, 0, 0
        ];
        var b1 = gl39[
            2 * x + 1, 2 * x + 1, 1,
            x, 2, 0,
            1, 0, 0
        ];

        var u3_3gl = Group.Generate("U3(3)gl", gl39, a1, b1);
        DisplayGroup.Head(u3_3gl); // |U3(3)| = 6048

        var s28 = new Sn(28);
        var a2 = s28[(1, 5, 7, 3, 12, 24, 11), (2, 23, 4, 27, 13, 14, 26), (6, 20, 18, 8, 25, 21, 28), (9, 10, 17, 15, 22, 16, 19)];
        var b2 = s28[(3, 4), (5, 17, 7, 16, 8, 20, 6, 13), (9, 19, 11, 14, 12, 18, 10, 15), (21, 23, 26, 28, 24, 22, 27, 25)];
        var u3_3pg = Group.Generate("U3(3)pg", s28, a2, b2);
        DisplayGroup.Head(u3_3pg);

        DisplayGroup.AreIsomorphics(u3_3gl, u3_3pg);
    }

    public static void U3_4()
    {
        var gl316 = new GLnq(3, 16);
        var x = gl316.Fq['x'];
        var x3 = x.Pow(3);
        var x11 = x.Pow(3) + x.Pow(2) + x;
        // gap> GeneratorsOfGroup(SU(3,4));
        var a = gl316[
            x, 0, 0,
            0, x3, 0,
            0, 0, x11
        ];
        var b = gl316[
            x, 1, 1,
            1, 1, 0,
            1, 0, 0
        ];

        var u3_4 = Group.Generate("U3(4)", gl316, a, b);
        DisplayGroup.Head(u3_4); // |U3(4)| = 62400
    }

    public static void U4_2()
    {
        var gl44 = new GLnq(4, 4);
        var x = gl44.Fq['x'];

        // gap> GeneratorsOfGroup(SU(4,2));
        var a = gl44[
            x, 0, 0, 0,
            0, x * x, 0, 0,
            0, 0, x * x, 0,
            0, 0, 0, x];
        var b = gl44[
            1, 0, 1, 0,
            1, 0, 0, 0,
            0, 1, 0, 1,
            0, 1, 0, 0];

        var u4_2 = Group.Generate("U4(2)", gl44, a, b);
        DisplayGroup.HeadOrders(u4_2);

        var gl43 = new GL(4, 3);
        var a1 = gl43[
            2, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 2
        ];
        var b1 = gl43[
            1, 0, 1, 0,
            1, 0, 0, 0,
            0, 1, 0, 1,
            0, 2, 0, 0
        ];

        var sp43 = Group.Generate("SP4(3)", gl43, a1, b1);
        DisplayGroup.Head(sp43);
        var z = Group.Zentrum(sp43);
        var psp43 = sp43.Over(z, "PSP4(3)");
        DisplayGroup.HeadOrders(psp43);

        var oa = u4_2.ElementsOrders[a]; // ord(a) = 3
        var ob = u4_2.ElementsOrders[b]; // ord(b) = 6
        var allea = psp43.Where(e => psp43.ElementsOrders[e] == oa).Order().ToArray();
        var eb = psp43.Where(e => psp43.ElementsOrders[e] == ob).Min();
        foreach (var ea in allea)
        {
            var pMap = Group.PartialMap((a, ea), (b, eb));
            var iso = Group.IsomorphismMap(u4_2, psp43, pMap);
            if (iso.Count != 0)
            {
                Console.WriteLine($"{u4_2} is isomorphic to {psp43}");
                break;
            }
        }
    }
}