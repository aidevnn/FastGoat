using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;

namespace FastGoat.Examples;

// H.E.Rose, A Course on Finite Groups, page 266 
public static class UnitaryGroup
{
    /// <summary>
    /// Performs operations related to the U3(2) the unitary group of 3x3 matrix in finite field F4 and
    /// then checks its isomorphism with (C3 x C3) x: Q8 group.
    /// </summary>
    public static void U3_2()
    {
        // Initialize GLnq(3, 4) and 'x' variable
        var gl34 = new GLnq(3, 4);
        var x = gl34.Fq['x'];

        // gap> GeneratorsOfGroup(SU(3,2));
        // Define matrix 'a' for SU3(2)
        var a = gl34[
            1, x, x,
            0, 1, x + 1,
            0, 0, 1
        ];

        // Define matrix 'b' for SU3(2)
        var b = gl34[
            x, 1, 1,
            1, 1, 0,
            1, 0, 0
        ];

        // Generate the SU3(2) special unitary group with 'a' and 'b'
        var su3_2 = Group.Generate("SU3(2)", gl34, a, b);
        DisplayGroup.Head(su3_2);

        // Calculate the Zentrum (central subgroup) of SU3(2)
        var z = Group.Zentrum(su3_2);
        DisplayGroup.Head(z);

        // Form the quotient group U3(2) over the Zentrum of SU3(2)
        var u3_2 = su3_2.Over(z, "U3(2)");
        DisplayGroup.Head(u3_2);

        // Create the semi-direct product group with Abelian(3,3) and Quaternion(8)
        var u3_2sdp = Group.SemiDirectProd(FG.Abelian(3, 3), FG.Quaternion(8));
        DisplayGroup.Head(u3_2sdp);

        // Check if U3(2) and the semi-direct product are isomorphic
        DisplayGroup.AreIsomorphics(u3_2, u3_2sdp);
    }

    /// <summary>
    /// Performs operations related to the U3(3) group, which represents the unitary group in its 3x3 matrix group over
    /// finite field F9 form, and also in its permutations group forms within Symmetric group 28, and checks their isomorphism.
    /// </summary>
    public static void U3_3()
    {
        // Initialize GLnq(3, 9) and 'x' variable
        var gl39 = new GLnq(3, 9);
        var x = gl39.Fq['x'];

        // H.E.Rose, A Course on Finite Groups, page 266
        // Define matrix 'a1' for U3(3)gl
        var a1 = gl39[
            x + 2, 1, 1,
            1, 2, 0,
            1, 0, 0
        ];

        // Define matrix 'b1' for U3(3)gl
        var b1 = gl39[
            2 * x + 1, 2 * x + 1, 1,
            x, 2, 0,
            1, 0, 0
        ];

        // Generate the U3(3)gl group with 'a1' and 'b1'
        var u3_3gl = Group.Generate("U3(3)gl", gl39, a1, b1);

        // Display the order of U3(3)gl
        DisplayGroup.Head(u3_3gl); // |U3(3)gl| = 6048

        // Initialize Symmetric group S28
        var s28 = new Sn(28);

        // Define permutations 'a2' and 'b2' for U3(3)pg in S28
        var a2 = s28[(1, 5, 7, 3, 12, 24, 11), (2, 23, 4, 27, 13, 14, 26), (6, 20, 18, 8, 25, 21, 28), (9, 10, 17, 15, 22, 16, 19)];
        var b2 = s28[(3, 4), (5, 17, 7, 16, 8, 20, 6, 13), (9, 19, 11, 14, 12, 18, 10, 15), (21, 23, 26, 28, 24, 22, 27, 25)];

        // Generate the U3(3)pg group with 'a2' and 'b2' in S28
        var u3_3pg = Group.Generate("U3(3)pg", s28, a2, b2);

        // Display the order of U3(3)pg
        DisplayGroup.Head(u3_3pg);

        // Check if U3(3)gl and U3(3)pg are isomorphic
        DisplayGroup.AreIsomorphics(u3_3gl, u3_3pg);
    }

    /// <summary>
    /// Performs operations related to the U3(4) group, which is the unitary group as a 3x3 matrix group over F16.
    /// </summary>
    public static void U3_4()
    {
        // Initialize GLnq(3, 16) and 'x' variable
        var gl316 = new GLnq(3, 16);
        var x = gl316.Fq['x'];
        var x3 = x.Pow(3);
        var x11 = x.Pow(3) + x.Pow(2) + x;

        // gap> GeneratorsOfGroup(SU(3,4));
        // Define matrix 'a' for U3(4)
        var a = gl316[
            x, 0, 0,
            0, x3, 0,
            0, 0, x11
        ];

        // Define matrix 'b' for U3(4)
        var b = gl316[
            x, 1, 1,
            1, 1, 0,
            1, 0, 0
        ];

        // Generate the U3(4) group with 'a' and 'b'
        var u3_4 = Group.Generate("U3(4)", gl316, a, b);

        // Display the order of U3(4)
        DisplayGroup.Head(u3_4); // |U3(4)| = 62400
    }

    /// <summary>
    /// Performs operations related to U4(2) the unitary group of 4x4 matrix over finite field F16 and PSP4(3) the projective
    /// sympletic group of 4x4 matrix over F3 and attempts to find an isomorphism between them.
    /// </summary>
    public static void U4_2()
    {
        // Initialize GLnq(4, 4) and 'x' variable
        var gl44 = new GLnq(4, 4);
        var x = gl44.Fq['x'];

        // gap> GeneratorsOfGroup(SU(4,2));
        // Define matrix 'a' for U4(2)
        var a = gl44[
            x, 0, 0, 0,
            0, x * x, 0, 0,
            0, 0, x * x, 0,
            0, 0, 0, x];

        // Define matrix 'b' for U4(2)
        var b = gl44[
            1, 0, 1, 0,
            1, 0, 0, 0,
            0, 1, 0, 1,
            0, 1, 0, 0];

        // Generate the U4(2) group with 'a' and 'b'
        var u4_2 = Group.Generate("U4(2)", gl44, a, b);
        DisplayGroup.HeadOrders(u4_2);

        // Initialize GL(4, 3) and matrices 'a1' and 'b1' for SP4(3)
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

        // Generate the SP4(3) group with 'a1' and 'b1'
        var sp43 = Group.Generate("SP4(3)", gl43, a1, b1);
        DisplayGroup.Head(sp43);

        // Calculate the Zentrum (central subgroup) of SP4(3)
        var z = Group.Zentrum(sp43);

        // Form the quotient group PSP4(3) over the Zentrum of SP4(3)
        var psp43 = sp43.Over(z, "PSP4(3)");
        DisplayGroup.HeadOrders(psp43);

        // Calculate the orders of elements 'a' and 'b' in U4(2)
        var oa = u4_2.ElementsOrders[a]; // ord(a) = 3
        var ob = u4_2.ElementsOrders[b]; // ord(b) = 6

        // Find elements in PSP4(3) with orders equal to 'oa' and 'ob'
        var allea = psp43.Where(e => psp43.ElementsOrders[e] == oa).Order().ToArray();
        var eb = psp43.Where(e => psp43.ElementsOrders[e] == ob).Min();

        // Check for isomorphism between U4(2) and PSP4(3)
        foreach (var ea in allea)
        {
            // Create a partial map and attempt to find an isomorphism
            var pMap = Group.PartialMap((a, ea), (b, eb));
            var iso = Group.IsomorphismMap(u4_2, psp43, pMap);

            // If an isomorphism is found, print the result and exit the loop
            if (iso.Count != 0)
            {
                Console.WriteLine($"{u4_2} is isomorphic to {psp43}");
                break;
            }
        }
    }
}