using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class GLnK
{
    public static void ExampleDihedralInReal()
    {
        var gl = FG.GLnK($"R", 2, Rational.KZero());

        var m0 = gl[0, 1, -1, 0];
        var m1 = gl[0, 1, 1, 0];

        var g0 = Group.Generate("R(π/2)", gl, m0, m1);
        DisplayGroup.HeadElements(g0);
        var g1 = Group.Generate("R(y=x)", gl, m1);
        DisplayGroup.HeadElements(g1);
        
        var g2 = Group.Generate("D8a", gl, m0, m1);
        DisplayGroup.HeadElements(g2);

        DisplayGroup.AreIsomorphics(g2, FG.Dihedral(4));
    }

    public static void ExampleRotation1()
    {
        var x = FG.QPoly();
        var a = FG.NumberFieldQ(x.Pow(2) - 2, "√2");
        var gl = FG.GLnK($"Q({a})", 2, a);

        var m0 = gl[a / 2, a / 2, -a / 2, a / 2];
    
        var c8 = Group.Generate("R(π/4)", gl, m0);
        DisplayGroup.HeadElements(c8);
    
        DisplayGroup.AreIsomorphics(c8, new Cn(8));
    }

    public static void ExampleRotation2()
    {
        var x = FG.QPoly();
        var (a, b) = FG.NumberFieldQ((x.Pow(2) - 2, "√2"), (x.Pow(2) - 3, "√3"));
        var gl = FG.GLnK($"Q({a},{b})", 2, a);

        var r1 = b.One / 2;
        var m0 = gl[a / 2, a / 2, -a / 2, a / 2];
        var m1 = gl[r1, b / 2, -b / 2, r1];

        var g0 = Group.Generate("R(π/3)", gl, m1);
        DisplayGroup.HeadElements(g0);
        var g1 = Group.Generate("R(π/12)", gl, m0, m1);
        DisplayGroup.HeadElements(g1);
    
        DisplayGroup.AreIsomorphics(g1, new Cn(24));
    }

    public static void ExampleDihedralComplex()
    {
        var x = FG.QPoly();
        var i = FG.NumberFieldQ(x.Pow(2) + 1, "i");
        var gl = FG.GLnK($"Q({i})", 2, i);
    
        var m0 = gl[i, 0, 0, -i];
        var m1 = gl[0, 1, 1, 0];
    
        var g0 = Group.Generate("C4", gl, m0);
        DisplayGroup.HeadElements(g0);
    
        var g1 = Group.Generate("D8a", gl, m0, m1);
        DisplayGroup.HeadElements(g1);
    
        DisplayGroup.AreIsomorphics(g1, FG.Dihedral(4));
    }
    
    public static void ExampleQuaternionComplex()
    {
        var x = FG.QPoly();
        var i = FG.NumberFieldQ(x.Pow(2) + 1, "i");
        var gl = FG.GLnK($"Q({i})", 2, i);
    
        var m0 = gl[i, 0, 0, -i];
        var m1 = gl[0, 1, -1, 0];
    
        var g1 = Group.Generate("Q8a", gl, m0, m1);
        DisplayGroup.HeadElements(g1);
    
        DisplayGroup.AreIsomorphics(g1, FG.Quaternion(8));
    }

    public static void ExampleDicyclic4()
    {
        var x = FG.QPoly();
        var (a, i) = FG.NumberFieldQ((x.Pow(2) - 2, "√2"), (x.Pow(2) + 1, "i"));
        var gl = FG.GLnK($"Q({i},{a})", 2, a);

        var m0 = gl[a / 2, a / 2, -a / 2, a / 2];
        var m1 = gl[i, 0, 0, -i];
        var m2 = gl[0, 1, -1, 0];

        var g0 = Group.Generate("C8", gl, m0);
        DisplayGroup.HeadElements(g0);
        var g1 = Group.Generate("Q8", gl, m1, m2);
        DisplayGroup.HeadElements(g1);
        var g2 = Group.Generate("Dic4a", gl, m0, m1, m2);
        DisplayGroup.HeadElements(g2);
    
        DisplayGroup.AreIsomorphics(g2, FG.DiCyclic(4));
    }

    public static void ExampleSemiDirectProd()
    {
        var x = FG.QPoly();
        var (a, i) = FG.NumberFieldQ((x.Pow(2) - 3, "√3"), (x.Pow(2) + 1, "i"));
        var gl = FG.GLnK($"Q({i},{a})", 2, a);
        
        var m0 = gl[(1 + i * a) / 2, 0, 0, (1 - i * a) / 2];
        var m1 = gl[i, 0, 0, -i];
        var m2 = gl[0, 1, -1, 0];
        
        var g0 = Group.Generate("C6", gl, m0);
        DisplayGroup.HeadElements(g0);
        var g1 = Group.Generate("Q8", gl, m1, m2);
        DisplayGroup.HeadElements(g1);
        var g2 = Group.Generate("G24", gl, m0, m1, m2);
        DisplayGroup.HeadElements(g2);
        
        var sdp = Group.SemiDirectProd(new Cn(3), FG.Quaternion(8));
        DisplayGroup.HeadOrders(g2);
        DisplayGroup.HeadSdpOrders(sdp);
        DisplayGroup.AreIsomorphics(g2, sdp);
    }
}
