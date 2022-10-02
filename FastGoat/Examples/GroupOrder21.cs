using FastGoat.Gp;
using FastGoat.UserGroup;

namespace FastGoat.Examples;

public static class GroupOrder21
{
    public static void Zn21()
    {
        var z21 = new Zn(21);
        var g21 = Group.Generate(z21[1]);
        DisplayGroup.HeadElements(g21);
    }

    public static void ZnDirectProduct()
    {
        var z7Xz3 = Product.Group(new Zn(7), new Zn(3));
        var g21 = Group.Generate(z7Xz3[1, 0], z7Xz3[0, 1]);
        DisplayGroup.HeadElements(g21);
    }

    public static void Symmetric7()
    {
        var s7 = new Sn(7);
        var a = s7[(1, 2, 3, 4, 5, 6, 7)];
        var allC3 = IntExt.GetPermutations(7).Select(s7.CreateElement).Where(p => (p ^ 3) == s7.Neutral());
        var b = allC3.First(p => Group.GenerateElements(s7, a, p).Count() == 21);
        var g21 = Group.Generate(a, b);
        DisplayGroup.HeadElements(g21);
    }

    public static void Symmetric7Fast()
    {
        var s7 = new Sn(7);
        var a = s7[(1, 2, 3, 4, 5, 6, 7)];
        var allC3 = IntExt.GetPermutations(7).Select(s7.CreateElement).Where(p => (p ^ 3) == s7.Neutral());
        var b = allC3.First(p => (a ^ 2) == p * a * (p ^ -1));
        var g21 = Group.Generate(a, b);
        DisplayGroup.HeadElements(g21);
    }
    public static void SemiDirectProduct()
    {
        var c7 = Group.Generate("C7", new Zn(7)[1]);
        var c3 = Group.Generate("C3", new Zn(3)[1]);
        var g21 = Group.SemiDirectProd(c7, c3);
        DisplayGroup.HeadElementsSdp(g21);
    }
}