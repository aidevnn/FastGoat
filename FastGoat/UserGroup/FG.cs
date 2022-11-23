using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words;

namespace FastGoat.UserGroup;

public static partial class FG
{
    public static ConcreteGroup<Perm> Symmetric(int n) => new Symm(n);

    public static ConcreteGroup<Perm> Alternate(int n)
    {
        if (n < 3)
            throw new GroupException(GroupExceptionType.GroupDef);
        
        var sn = new Symm(n);
        var gi = (n - 2).Range(3).Select(i => sn[(1, 2, i)]).ToArray();
        return Group.Generate($"Alt{n}", sn, gi);
    }

    public static ConcreteGroup<Perm> Dihedral(int n)
    {
        var m = (n % 2) == 0 ? 1 : 2;
        var sn = new Sn(n);
        var an = Enumerable.Range(1, n).ToArray();
        var a2 = Enumerable.Range(m, n / 2).Select(i => (Tuple2Array)(i, n + m - i)).ToArray();
        var cn = sn.Cycle(an);
        var c2 = sn.ComposesCycles(a2);
        var d2n = Group.Generate($"D{2 * n}", sn, c2, cn);
        return d2n;
    }

    public static ConcreteGroup<Perm> PermGroup(string name, int n, params ValueType[] generators)
    {
        var sn = new Sn(n);
        var gi = generators.Select(g => sn.ComposesCycles(Tuple2Array.ComplexTuples(g))).ToArray();
        return Group.Generate(name, sn, gi);
    }

    public static ConcreteGroup<Perm> PermGroup(int n, params ValueType[] generators) => PermGroup("G", n, generators);

    public static ConcreteGroup<Ep<ZnInt>> Abelian(string name, params int[] seq)
    {
        return Group.Create(name, new Gp<ZnInt>(seq.Select(i => new Zn(i)).Cast<IGroup<ZnInt>>().ToArray()));
    }

    public static ConcreteGroup<Ep<ZnInt>> Abelian(params int[] seq)
    {
        return Product.GpGenerate(seq.Select(i => new Cn(i)).Cast<IGroup<ZnInt>>().ToArray());
    }

    public static SemiDirectProduct<ZnInt, ZnInt> DihedralSdp(int n)
    {
        var cn = new Cn(n);
        var autCn = Group.AutBase(cn);
        var y = autCn[(cn[1], cn[n - 1])];
        var aut = Group.Generate(autCn, y);
        
        var c2 = new Cn(2);
        var pMap = Group.PartialMap((c2[1], y));
        var theta =Group. Hom(c2, Group.HomomorphismMap(c2, aut, pMap));
        return Group.SemiDirectProd($"D{2 * n}", cn, theta, c2);
    }

    public static WordGroup DihedralWg(int n) => new($"D{2 * n}", $"a{n}, b2, abab");
}