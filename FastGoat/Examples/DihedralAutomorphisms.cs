using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;

namespace FastGoat.Examples;

public static class DihedralAutomorphisms
{
    static ConcreteGroup<Perm> DnPerm(int n = 16)
    {
        int m = (n % 2) == 0 ? 1 : 2;
        var sn = new Sn(n);
        var an = Enumerable.Range(1, n).ToArray();
        var a2 = Enumerable.Range(m, n / 2).Select(i => (Tuple2Array)(i, n + m - i)).ToArray();
        var cn = sn.Cycle(an);
        var c2 = sn.ComposesCycles(a2);
        var d2n = Group.Generate("D2n", sn, c2, cn);
        return d2n;
    }

    static SemiDirectProduct<ZnInt, ZnInt> DnSdp(int n = 16)
    {
        var cn = new Cn(n);
        var c2 = new Cn(2);
        var un = new Un(n);
        var theta = new Dictionary<ZnInt, Automorphism<ZnInt>>()
        {
            [c2[0]] = un.Neutral(),
            [c2[1]] = un[(cn[1], cn[n - 1])]
        };

        var homTheta = Group.Hom(c2, theta);
        var d2n = Group.SemiDirectProd("D2n", cn, homTheta, c2);
        return d2n;
    }

    public static void Dn()
    {
        int n = 3;
        var d2n = DnSdp(n);
        var d2n2 = DnPerm(n);
        DisplayGroup.HeadSdp(d2n);
        DisplayGroup.Head(d2n2);
        Console.WriteLine("IsIsomorphic : {0}", d2n2.IsIsomorphicTo(d2n));
    }

    public static void AutDn()
    {
        for (int n = 6; n <= 32; n += 2)
        {
            var d2n = DnSdp(n / 2);
            var autD2N = Group.AllAutomorphisms(d2n);

            var a = d2n[1, 0];
            var b = d2n[0, 1];
            var autH = autD2N.Where(aut => aut[a].Equals(a)).ToArray();
            var autK = autD2N.Where(aut => aut[b].Equals(b)).ToArray();

            var phi = new Un(n / 2).Count();
            Console.WriteLine("|Aut(D{0})| = {1}; phi({2}) = {3}", n, autD2N.Count, n / 2, phi);
            Console.WriteLine($"a[{d2n.ElementsOrders[a]}]={a} b[{d2n.ElementsOrders[b]}]={b}");
            Console.WriteLine("|y in Aut(D{0}), y(a)=a| = {1,-3}; |y in Aut(D{0}), y(b) = b| = {2,-3}", n, autH.Count(),
                autK.Count());
            Console.WriteLine();
        }

        // gap> for i in [2..16] do
        // > dn:=DihedralGroup(2*i);
        // > Print("Nb Aut(D", 2*i, ") = ", Size(AllAutomorphisms(dn)),"\n");
        // > od;
        // Nb Aut(D4) = 6
        // Nb Aut(D6) = 6
        // Nb Aut(D8) = 8
        // Nb Aut(D10) = 20
        // Nb Aut(D12) = 12
        // Nb Aut(D14) = 42
        // Nb Aut(D16) = 32
        // Nb Aut(D18) = 54
        // Nb Aut(D20) = 40
        // Nb Aut(D22) = 110
        // Nb Aut(D24) = 48
        // Nb Aut(D26) = 156
        // Nb Aut(D28) = 84
        // Nb Aut(D30) = 120
        // Nb Aut(D32) = 128
    }

    public static void AutDnPerm()
    {
        for (int n = 6; n <= 32; n += 2)
        {
            var d2n = DnPerm(n / 2);
            var autD2N = Group.AllAutomorphisms(d2n);

            var gens = d2n.GetGenerators().ToArray();
            var a = gens.First(e => d2n.ElementsOrders[e] == n / 2);
            var b = gens.First(e => d2n.ElementsOrders[e] == 2);
            var autH = autD2N.Where(aut => aut[a].Equals(a)).ToArray();
            var autK = autD2N.Where(aut => aut[b].Equals(b)).ToArray();

            var phi = new Un(n / 2).Count();
            Console.WriteLine("|Aut(D{0})| = {1}; phi({2}) = {3}", n, autD2N.Count, n / 2, phi);
            Console.WriteLine($"a[{d2n.ElementsOrders[a]}]={a} b[{d2n.ElementsOrders[b]}]={b}");
            Console.WriteLine("|y in Aut(D{0}), y(a)=a| = {1,-3}; |y in Aut(D{0}), y(b) = b| = {2,-3}", n, autH.Count(),
                autK.Count());
            Console.WriteLine();
        }
    }
}