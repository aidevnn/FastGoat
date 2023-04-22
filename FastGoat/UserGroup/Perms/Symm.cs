using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.UserGroup.Perms;

public class Symm : ConcreteGroup<Perm>
{
    public Symm(int n) : base($"Symm{n}",
        n is <= 8 and > 1 ? new Sn(n) : throw new GroupException(GroupExceptionType.GroupDef))
    {
        N = n;
    }

    public int N { get; }

    public IEnumerable<(string name, Perm[] gens)> TransitiveSubGroups()
    {
        if (N == 2)
        {
            yield return ("S2", new[] { this[(1, 2)] }); // order = 2
        }

        if (N == 3)
        {
            yield return ("A3", new[] { this[(1, 2, 3)] }); // order = 3
            yield return ("S3", new[] { this[(1, 2, 3)], this[(1, 2)] }); // order = 6
        }

        if (N == 4)
        {
            yield return ("C(4) = 4", new[] { this[(1, 2, 3, 4)] }); // order = 4
            yield return ("E(4) = 2[x]2", new[] { this[(1, 4), (2, 3)], this[(1, 2), (3, 4)] }); // order = 4
            yield return ("D(4)", new[] { this[(1, 2, 3, 4)], this[(1, 3)] }); // order = 8
            yield return ("A4", new[] { this[(1, 2, 3)], this[(2, 3, 4)] }); // order = 12
            yield return ("S4", new[] { this[(1, 2, 3, 4)], this[(1, 2)] }); // order = 24
        }

        if (N == 5)
        {
            yield return ("C(5) = 5", new[] { this[(1, 2, 3, 4, 5)] }); // order = 5
            yield return ("D(5) = 5:2", new[] { this[(1, 2, 3, 4, 5)], this[(1, 4), (2, 3)] }); // order = 10
            yield return ("F(5) = 5:4", new[] { this[(1, 2, 3, 4, 5)], this[(1, 2, 4, 3)] }); // order = 20
            yield return ("A5", new[] { this[(1, 2, 3, 4, 5)], this[(3, 4, 5)] }); // order = 60
            yield return ("S5", new[] { this[(1, 2, 3, 4, 5)], this[(1, 2)] }); // order = 120
        }

        if (N == 6)
        {
            yield return ("C(6) = 6 = 3[x]2", new[] { this[(1, 2, 3, 4, 5, 6)] }); // order = 6
            yield return ("D_6(6) = [3]2", new[] { this[(1, 3, 5), (2, 4, 6)], this[(1, 4), (2, 3), (5, 6)] }); // order = 6
            yield return ("D(6) = S(3)[x]2", new[] { this[(1, 2, 3, 4, 5, 6)], this[(1, 4), (2, 3), (5, 6)] }); // order = 12
            yield return ("A_4(6) = [2^2]3", new[] { this[(1, 4), (2, 5)], this[(1, 3, 5), (2, 4, 6)] }); // order = 12
            yield return ("F_18(6) = [3^2]2 = 3 wr 2", new[] { this[(2, 4, 6)], this[(1, 4), (2, 5), (3, 6)] }); // order = 18
            yield return ("2A_4(6) = [2^3]3 = 2 wr 3", new[] { this[(3, 6)], this[(1, 3, 5), (2, 4, 6)] }); // order = 24
            yield return ("S_4(6d) = [2^2]S(3)",
                new[] { this[(1, 4), (2, 5)], this[(1, 3, 5), (2, 4, 6)], this[(1, 5), (2, 4)] }); // order = 24
            yield return ("S_4(6c) = 1/2[2^3]S(3)",
                new[] { this[(1, 4), (2, 5)], this[(1, 3, 5), (2, 4, 6)], this[(1, 5), (2, 4), (3, 6)] }); // order = 24
            yield return ("F_18(6):2 = [1/2.S(3)^2]2",
                new[] { this[(2, 4, 6)], this[(1, 5), (2, 4)], this[(1, 4), (2, 5), (3, 6)] }); // order = 36
            yield return ("F_36(6) = 1/2[S(3)^2]2",
                new[] { this[(2, 4, 6)], this[(1, 5), (2, 4)], this[(1, 4, 5, 2), (3, 6)] }); // order = 36
            yield return ("2S_4(6) = [2^3]S(3) = 2 wr S(3)",
                new[] { this[(3, 6)], this[(1, 3, 5), (2, 4, 6)], this[(1, 5), (2, 4)] }); // order = 48
            yield return ("L(6) = PSL(2,5) = A_5(6)", new[] { this[(1, 2, 3, 4, 6)], this[(1, 4), (5, 6)] }); // order = 60
            yield return ("F_36(6):2 = [S(3)^2]2 = S(3) wr 2",
                new[] { this[(2, 4, 6)], this[(2, 4)], this[(1, 4), (2, 5), (3, 6)] }); // order = 72
            yield return ("L(6):2 = PGL(2,5) = S_5(6)", new[] { this[(1, 2, 3, 4, 6)], this[(1, 2), (3, 4), (5, 6)] }); // order = 120
            yield return ("A6", new[] { this[(1, 2, 3, 4, 5)], this[(4, 5, 6)] }); // order = 360
            yield return ("S6", new[] { this[(1, 2, 3, 4, 5, 6)], this[(1, 2)] }); // order = 720
        }

        if (N == 7)
        {
            yield return ("C(7) = 7", new[] { this[(1, 2, 3, 4, 5, 6, 7)] }); // order = 7
            yield return ("D(7) = 7:2", new[] { this[(1, 2, 3, 4, 5, 6, 7)], this[(1, 6), (2, 5), (3, 4)] }); // order = 14
            yield return ("F_21(7) = 7:3", new[] { this[(1, 2, 3, 4, 5, 6, 7)], this[(1, 2, 4), (3, 6, 5)] }); // order = 21
            yield return ("F_42(7) = 7:6", new[] { this[(1, 2, 3, 4, 5, 6, 7)], this[(1, 3, 2, 6, 4, 5)] }); // order = 42
            yield return ("L(7) = L(3,2)", new[] { this[(1, 2, 3, 4, 5, 6, 7)], this[(1, 2), (3, 6)] }); // order = 168
            yield return ("A7", new[] { this[(1, 2, 3, 4, 5, 6, 7)], this[(5, 6, 7)] }); // order = 2520
            yield return ("S7", new[] { this[(1, 2, 3, 4, 5, 6, 7)], this[(1, 2)] }); // order = 5040
        }

        if (N == 8)
        {
            yield return ("C(8)=8", new[] { this[(1, 2, 3, 4, 5, 6, 7, 8)] }); // order = 8
            yield return ("4[x]2", new[] { this[(1, 2, 3, 8), (4, 5, 6, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)] }); // order = 8
            yield return ("E(8)=2[x]2[x]2",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)]
                }); // order = 8
            yield return ("D_8(8)=[4]2",
                new[] { this[(1, 2, 3, 8), (4, 5, 6, 7)], this[(1, 6), (2, 5), (3, 4), (7, 8)] }); // order = 8
            yield return ("Q_8(8)", new[] { this[(1, 2, 3, 8), (4, 5, 6, 7)], this[(1, 7, 3, 5), (2, 6, 8, 4)] }); // order = 8
            yield return ("D(8)", new[] { this[(1, 2, 3, 4, 5, 6, 7, 8)], this[(1, 6), (2, 5), (3, 4), (7, 8)] }); // order = 16
            yield return ("1/2[2^3]4", new[] { this[(1, 2, 3, 4, 5, 6, 7, 8)], this[(1, 5), (3, 7)] }); // order = 16
            yield return ("2D_8(8)=[D(4)]2", new[] { this[(1, 2, 3, 4, 5, 6, 7, 8)], this[(1, 3), (2, 6), (5, 7)] }); // order = 16
            yield return ("E(8):2=D(4)[x]2",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(4, 5), (6, 7)]
                }); // order = 16
            yield return ("[2^2]4", new[] { this[(1, 5), (3, 7)], this[(1, 2, 3, 8), (4, 5, 6, 7)] }); // order = 16
            yield return ("1/2[2^3]E(4)=Q_8:2",
                new[] { this[(1, 5), (3, 7)], this[(1, 3, 5, 7), (2, 4, 6, 8)], this[(1, 4, 5, 8), (2, 3, 6, 7)] }); // order = 16
            yield return ("2A_4(8)=[2]A(4)=SL(2,3)",
                new[] { this[(1, 3, 5, 7), (2, 4, 6, 8)], this[(1, 3, 8), (4, 5, 7)] }); // order = 24
            yield return ("E(8):3=A(4)[x]2",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 2, 3), (4, 6, 5)]
                }); // order = 24
            yield return ("S(4)[1/2]2=1/2(S_4[x]2)",
                new[]
                {
                    this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 2, 3), (5, 6, 7)], this[(1, 4), (2, 6), (3, 7), (5, 8)]
                }); // order = 24
            yield return ("[1/4.cD(4)^2]2",
                new[] { this[(1, 2, 3, 4, 5, 6, 7, 8)], this[(1, 5), (3, 7)], this[(1, 6), (2, 5), (3, 4), (7, 8)] }); // order = 32
            yield return ("1/2[2^4]4", new[] { this[(2, 6), (3, 7)], this[(1, 2, 3, 4, 5, 6, 7, 8)] }); // order = 32
            yield return ("[4^2]2", new[] { this[(1, 2, 3, 8)], this[(1, 5), (2, 6), (3, 7), (4, 8)] }); // order = 32
            yield return ("E(8):E_4=[2^2]D(4)",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(4, 5), (6, 7)], this[(4, 6), (5, 7)]
                }); // order = 32
            yield return ("E(8):4=[1/4.eD(4)^2]2",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 3), (4, 5, 6, 7)]
                }); // order = 32
            yield return ("[2^3]4", new[] { this[(2, 6), (3, 7)], this[(1, 2, 3, 8), (4, 5, 6, 7)] }); // order = 32
            yield return ("1/2[2^4]E(4)=[1/4.dD(4)^2]2",
                new[]
                {
                    this[(1, 5), (3, 7)], this[(1, 4, 5, 8), (2, 3), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)]
                }); // order = 32
            yield return ("E(8):D_4=[2^3]2^2",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(2, 3), (4, 5)], this[(2, 3), (6, 7)]
                }); // order = 32
            yield return ("2S_4(8)=GL(2,3)", new[] { this[(1, 2, 3, 4, 5, 6, 7, 8)], this[(1, 3, 8), (4, 5, 7)] }); // order = 48
            yield return ("E(8):D_6=S(4)[x]2",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 2, 3), (4, 6, 5)], this[(2, 3), (4, 5)]
                }); // order = 48
            yield return ("E(8):7=F_56(8)",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 2, 6, 3, 4, 5, 7)]
                }); // order = 56
            yield return ("1/2[2^4]eD(4)",
                new[] { this[(1, 2, 3, 4, 5, 6, 7, 8)], this[(1, 5), (4, 8)], this[(1, 7), (3, 5), (4, 8)] }); // order = 64
            yield return ("[2^4]4", new[] { this[(4, 8)], this[(1, 2, 3, 8), (4, 5, 6, 7)] }); // order = 64
            yield return ("1/2[2^4]dD(4)",
                new[] { this[(2, 6), (3, 7)], this[(1, 3), (5, 7)], this[(1, 2, 3, 4, 5, 6, 7, 8)] }); // order = 64
            yield return ("E(8):D_8=[2^3]D(4)",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 3), (4, 5, 6, 7)], this[(1, 3), (5, 7)]
                }); // order = 64
            yield return ("1/2[2^4]cD(4)",
                new[] { this[(2, 6), (3, 7)], this[(1, 3), (4, 8), (5, 7)], this[(1, 2, 3, 8), (4, 5, 6, 7)] }); // order = 64
            yield return ("[2^4]E(4)",
                new[] { this[(4, 8)], this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)] }); // order = 64
            yield return ("[2^3]A(4)",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 2, 3), (4, 6, 5)], this[(2, 5), (3, 4)]
                }); // order = 96
            yield return ("E(8):A_4=[1/3.A(4)^2]2=E(4):6",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 2, 3), (4, 6, 5)], this[(4, 6), (5, 7)]
                }); // order = 96
            yield return ("1/2[E(4)^2:S_3]2=E(4)^2:D_6",
                new[] { this[(1, 8), (2, 3)], this[(1, 2, 3), (5, 6, 7)], this[(1, 5), (2, 7), (3, 6), (4, 8)] }); // order = 96
            yield return ("[2^4]D(4)", new[] { this[(4, 8)], this[(1, 3), (5, 7)], this[(1, 2, 3, 8), (4, 5, 6, 7)] }); // order = 128
            yield return ("E(8):F_21",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 2, 6, 3, 4, 5, 7)], this[(1, 2, 3), (4, 6, 5)]
                }); // order = 168
            yield return ("L(8)=PSL(2,7)",
                new[]
                {
                    this[(1, 2, 3, 4, 5, 6, 8)], this[(1, 2, 4), (3, 6, 5)], this[(1, 6), (2, 3), (4, 5), (7, 8)]
                }); // order = 168
            yield return ("[2^4]A(4)",
                new[] { this[(4, 8)], this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 2, 3), (5, 6, 7)] }); // order = 192
            yield return ("[2^3]S(4)",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 2, 3), (4, 6, 5)], this[(1, 6), (2, 3, 5, 4)]
                }); // order = 192
            yield return ("1/2[2^4]S(4)",
                new[]
                {
                    this[(1, 5), (4, 8)], this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 2, 3), (5, 6, 7)],
                    this[(2, 3), (4, 8), (6, 7)]
                }); // order = 192
            yield return ("E(8):S_4=[E(4)^2:S_3]2=E(4)^2:D_12",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 2, 3), (4, 6, 5)], this[(1, 3), (4, 5, 6, 7)]
                }); // order = 192
            yield return ("[A(4)^2]2",
                new[] { this[(1, 3), (2, 8)], this[(1, 2, 3)], this[(1, 5), (2, 6), (3, 7), (4, 8)] }); // order = 288
            yield return ("L(8):2=PGL(2,7)",
                new[] { this[(1, 2, 3, 4, 5, 6, 8)], this[(1, 3, 2, 6, 4, 5)], this[(1, 6), (2, 3), (4, 5), (7, 8)] }); // order = 336
            yield return ("[2^4]S(4)", new[] { this[(4, 8)], this[(1, 8), (4, 5)], this[(1, 2, 3, 8), (4, 5, 6, 7)] }); // order = 384
            yield return ("[1/2.S(4)^2]2",
                new[]
                {
                    this[(1, 3), (2, 8)], this[(1, 2, 3)], this[(1, 8), (4, 5)], this[(1, 5), (2, 6), (3, 7), (4, 8)]
                }); // order = 576
            yield return ("1/2[S(4)^2]2",
                new[]
                {
                    this[(1, 3), (2, 8)], this[(1, 2, 3)], this[(1, 8), (4, 5)], this[(1, 5), (2, 7, 3, 6), (4, 8)]
                }); // order = 576
            yield return ("[S(4)^2]2",
                new[] { this[(1, 2, 3, 8)], this[(2, 3)], this[(1, 5), (2, 6), (3, 7), (4, 8)] }); // order = 1152
            yield return ("E(8):L_7=AL(8)",
                new[]
                {
                    this[(1, 8), (2, 3), (4, 5), (6, 7)], this[(1, 3), (2, 8), (4, 6), (5, 7)], this[(1, 5), (2, 6), (3, 7), (4, 8)],
                    this[(1, 2, 6, 3, 4, 5, 7)], this[(1, 2, 3), (4, 6, 5)], this[(1, 2), (5, 6)]
                }); // order = 1344
            yield return ("A8", new[] { this[(1, 2, 3, 4, 5, 6, 7)], this[(6, 7, 8)] }); // order = 20160
            yield return ("S8", new[] { this[(1, 2, 3, 4, 5, 6, 7, 8)], this[(1, 2)] }); // order = 40320
        }
    }
}