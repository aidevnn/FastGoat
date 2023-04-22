using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly('X');
    var P = x.Pow(5) - 2;
    var gal = GaloisApplicationsPart2.GaloisGroupChebotarev(x.Pow(5) - 2, detail: true);

    var (X0, y0) = FG.EPolyXc(P, 'y');
    var P0 = P.Substitute(X0) / (X0 - y0);
    var (r, a, b) = IntFactorisation.PrimitiveElt(P0);
    var minPol = r.SubstituteChar('X');
    
    var (X1, y1) = FG.EPolyXc(minPol, 'y');
    Console.WriteLine(minPol);
    var def = IntFactorisation.Deflate(minPol, 5);
    Console.WriteLine(def);
    var facts = IntFactorisation.AlgebraicFactors(def.nf.Substitute(X1), true);
    var roots = facts.SelectMany(f => IntFactorisation.AlgebraicRoots(f.Substitute(X1.Pow(5)), true)).ToList();

    var subFields = GaloisTheory.SubFields(roots, 3).ToArray();
    var extTowers = GaloisApplications.ExtensionsTower(subFields);
    GaloisApplications.GaloisCorrespondence(extTowers);
}
/***
With P = X^20 + 156*X^15 + 39376*X^10 + 421776*X^5 + 2576816
Factorization in Q(y)/Q
    X + -y
    X + 1/272800*y^16 + 391/682000*y^11 + 47339/341000*y^6 + 226429/170500*y
    X + 3/426250*y^16 + 3673/3410000*y^11 + 231933/852500*y^6 + 2088549/852500*y
    X + -37/6820000*y^16 + -1301/1705000*y^11 + -344989/1705000*y^6 + 71221/213125*y
    X + 9/13640000*y^16 + 369/6820000*y^11 + 71843/3410000*y^6 + -3002093/1705000*y
    X + -3/682000*y^16 + -27/42625*y^11 + -7152/42625*y^6 + -76201/85250*y
    X + -111/13640000*y^16 + -8271/6820000*y^11 + -1072477/3410000*y^6 + -2640133/1705000*y
    X + 1/426250*y^16 + 29/77500*y^11 + 84131/852500*y^6 + 52309/38750*y
    X + 59/13640000*y^16 + 389/620000*y^11 + 545233/3410000*y^6 + 87927/155000*y
    X + -7/2728000*y^16 + -43/124000*y^11 + -67093/682000*y^6 + 67/248*y
    X + 57/13640000*y^16 + 4507/6820000*y^11 + 573219/3410000*y^6 + 4286241/1705000*y
    X + -19/13640000*y^16 + -1399/6820000*y^11 + -163793/3410000*y^6 + 1185483/1705000*y
    X + -17/13640000*y^16 + -1627/6820000*y^11 + -191779/3410000*y^6 + -2133561/1705000*y
    X + 17/2728000*y^16 + 251/272800*y^11 + 161771/682000*y^6 + 19733/341000*y
    X + 23/6820000*y^16 + 859/1705000*y^11 + 227171/1705000*y^6 + 25976/213125*y
    X + -7/3410000*y^16 + -221/852500*y^11 + -58909/852500*y^6 + 310322/213125*y
    X + -21/6820000*y^16 + -1791/3410000*y^11 + -214237/1705000*y^6 + -2059103/852500*y
    X + -7/1705000*y^16 + -203/310000*y^11 + -136573/852500*y^6 + -92019/77500*y
    X + -51/6820000*y^16 + -3951/3410000*y^11 + -500317/1705000*y^6 + -1968613/852500*y
    X + 1/124000*y^16 + 823/682000*y^11 + 1901/6200*y^6 + 208331/170500*y

Galois Group
|G1| = 20
Type        NonAbelianGroup
BaseGroup   S20
SuperGroup  |Gal( Q(y)/Q )| = 20

Elements
( 1)[1] = []
( 2)[2] = [(1 16)(2 20)(3 19)(4 17)(5 18)(6 14)(7 12)(8 15)(9 11)(10 13)]
( 3)[2] = [(1 17)(2 16)(3 20)(4 18)(5 19)(6 11)(7 13)(8 12)(9 15)(10 14)]
( 4)[2] = [(1 18)(2 17)(3 16)(4 19)(5 20)(6 15)(7 14)(8 13)(9 12)(10 11)]
( 5)[2] = [(1 19)(2 18)(3 17)(4 20)(5 16)(6 12)(7 11)(8 14)(9 13)(10 15)]
( 6)[2] = [(1 20)(2 19)(3 18)(4 16)(5 17)(6 13)(7 15)(8 11)(9 14)(10 12)]
( 7)[4] = [(1 6 17 11)(2 7 16 13)(3 9 20 15)(4 8 18 12)(5 10 19 14)]
( 8)[4] = [(1 7 18 14)(2 9 17 12)(3 10 16 11)(4 6 19 15)(5 8 20 13)]
( 9)[4] = [(1 8 16 15)(2 6 20 14)(3 7 19 12)(4 10 17 13)(5 9 18 11)]
(10)[4] = [(1 9 19 13)(2 10 18 15)(3 8 17 14)(4 7 20 11)(5 6 16 12)]
(11)[4] = [(1 10 20 12)(2 8 19 11)(3 6 18 13)(4 9 16 14)(5 7 17 15)]
(12)[4] = [(1 11 17 6)(2 13 16 7)(3 15 20 9)(4 12 18 8)(5 14 19 10)]
(13)[4] = [(1 12 20 10)(2 11 19 8)(3 13 18 6)(4 14 16 9)(5 15 17 7)]
(14)[4] = [(1 13 19 9)(2 15 18 10)(3 14 17 8)(4 11 20 7)(5 12 16 6)]
(15)[4] = [(1 14 18 7)(2 12 17 9)(3 11 16 10)(4 15 19 6)(5 13 20 8)]
(16)[4] = [(1 15 16 8)(2 14 20 6)(3 12 19 7)(4 13 17 10)(5 11 18 9)]
(17)[5] = [(1 2 3 5 4)(6 10 7 8 9)(11 15 12 13 14)(16 17 18 19 20)]
(18)[5] = [(1 3 4 2 5)(6 7 9 10 8)(11 12 14 15 13)(16 18 20 17 19)]
(19)[5] = [(1 4 5 3 2)(6 9 8 7 10)(11 14 13 12 15)(16 20 19 18 17)]
(20)[5] = [(1 5 2 4 3)(6 8 10 9 7)(11 13 15 14 12)(16 19 17 20 18)]

Tower 1
  |G14| = 1  => [Q(y):Q] = 20
  |G9| = 2   => [Q(c):Q] = 10 with c=y^11 + -22*y^6 + -4*y
  |G8| = 4   => [Q(l):Q] = 5 with l=y^16 + 144*y^11 + 38144*y^6 + 89536*y
  |G1| = 20  => [Q:Q]    = 1
Tower 2
  |G14| = 1  => [Q(y):Q] = 20
  |G10| = 2  => [Q(d):Q] = 10 with d=y^11 + 253*y^6 + -1654*y
  |G7| = 4   => [Q(k):Q] = 5 with k=y^16 + 144*y^11 + 38144*y^6 + -251464*y
  |G1| = 20  => [Q:Q]    = 1
Tower 3
  |G14| = 1  => [Q(y):Q] = 20
  |G11| = 2  => [Q(e):Q] = 10 with e=y^11 + 968*y^6 + 6156*y
  |G5| = 4   => [Q(i):Q] = 5 with i=y^16 + 1358/9*y^11 + 341932/9*y^6 + 241064*y
  |G1| = 20  => [Q:Q]    = 1
Tower 4
  |G14| = 1  => [Q(y):Q] = 20
  |G12| = 2  => [Q(f):Q] = 10 with f=y^11 + -572*y^6 + -24204*y
  |G6| = 4   => [Q(j):Q] = 5 with j=y^16 + 556/3*y^11 + 145804/3*y^6 + 1953504*y
  |G1| = 20  => [Q:Q]    = 1
Tower 5
  |G14| = 1  => [Q(y):Q] = 20
  |G13| = 2  => [Q(g):Q] = 10 with g=y^11 + 902/9*y^6 + 1096*y
  |G4| = 4   => [Q(h):Q] = 5 with h=y^16 + 162*y^11 + 39948*y^6 + 450264*y
  |G1| = 20  => [Q:Q]    = 1
Tower 6
  |G14| = 1  => [Q(y):Q] = 20
  |G9| = 2   => [Q(c):Q] = 10 with c=y^11 + -22*y^6 + -4*y
  |G2| = 10  => [Q(b):Q] = 2 with b=y^15 + -374/3*y^10 + -1412*y^5
  |G1| = 20  => [Q:Q]    = 1
Tower 7
  |G14| = 1  => [Q(y):Q] = 20
  |G10| = 2  => [Q(d):Q] = 10 with d=y^11 + 253*y^6 + -1654*y
  |G2| = 10  => [Q(b):Q] = 2 with b=y^15 + -374/3*y^10 + -1412*y^5
  |G1| = 20  => [Q:Q]    = 1
Tower 8
  |G14| = 1  => [Q(y):Q] = 20
  |G11| = 2  => [Q(e):Q] = 10 with e=y^11 + 968*y^6 + 6156*y
  |G2| = 10  => [Q(b):Q] = 2 with b=y^15 + -374/3*y^10 + -1412*y^5
  |G1| = 20  => [Q:Q]    = 1
Tower 9
  |G14| = 1  => [Q(y):Q] = 20
  |G12| = 2  => [Q(f):Q] = 10 with f=y^11 + -572*y^6 + -24204*y
  |G2| = 10  => [Q(b):Q] = 2 with b=y^15 + -374/3*y^10 + -1412*y^5
  |G1| = 20  => [Q:Q]    = 1
Tower 10
  |G14| = 1  => [Q(y):Q] = 20
  |G13| = 2  => [Q(g):Q] = 10 with g=y^11 + 902/9*y^6 + 1096*y
  |G2| = 10  => [Q(b):Q] = 2 with b=y^15 + -374/3*y^10 + -1412*y^5
  |G1| = 20  => [Q:Q]    = 1
Tower 11
  |G14| = 1  => [Q(y):Q] = 20
  |G3| = 5   => [Q(a):Q] = 4 with a=y^5
  |G2| = 10  => [Q(b):Q] = 2 with b=y^15 + -374/3*y^10 + -1412*y^5
  |G1| = 20  => [Q:Q]    = 1
*/