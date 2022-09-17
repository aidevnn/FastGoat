using FastGoat;

var z2 = new Zn(2);
var z3 = new Zn(3);
var z4 = new Zn(4);
var z5 = new Zn(5);
var z6 = new Zn(6);
var z7 = new Zn(7);
var z8 = new Zn(8);
var z9 = new Zn(9);
var z10 = new Zn(10);
var z12 = new Zn(12);
var z15 = new Zn(15);
var z18 = new Zn(18);
var z20 = new Zn(20);
var z30 = new Zn(30);

var s2 = new Sn(2);
var s3 = new Sn(3);
var s4 = new Sn(4);
var s5 = new Sn(5);
var s6 = new Sn(6);
var s7 = new Sn(7);

// {
//     var sn = new Sn(4);
//     Perm p0 = (sn, (1, 2));
//     Perm p1 = (sn, (3, 4));
//     var wgr = sn.Generate(p0, p1);
//     wgr.DisplayDetails();
// }

// {
//     var z2xz2 = Group.CartesianProduct(z2, z2);
//     var wgr = z2xz2.Generate(z2xz2[0, 1], z2xz2[1, 0]);
//     wgr.DisplayDetails();
// }

// {
//     var wgr = z4.Generate(z4[1]);
//     wgr.DisplayDetails();
// }

// {
//     var sn = new Sn(3);
//     Perm p0 = (sn, (1, 2));
//     Perm p1 = (sn, (1, 2, 3));
//     var wgr = sn.Generate(p0, p1);
//     wgr.DisplayDetails();
// }

// {
//     var z12 = new Zn(12);
//     var wgr = z12.Generate(z12[3], z12[4]);
//     wgr.DisplayDetails();
// }

// {
//     var z4xz3 = Group.CartesianProduct(z4, z3);
//     var wgr = z4xz3.Generate(z4xz3[1, 0], z4xz3[0, 1]);
//     wgr.DisplayDetails();
// }

// {
//     var z6xz2 = Group.CartesianProduct(z6, z2);
//     var wgr = z6xz2.Generate(z6xz2[1, 0], z6xz2[0, 1]);
//     wgr.DisplayDetails();
// }

// {
//     var z3xz2xz2 = Group.CartesianProduct(z3, z2, z2);
//     var wgr = z3xz2xz2.Generate(z3xz2xz2[1, 0, 0], z3xz2xz2[0, 1, 0], z3xz2xz2[0, 0, 1]);
//     wgr.DisplayDetails();
// }

// {
//     var sn = new Sn(4);
//     Perm p0 = (sn, (1, 2));
//     Perm p1 = (sn, (3, 4));
//     Perm p2 = (sn, (1, 2, 3));
//     var wgr = sn.Generate(p0 * p1, p2);
//     wgr.DisplayDetails();
// }

// {
//     var G = z20.Generate(z20[1]);
//     var H = z20.Generate(z20[4]);
//     G.DisplayDetails(SortElements.ByValue);
//     H.DisplayDetails(SortElements.ByValue);
//     var G_H = G.Over(H);
//     G_H.DisplayDetails(SortElements.ByValue);
//     G_H.DisplayCosets();
// }

// {
//     var sn = new Sn(4);
//     Perm p0 = (sn, (1, 2));
//     Perm p1 = (sn, (3, 4));
//     Perm p2 = (sn, (1, 3));
//     Perm p3 = (sn, (2, 4));
//     Perm c3 = (sn, (1, 2, 3));
//     var G = sn.Generate(p0 * p1, c3);
//     var H = sn.Generate(p0 * p1, p2 * p3);
//     G.DisplayDetails();
//     H.DisplayDetails();
//     var G_H = G.Over(H);
//     G_H.DisplayDetails();
//     G_H.DisplayCosets();
// }

// {
//     var z8x18x30 = Group.CartesianProduct(z8, z18, z30);
//     var G = z8x18x30.Generate(z8x18x30[1, 0, 0], z8x18x30[0, 1, 0], z8x18x30[0, 0, 1]);
//     G.DisplayHead();
//     var eMax = G.OrderByDescending(G.GetOrderOf).ThenAscending().First();
//     G.DisplayElement(eMax, "eMax");
//     var H = z8x18x30.Generate(z8x18x30[1]);
//     H.DisplayHead();
//     var G_H = G.Over(H);
//     G_H.DisplayHead();
// }

// {
//     var z6x8x18x20 = Group.CartesianProduct(z4, z6, z8, z10);
//     var wgr = z6x8x18x20.Generate(z6x8x18x20[1, 0, 0, 0], z6x8x18x20[0, 1, 0, 0], z6x8x18x20[0, 0, 1, 0], z6x8x18x20[0, 0, 0, 1]);
//     wgr.InvariantFactors();
// }

// {
//     var z8x18x30 = Group.CartesianProduct(z8, z18, z30);
//     var G = z8x18x30.Generate(z8x18x30[1, 0, 0], z8x18x30[0, 1, 0], z8x18x30[0, 0, 1]);
//     G.InvariantFactors();
// }

// {
//     var sn = new Sn(6);
//     Perm p0 = (sn, (1, 2));
//     Perm p1 = (sn, (3, 4, 5, 6));
//     var G = sn.Generate(p0, p1);
//     G.InvariantFactors();
// }

// {
//     var s6 = new Sn(6);
//     var s6xs6 = Group.CartesianProduct(s6, s6);
//     Perm id = s6.Neutral();
//     Perm p0 = (s6, (1, 2));
//     Perm p1 = (s6, (3, 4, 5));
//     Perm p2 = (s6, (3, 4, 5, 6));
//     var G = s6xs6.Generate((p0 * p1, id), (id, p2));
//     G.InvariantFactors();
//     G.DisplayDetails();
// }

{
    // Creating the BaseGroup generator. These classes are created by users
    // by implementing the interface IGroup. They dont contains data structures
    // for storing the elements.
    // var z6 = new Zn(6);
    // var z10 = new Zn(10);
    // var z18 = new Zn(18);

    // Implicit definition, the tuple will be processed to a group element.
    // Base Elements classes are created by users by implementing the
    // interface IElt and they also contain a property named Group
    // for the BaseGroup they belong to.
    // ZnInt e1 = (z6, 1);
    // Console.WriteLine(e1);

    // Dotnet Tuple elements, nothing will be done. They keyword var
    // will not interprete the tuple as a Group element
    // var v1 = (z6, 1);
    // Console.WriteLine("{0} != {1}", e1, v1);

    // G[i] is another way to create the element i from the group G,
    // but it can be hard to predict its value depending of the group,
    // Console.WriteLine("{0} == {1}", e1, z6[1]);

    // ZnInt e2 = (z10, 2);
    // // Tuple of integers from Z6 x Z10
    // Ep<ZnInt, ZnInt> ep = (e1, e2);
    // Console.WriteLine("Tuple {0} of {1}", ep, ep.Group);
    // Console.WriteLine();
}

// {
//     ZnInt e1 = (z6, 2);
//     ZnInt e2 = (z18, 3);
//     Group.Generate(z6[2]).DisplayDetails();
//     Group.Generate(z18[3]).DisplayDetails();
//     Ep<ZnInt, ZnInt> ep = (z6[2], z18[3]);
//     Group.Generate(ep).DisplayDetails();

//     Group.Generate<Ep<ZnInt, ZnInt>>(((z6, 2), (z18, 0)), ((z6, 0), (z18, 3))).DisplayDetails();
// }

// {
//     GlobalStopWatch.Restart();
//     Perm c7 = (s7, (1, 2, 3, 4, 5, 6, 7));
//     Perm c2 = (s7, (1, 2));
//     var S7 = Group.Generate(c2, c7);

//     var allC3 = S7.Where(e => S7.GetOrderOf(e) == 3);
//     Console.WriteLine("|S7|={0}, |{{K in S7 with K~C3}}| = {1}", S7.Count(), allC3.Count());
//     Console.WriteLine();

//     var sc3 = allC3.First(c3 => Group.Generate(c7, c3).Count() == 21);
//     Console.WriteLine("First Solution |HK| = 21 : h = {0} and k = {1}", c7, sc3);

//     Console.WriteLine();

//     var H = Group.Generate(c7);
//     var K = Group.Generate(sc3);
//     var G = Group.Generate(c7, sc3);
//     var GoH = G.Over(H);

//     H.DisplayDetails("H");
//     K.DisplayDetails("K");
//     G.DisplayDetails("G=HK");
//     GoH.DisplayDetails("G/H");
//     GoH.DisplayCosets();
//     GlobalStopWatch.Show("Example");
// }

// {
//     GlobalStopWatch.Restart();
//     var s7 = new Sn(7);
//     Perm c7 = (s7, (1, 2, 3, 4, 5, 6, 7));
//     Perm c2 = (s7, (1, 2));
//     var S7 = Group.Generate(c2, c7);

//     var allC3 = S7.Where(e => S7.GetOrderOf(e) == 3);
//     Console.WriteLine("|S7|={0}, |{{b in S7 with b^3 = 1}}| = {1}", S7.Count(), allC3.Count());
//     Console.WriteLine();

//     var sc3 = allC3.First(c3 => (c7 ^ 2) * c3 == c3 * c7);
//     var allSols = allC3.Count(c3 => (c7 ^ 2) * c3 == c3 * c7);
//     Console.WriteLine("First Solution a^7 = b^3 = 1 and a^2 * b = b * a : a = {0} and b = {1}", c7, sc3);
//     Console.WriteLine("All Solutions : {0}", allSols);

//     Console.WriteLine();

//     var H = Group.Generate(c7);
//     var K = Group.Generate(sc3);
//     var G = Group.Generate(c7, sc3);
//     var GoH = G.Over(H);

//     H.DisplayDetails("H");
//     K.DisplayDetails("K");
//     G.DisplayDetails("G=HK");
//     GoH.DisplayDetails("G/H");
//     GoH.DisplayCosets();
//     GlobalStopWatch.Show("Example");
// }

// {
//     GlobalStopWatch.Restart();
//     var s7 = new Sn(7);
//     Perm c7 = (s7, (1, 2, 3, 4, 5, 6, 7));
//     Perm c2 = (s7, (1, 2));
//     var S7 = Group.Generate(c2, c7);

//     var allC2 = S7.Where(e => S7.GetOrderOf(e) == 2);
//     var allC3 = S7.Where(e => S7.GetOrderOf(e) == 3);
//     var set = allC3.SelectMany(a => allC2.Select(b => (a, b)));
//     var ne = S7.Neutral();
//     var filter1 = from e in set
//                   where ((e.a * e.b) ^ 2) != ne && ((e.a * e.b) ^ 4) == ne
//                   select e;

//     var filter2 = from e in filter1
//                   from ec in allC2
//                   where ((e.a * ec) ^ 2) == ne
//                   && ((e.b * ec) ^ 3) == ne
//                   && !Group.Generate(e.a, e.b).Contains(ec)
//                   select (e.a, e.b, ec);

//     var (a, b, c) = filter2.First();

//     Console.WriteLine("All C2 = {0}", allC2.Count());
//     Console.WriteLine("All C3 = {0}", allC3.Count());
//     Console.WriteLine("|S7|={0}, |{{(a,b) in S7xS7 with a^3 = b^2 = 1}}| = {1}", S7.Count(), set.Count());
//     Console.WriteLine("Filter (ab)^4 = 1. Count = {0}", filter1.Count());
//     Console.WriteLine("First Solution (ac)^2 = (bc)^3 = 1 : a = {0} and b = {1} and c = {2}", a, b, c);
//     // Console.WriteLine();

//     var G = Group.Generate(a, b, c);
//     G.DisplayHead("G");

//     GlobalStopWatch.Show("Example");
// }

// {
//     ZnInt cn = (z4, 1);
//     ZnInt cg = (z2, 1);
//     var Cg = Group.Generate(cg);
//     var Cn = Group.Generate(cn);
//     Func<ZnInt, ZnInt, ZnInt> action = (ZnInt g, ZnInt x) => g == Cg.Neutral() ? x : x ^ -1;
//     var Cn_sp_Cg = new SemiDirectProduct<ZnInt, ZnInt>(Cn, Cg, action);
//     Cn_sp_Cg.DisplayDetails();
//     Cn_sp_Cg.Over(Cn_sp_Cg.Ncan).DisplayDetails();
// }

// {
//     GlobalStopWatch.Restart();
//     ZnInt cn = (z7, 1);
//     ZnInt cg = (z3, 1);
//     var Cg = Group.Generate(cg);
//     var Cn = Group.Generate(cn);
//     var Cn_sp_Cg = Group.SemiDirectProd(Cn, Cg, (g, x) => g == (z3, 0) ? x : g == cg ? x ^ 2 : x ^ 4);
//     Cn_sp_Cg.DisplayDetails("G=HK");
//     var Cn_sp_CgoNcan = Cn_sp_Cg.Over(Cn_sp_Cg.Ncan);
//     Cn_sp_CgoNcan.DisplayDetails("G/H");
//     Cn_sp_CgoNcan.DisplayCosets();
//     GlobalStopWatch.Show("SemiDirectProduct");
// }

// {
//     Perm cn = (s7, (1, 2, 3, 4, 5, 6, 7));
//     Perm cg = (s3, (1, 2, 3));
//     var Cg = Group.Generate(cg);
//     var Cn = Group.Generate(cn);
//     Func<Perm, Perm, Perm> action = (Perm g, Perm x) => g == Cg.Neutral() ? x : g == cg ? x ^ 2 : x ^ 4;
//     var Cn_sp_Cg = new SemiDirectProduct<Perm, Perm>(Cn, Cg, action);
//     Cn_sp_Cg.DisplayDetails();
//     Cn_sp_Cg.Over(Cn_sp_Cg.Ncan).DisplayDetails();
// }

// {
//     var z0 = z12;
//     var z1 = z8;
//     var zx = Group.CartesianProduct(z0, z1);
//     ZnInt c0 = (z0, 1);
//     ZnInt c1 = (z1, 1);
//     var C0 = Group.Generate(c0);
//     var C1 = Group.Generate(c1);
//     var C2 = Group.Generate(zx[0, 1], zx[1, 0]);
//     var C2sg = Group.Generate(zx[0, 2], zx[3, 0]).ToList();
//     C2sg.ForEach(e => Console.WriteLine(e));
//     Console.WriteLine("IsSubGroup : {0}", C2.VerifySubGroup(C2sg));
//     Group.Homomorphism(C0, C0, z0[1], z0[2]).DisplayFullGraph();
//     Group.Homomorphism(C0, C0, z0[1], z0[3]).DisplayFullGraph();
//     Group.Homomorphism(C0, C1, z0[1], z1[3]).DisplayFullGraph();
//     Group.Homomorphism(C1, C0, z1[1], z0[1]).DisplayFullGraph();
//     Group.Homomorphism(C1, C0, z1[1], z0[2]).DisplayFullGraph();
//     Group.Homomorphism(C2, C1, (zx[0, 1], z1[1]), (zx[1, 0], z1[0])).DisplayFullGraph();
//     Group.Homomorphism(C1, C1, (z1[1], z1[1])).DisplayFullGraph();
//     Group.Homomorphism(C1, C1, (z1[1], z1[2])).DisplayFullGraph();
//     Group.Homomorphism(C1, C1, (z1[1], z1[4])).DisplayFullGraph();
//     Group.Homomorphism(C2, C1, (zx[0, 1], z1[1]), (zx[1, 0], z1[2])).DisplayFullGraph();
// }

// {
//     GlobalStopWatch.Restart();
//     var z7xz3 = Group.CartesianProduct(z7, z3);
//     var Cg = Group.Generate(z7xz3[0, 1]);
//     var Cn = Group.Generate(z7xz3[1, 0]);
//     var Cn_sp_Cg = Group.SemiDirectProd(Cn, Cg, (g, x) => z7xz3.Op(z7xz3.Op(g, x), z7xz3.Invert(g)));
//     Cn_sp_Cg.DisplayDetails("G=HK");
//     var Cn_sp_CgoNcan = Cn_sp_Cg.Over(Cn_sp_Cg.Ncan);
//     Cn_sp_CgoNcan.DisplayDetails("G/H");
//     Cn_sp_CgoNcan.DisplayCosets();
//     GlobalStopWatch.Show("SemiDirectProduct");
// }

// {
//     var zx = z7;
//     var zg = z2;
//     var n = (zx.mod - 1) / (zg.mod - 1);
//     var str = Enumerable.Range(0, zg.mod).Select(i => i == 0 ? zx[1] : zx[1] ^ (n * i)).Glue(" | ");
//     Console.WriteLine(str);
// }

// for (int x = 3; x < 10; ++x)
// {
//     var zx = new Zn(x);
//     for (int g = 2; g < x; ++g)
//     {
//         var zg = new Zn(g);
//         var str = Enumerable.Range(0, zg.mod).Select(i => i == 0 ? zx[1] : zx[1] ^ -i).Glue(" | ", "{0,3}");
//         Console.WriteLine("{0,5} {1,5}; {2}", zx, zg, str);
//     }
//     Console.WriteLine();
// }

{
    GlobalStopWatch.Restart();
    var zx = z4;
    var zg = z2;
    var Cn = Group.Generate(zx[1]);
    var Cg = Group.Generate(zg[1]);
    var Cn_sp_Cg = Group.SemiDirectProd(Cn, Cg);
    Cn_sp_Cg.DisplayDetails("G = H ⋊  K");
    var Cn_sp_CgoNcan = Cn_sp_Cg.Over(Cn_sp_Cg.Ncan);
    Cn_sp_CgoNcan.DisplayDetails("G/H");
    Cn_sp_CgoNcan.DisplayCosets();
    GlobalStopWatch.Show("SemiDirectProduct");
}

{
    GlobalStopWatch.Restart();
    var zx = new Zn(7);
    var zg = z3;
    var Cn = Group.Generate(zx[1]);
    var Cg = Group.Generate(zg[1]);
    var Cn_sp_Cg = Group.SemiDirectProd(Cn, Cg);
    Cn_sp_Cg.DisplayDetails("G = H ⋊  K");
    var Cn_sp_CgoNcan = Cn_sp_Cg.Over(Cn_sp_Cg.Ncan);
    Cn_sp_CgoNcan.DisplayDetails("G/H");
    Cn_sp_CgoNcan.DisplayCosets();
    GlobalStopWatch.Show("SemiDirectProduct");
}

{
    GlobalStopWatch.Restart();
    var zx = z5;
    var zg = z4;
    var Cn = Group.Generate(zx[1]);
    var Cg = Group.Generate(zg[1]);
    var Cn_sp_Cg = Group.SemiDirectProd(Cn, Cg);
    Cn_sp_Cg.DisplayDetails("G = H ⋊  K");
    var Cn_sp_CgoNcan = Cn_sp_Cg.Over(Cn_sp_Cg.Ncan);
    Cn_sp_CgoNcan.DisplayDetails("G/H");
    Cn_sp_CgoNcan.DisplayCosets();
    GlobalStopWatch.Show("SemiDirectProduct");
}
