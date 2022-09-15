using FastGoat;

var z2 = new Zn(2);
var z3 = new Zn(3);
var z4 = new Zn(4);
var z5 = new Zn(5);
var z6 = new Zn(6);
var z8 = new Zn(8);
var z10 = new Zn(10);
var z15 = new Zn(15);
var z18 = new Zn(18);
var z20 = new Zn(20);
var z30 = new Zn(30);

var s2 = new Sn(2);
var s3 = new Sn(3);
var s4 = new Sn(4);
var s5 = new Sn(5);
var s6 = new Sn(6);

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

{
    GlobalStopWatch.Restart();
    var s7 = new Sn(7);
    Perm c7 = (s7, (1, 2, 3, 4, 5, 6, 7));
    Perm c2 = (s7, (1, 2));
    var S7 = Group.Generate(c2, c7);

    var allC7 = S7.Where(e => S7.GetOrderOf(e) == 7);
    var allC3 = S7.Where(e => S7.GetOrderOf(e) == 3);
    var set = allC7.SelectMany(h => allC3.Select(k => (h, k)));
    Console.WriteLine("|S7|={0}, |{{HK in S7 with H~C7 and K~C3}}| = {1}", S7.Count(), set.Count());
    Console.WriteLine();

    var s = set.First(e => Group.Generate(e.h, e.k).Count() == 21);
    Console.WriteLine("First Solution |HK| = 21 : h = {0} and k = {1}", s.h, s.k);

    Console.WriteLine();

    var H = Group.Generate(s.h);
    var K = Group.Generate(s.k);
    var G = Group.Generate(s.h, s.k);
    var GoH = G.Over(H);

    H.DisplayDetails("H");
    K.DisplayDetails("K");
    G.DisplayDetails("G=HK");
    GoH.DisplayDetails("G/H");
    GoH.DisplayCosets();

    Console.WriteLine();
    GlobalStopWatch.Show("Example");
}

{
    GlobalStopWatch.Restart();
    var s7 = new Sn(7);
    Perm c7 = (s7, (1, 2, 3, 4, 5, 6, 7));
    Perm c2 = (s7, (1, 2));
    var S7 = Group.Generate(c2, c7);

    var allC7 = S7.Where(e => S7.GetOrderOf(e) == 7);
    var allC3 = S7.Where(e => S7.GetOrderOf(e) == 3);
    var set = allC7.SelectMany(h => allC3.Select(k => (h, k)));
    Console.WriteLine("|S7|={0}, |{{HK in S7 with H~C7 and K~C3}}| = {1}", S7.Count(), set.Count());
    Console.WriteLine();

    var s = set.First(e => Group.Generate(e.h, e.k).Count() == 21);
    Console.WriteLine("First Solution |HK| = 21 : h = {0} and k = {1}", s.h, s.k);

    Console.WriteLine();

    var H = Group.Generate(s.h);
    var K = Group.Generate(s.k);
    var G = Group.Generate(s.h, s.k);
    var GoH = G.Over(H);

    H.DisplayDetails("H");
    K.DisplayDetails("K");
    G.DisplayDetails("G=HK");
    GoH.DisplayDetails("G/H");
    GoH.DisplayCosets();

    Console.WriteLine();
    GlobalStopWatch.Show("Example");
}

{
    GlobalStopWatch.Restart();
    var s7 = new Sn(7);
    Perm c7 = (s7, (1, 2, 3, 4, 5, 6, 7));
    Perm c2 = (s7, (1, 2));
    var S7 = Group.Generate(c2, c7);

    var allC7 = S7.Where(e => S7.GetOrderOf(e) == 7);
    var allC3 = S7.Where(e => S7.GetOrderOf(e) == 3);
    var set = allC7.SelectMany(h => allC3.Select(k => (h, k)));
    Console.WriteLine("|S7|={0}, |{{HK in S7 with H~C7 and K~C3}}| = {1}", S7.Count(), set.Count());
    Console.WriteLine();

    var s = set.First(e => Group.Generate(e.h, e.k).Count() == 21);
    Console.WriteLine("First Solution |HK| = 21 : h = {0} and k = {1}", s.h, s.k);

    Console.WriteLine();

    var H = Group.Generate(s.h);
    var K = Group.Generate(s.k);
    var G = Group.Generate(s.h, s.k);
    var GoH = G.Over(H);

    H.DisplayDetails("H");
    K.DisplayDetails("K");
    G.DisplayDetails("G=HK");
    GoH.DisplayDetails("G/H");
    GoH.DisplayCosets();

    Console.WriteLine();
    GlobalStopWatch.Show("Example");
}
