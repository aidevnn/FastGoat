using FastGoat;
using System.Diagnostics;

var sw = Stopwatch.StartNew();

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

{
    var s6 = new Sn(6);
    var s6xs6 = Group.CartesianProduct(s6, s6);
    Perm id = s6.Neutral();
    Perm p0 = (s6, (1, 2));
    Perm p1 = (s6, (3, 4, 5));
    Perm p2 = (s6, (3, 4, 5, 6));
    var G = s6xs6.Generate((p0 * p1, id), (id, p2));
    G.InvariantFactors();
    G.DisplayDetails();
}
