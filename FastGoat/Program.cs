using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void MatrixFormTinyGroups(int maxOrder)
{
    var total = 0;
    var missing = new List<ANameElt>();
    foreach (var g in maxOrder.Range(1).SelectMany(o => FG.AllGroupsOfOrder(o)))
    {
        var found = false;
        ++total;
        var gSubgrs = g.AllSubgroups().ToGroupWrapper();
        var names = NamesTree.BuildName(gSubgrs);
        if (names[0] is Leaf leaf)
        {
            var (pr, coefs) = leaf.LeafDetails();
            if (pr == "C")
            {
                found = true;
                var abMat = FG.AbelianMat(coefs);
                FG.DisplayBox(gSubgrs, -1, false, false, 16);
                DisplayGroup.Generators(abMat);
                names.Println("Group names");

                Console.WriteLine();
            }
            else if (pr == "Q" || pr == "Dic")
            {
                found = true;
                var m = pr == "Dic" ? coefs[0] : coefs[0] / 4;
                var dicMat = FG.DicyclicGL2p(m);
                FG.DisplayBox(gSubgrs, -1, false, false, 16);
                DisplayGroup.Generators(dicMat);
                names.Println("Group names");

                Console.WriteLine();
            }
            else if (pr == "A" || pr == "S")
            {
                if (coefs[0] == 4)
                {
                    found = true;
                    var gl = new GL(3, 3);
                    var (a, b) = pr == "S"
                        ? (gl[0, 1, 0, 0, 0, 1, 1, 0, 0], gl[1, 0, 0, 0, 0, 1, 0, 2, 0])
                        : (gl[0, 1, 0, 0, 0, 1, 1, 0, 0], gl[0, 0, 1, 2, 0, 0, 0, 2, 0]);
                    var gMat = Group.Generate(g.Name, gl, a, b);
                    FG.DisplayBox(gSubgrs, -1, false, false, 16);
                    DisplayGroup.Generators(gMat);
                    names.Println("Group names");

                    Console.WriteLine();
                }
            }
            else if (pr == "GL" || pr == "SL")
            {
                found = true;
                var gMat = pr == "GL" ? FG.GLnp(coefs[0], coefs[1]) : FG.SLnp(coefs[0], coefs[1]);
                FG.DisplayBox(gSubgrs, -1, false, false, 16);
                DisplayGroup.Generators(gMat);
                names.Println("Group names");

                Console.WriteLine();
            }
        }
        else if (names[0] is SemiDirectProductOp sdpOp)
        {
            var mtCyc = sdpOp.MetaCyclicDetails();
            if (mtCyc.Length == 3)
            {
                found = true;
                var mtMat = FG.MetaCyclicSdpMat(mtCyc[0], mtCyc[1], mtCyc[2]);
                FG.DisplayBox(gSubgrs, -1, false, false, 16);
                DisplayGroup.Generators(mtMat);
                names.Println("Group names");

                Console.WriteLine();
            }
        }
        
        if(!found)
            missing.Add(names[0]);
    }

    Console.WriteLine($"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    missing.Where(e=>e.ContentType == ANameElt.NodeType.DirectProduct).Println("Missing Direct Products");
}

{
    MatrixFormTinyGroups(24);
}