using FastGoat;
using FastGoat.Gp;
using FastGoat.UserGroup;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

var n = Group.Generate("N", new Zn(9)[1]);
var g = Group.Generate("G", new Zn(3)[1]);
var sdp = Group.SemiDirectProd(n, g);
DisplayGroup.HeadElementsTableSdp(sdp);