using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
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
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    FG.Abelian(8).ToCGW().AllSubgroups().DisplayLattice();
    FG.Abelian(4, 2).ToCGW().AllSubgroups().DisplayLattice();
    FG.Abelian(2, 2, 2).ToCGW().AllSubgroups().DisplayLattice();
    FG.Dihedral(4).ToCGW().AllSubgroups().DisplayLattice();
    FG.Quaternion(8).ToCGW().AllSubgroups().DisplayLattice();
    FG.Abelian(12).ToCGW().AllSubgroups().DisplayLattice();
    FG.Abelian(6, 2).ToCGW().AllSubgroups().DisplayLattice();
    FG.Alternate(4).ToCGW().AllSubgroups().DisplayLattice();
    FG.Dihedral(6).ToCGW().AllSubgroups().DisplayLattice();
    FG.DiCyclic(3).ToCGW().AllSubgroups().DisplayLattice();
    FG.Abelian(8, 2).ToCGW().AllSubgroups().DisplayLattice();
    
    FG.Quaternion(16).ToCGW().AllSubgroups().DisplayLattice();
    FG.SL2p(3).ToCGW().AllSubgroups().DisplayLattice();
    FG.GL2p(3).ToCGW().AllSubgroups().DisplayLattice();
    FG.Alternate(5).ToCGW().AllSubgroups().DisplayLattice();
    FG.Symmetric(5).ToCGW().AllSubgroups().DisplayLattice();
}

/* 
   |SL(2,3)| = 24
   SubGroupsInfos { AllSubGr = 15, AllConjsCl = 7, AllNorms = 4 }
   
   Lattice Subgroups
       SL(2,3) --> C6 -------> C2 -------> C1      
       SL(2,3) --> C6 -------> C3 -------> C1      
       SL(2,3) --> Q8 -------> C4 -------> C2 -------> C1      
   Total:3
   
   Chief Series
       C1       ⊲  C2       ⊲  Q8       ⊲  SL(2,3) 
   Composition Series
       C1       ⊲  C2       ⊲  C4       ⊲  Q8       ⊲  SL(2,3) 
   
   |GL(2,3)| = 48
   SubGroupsInfos { AllSubGr = 55, AllConjsCl = 16, AllNorms = 5 }
   
   Lattice Subgroups
       GL(2,3) --> D12 ------> C2 x C2 --> C2₁ ------> C1      
       GL(2,3) --> D12 ------> C2 x C2 --> C2₂ ------> C1      
       GL(2,3) --> D12 ------> C6 -------> C2₁ ------> C1      
       GL(2,3) --> D12 ------> C6 -------> C3 -------> C1      
       GL(2,3) --> D12 ------> S3₁ ------> C2₂ ------> C1      
       GL(2,3) --> D12 ------> S3₁ ------> C3 -------> C1      
       GL(2,3) --> D12 ------> S3₂ ------> C2₂ ------> C1      
       GL(2,3) --> D12 ------> S3₂ ------> C3 -------> C1      
       GL(2,3) --> SL(2,3) --> C6 -------> C2₁ ------> C1      
       GL(2,3) --> SL(2,3) --> C6 -------> C3 -------> C1      
       GL(2,3) --> QD16 -----> Q8 -------> C4 -------> C2₁ ------> C1      
       GL(2,3) --> QD16 -----> C8 -------> C4 -------> C2₁ ------> C1      
       GL(2,3) --> QD16 -----> D8 -------> C4 -------> C2₁ ------> C1      
       GL(2,3) --> QD16 -----> D8 -------> C2 x C2 --> C2₁ ------> C1      
       GL(2,3) --> QD16 -----> D8 -------> C2 x C2 --> C2₂ ------> C1      
       GL(2,3) --> SL(2,3) --> Q8 -------> C4 -------> C2₁ ------> C1      
   Total:16
   
   Chief Series
       C1       ⊲  C2       ⊲  Q8       ⊲  SL(2,3)  ⊲  GL(2,3) 
   Composition Series
       C1       ⊲  C2       ⊲  C4       ⊲  Q8       ⊲  SL(2,3)  ⊲  GL(2,3) 
*/