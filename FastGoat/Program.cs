using FastGoat;
using FastGoat.Gp;
using FastGoat.UserGroup;
using static FastGoat.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    for (int k = 2; k < 32; ++k)
    {
        Console.WriteLine($"############# U{k,-2} #############");
        DisplayGroup.HeadElements(new Un(k));
    }
}
