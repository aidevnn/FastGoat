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
    var g = Product.Generate(new Cn(4), new Cn(8));
    var allHoms = Group.AllHomomorphisms(g, g);
    int k = 1;
    HashSet<int> set = new HashSet<int>();
    List<ConcreteGroup<Ep2<ZnInt, ZnInt>>> allSubGroups = new List<ConcreteGroup<Ep2<ZnInt, ZnInt>>>();
    foreach (var kp in allHoms.OrderBy(e => e.Values.Distinct().Count()))
    {
        var sub = kp.Values.Distinct().Ascending().ToArray();
        var hash = sub.Select(e => g.ElementsOrders[e]).Aggregate(1, (acc, a) => (acc, a).GetHashCode());
        if (!set.Add(hash))
            continue;

        var subGroup = Group.Generate($"S{k++}", g, sub);
        allSubGroups.Add(subGroup);
    }

    foreach (var subGroup in allSubGroups)
    {
        DisplayGroup.Head(subGroup);
        Console.WriteLine(subGroup.GetGenerators().Glue(", "));
    }

    Console.WriteLine($"Nb SubGroups = {k}");
}