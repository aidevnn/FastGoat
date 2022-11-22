using System.Diagnostics;
using System.Runtime.Intrinsics.X86;
using System.Threading.Channels;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

ConcreteGroup<Perm> PermGroup(string name, int n, params ValueType[] generators)
{
    var sn = new Sn(n);
    return Group.Generate(name, sn, generators.Select(v => sn.ComposesCycles(Tuple2Array.ComplexTuples(v))).ToArray());
}

IEnumerable<ConcreteGroup<Perm>> PermGroups()
{
    // 0
    yield return PermGroup("C4", 4, (1, 2, 3, 4));
    yield return PermGroup("S4", 4, (1, 2, 3, 4), (1, 2));
    yield return PermGroup("D8", 4, (1, 2, 3, 4), (1, 3));
    yield return PermGroup("V", 4, ((1, 2), (3, 4)), ((1, 3), (2, 4)));
    yield return PermGroup("A4", 4, (1, 2, 3), ((1, 2), (3, 4)));

    // 5
    yield return PermGroup("C5", 5, (1, 2, 3, 4, 5));
    yield return PermGroup("S5", 5, (1, 2, 3, 4, 5), (1, 2));
    yield return PermGroup("D10", 5, (1, 2, 3, 4, 5), ((2, 5), (3, 4)));
    yield return PermGroup("A5", 5, (1, 2, 3, 4, 5), (1, 2, 3));
    yield return PermGroup("C5 : C4", 5, (1, 2, 3, 4, 5), (2, 3, 5, 4));
    
    // 10
    yield return PermGroup("C6", 6, (1, 2, 3, 4, 5, 6));
    yield return PermGroup("H", 6, ((1, 5), (2, 4), (3, 6)), ((1, 6), (2, 5), (3, 4)));
    yield return PermGroup("D12", 6, (1, 2, 3, 4, 5, 6), ((2, 6), (3, 5)));
    yield return PermGroup("A4", 6, ((1, 3, 5), (2, 4, 6)), ((1, 2), (5, 6)));
    yield return PermGroup("S3 x C3", 6, (1, 2, 3), ((1, 4), (2, 5), (3, 6)));
    yield return PermGroup("A4 x C2 ", 6, ((1, 3, 5), (2, 4, 6)), (1, 2));
    yield return PermGroup("S4a", 6, ((1, 3, 5), (2, 4, 6)), ((1, 6), (2, 5)));
    yield return PermGroup("S4b", 6, ((1, 2), (3, 4), (5, 6)), ((1, 2, 3), (4, 5, 6)));
    yield return PermGroup("S3 x S3", 6, (1, 2, 3, 4, 5, 6), ((1, 3), (2, 4)));
    yield return PermGroup("(C3 x C3) : C4", 6, (1, 2, 3), ((1, 5, 2, 4), (3, 6)));
    yield return PermGroup("H", 6, (1, 2, 3, 4), ((1, 5), (3, 6)));
    yield return PermGroup("H", 6, (1, 2, 3, 4, 5), ((1, 6), (2, 5)));
    yield return PermGroup("H", 6, (1, 2, 3, 4, 5, 6), (1, 3));
    yield return PermGroup("H", 6, (1, 2, 3, 4, 5), ((1, 6), (2, 3), (4, 5)));
    yield return PermGroup("H", 6, (1, 2, 3), (1, 2, 4), (1, 2, 5), (1, 2, 6));
    yield return PermGroup("H", 6, (1, 2, 3, 4, 5, 6), (4, 5, 6));
    yield return PermGroup("H", 6, (1, 2, 3, 4, 5, 6), (1, 2));

    // 27
    yield return PermGroup("H", 7, (1, 2, 3, 4, 5, 6, 7));
    yield return PermGroup("H", 7, (1, 2, 3, 4, 5, 6, 7), ((2, 7), (3, 6), (4, 5)));
    yield return PermGroup("H", 7, (1, 2, 3, 4, 5, 6, 7), ((2, 3, 5), (4, 7, 6)));
    yield return PermGroup("H", 7, (1, 2, 3, 4, 5, 6, 7), (2, 4, 3, 7, 5, 6));
    yield return PermGroup("H", 7, (1, 2, 3, 4, 5, 6, 7), ((2, 3), (4, 7)));
    yield return PermGroup("H", 7, (1, 2, 3, 4, 5, 6, 7), (1, 2, 3));
    yield return PermGroup("H", 7, (1, 2, 3), (1, 2, 4), (1, 2, 5), (1, 2, 6), (1, 2, 7));
    yield return PermGroup("H", 7, (1, 2, 3, 4, 5, 6, 7), (1, 2));
    
    // 35
    yield return PermGroup("M11", 11, (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), ((3, 7, 11, 8), (4, 10, 5, 6)));
    yield return PermGroup("M12", 12, (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), ((3, 7, 11, 8), (4, 10, 5, 6)),
        ((1, 12), (2, 11), (3, 6), (4, 8), (5, 9), (7, 10)));
}

Xi ImgPerm(Perm p, Xi i) =>
    p.Sn.N > i.xi - 'a'
        ? new Xi(p.Table[i.xi - 'a'])
        : throw new IndexOutOfRangeException();

Ep<Xi> ImgPermArr(Perm p, Ep<Xi> xk) => new(xk.Ei.Select(i => ImgPerm(p, i)).ToArray());

void KTransitivity(ConcreteGroup<Perm> g, int k, int n)
{
    var set = n.Range().Select(i => new Xi(i)).MultiLoop(k)
        .Where(e => e.Distinct().Count() == k)
        .Select(e => new Ep<Xi>(e.ToArray()))
        .ToArray();

    Console.WriteLine("Set Xkn : {0}", set.Length);

    var classes = Group.AllOrbits(g, set, ImgPermArr);
    Console.WriteLine("{0, -20} is {1}-transitive of degree {2} : {3}", $"{g} of order {g.Count()}", k, n, classes.Count == 1);
    Group.DisplayOrbx(classes);
    Console.WriteLine();
}
