using System.Security.Cryptography;
using System.Text;
using Craft;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.Tools;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Group.ActivedStorage(false);

Perm.Style = DisplayPerm.CyclesComma;

Perm ConcatPerm(params Perm[] perms)
{
    var Ns = perms.Select(p => p.Sn.N).ToArray();
    var sn = new Sn(Ns.Sum());
    var cum = Ns.Aggregate(new[] { 0 }, (acc, e) => acc.Append(acc.Last() + e).ToArray());
    var table = perms.Zip(cum).SelectMany(e => e.First.Table.Select(i => i + e.Second)).ToArray();
    return sn.CreateElementTable(table);
}

Perm PaddingRight(Perm perm, int pad) => pad == 0 ? perm : ConcatPerm(perm, (new Sn(pad)).Neutral());
Perm PaddingLeft(Perm perm, int pad) => pad == 0 ? perm : ConcatPerm((new Sn(pad)).Neutral(), perm);

Perm[] CyclesSplit(int m)
{
    var dec = IntExt.PrimesDec(m).Select(e => e.Key.Pow(e.Value)).Order().ToArray();
    return dec.Select(e => new Sn(e).Cycle(e.Range(1))).ToArray();
}

Perm Cycles(int m)
{
    if (m == 1)
        return new Sn(1).Neutral();

    return ConcatPerm(CyclesSplit(m));
}

IEnumerable<int[]> CycleWalk(int[] c)
{
    var n = c.Length;
    return n.SeqLazy(1).Select(i => n.SeqLazy(i).Select(j => c[j % n]).ToArray());
}

IEnumerable<int[][]> BlocksPermutes(int[][] cycles)
{
    return cycles.GroupBy(c => c.Length).OrderBy(e => e.Key)
        .Select(a => (ord: a.Key, block: a.ToArray()))
        .Select(a => new Sn(a.block.Length).Select(p => (a.ord, block: p.Apply(a.block))))
        .ToArray()
        .MultiLoop()
        .Select(l => l.SelectMany(m => m.block).ToArray());
}

IEnumerable<Perm> InnerAut(Perm a, Perm b)
{
    var sn = a.Sn;
    if (sn.N != b.Sn.N || !a.PermType.SequenceEqual(b.PermType))
        throw new($"a:{sn} a:{a.PermTypeStr} b:{b.Sn} b:{b.PermTypeStr}");

    var aOrbits = a.Orbits.OrderBy(e => e.Length).ToArray();
    var bOrbits = b.Orbits.OrderBy(e => e.Length).ToArray();
    foreach (var blocks in BlocksPermutes(bOrbits))
    {
        foreach (var cycle in blocks.Select(c => CycleWalk(c)).MultiLoop().Select(l => l.ToArray()))
        {
            var autTable = new int[sn.N];
            foreach (var (bOrb, aOrb) in cycle.Zip(aOrbits))
                for (int i = 0; i < bOrb.Length; i++)
                    autTable[bOrb[i]] = aOrb[i];

            yield return sn.CreateElementTable(autTable);
        }
    }
}

IEnumerable<Perm> InnerAutCn(Perm a)
{
    return IntExt.Coprimes(a.Order).Select(i => a ^ i).SelectMany(ar => InnerAut(a, ar));
}

void TestInnerAutCn(Perm a)
{
    var sn = a.Sn;
    var allConjs1 = InnerAutCn(a).ToArray();
    var act = Group.ByConjugate(sn);
    var H = Group.Cycle(sn, a);
    var allConjs2 = sn.Where(g => H.Keys.Contains(act(g, a))).ToHashSet();

    Console.WriteLine($"a:{a} conj:{allConjs1.Length}");
    foreach (var x in allConjs1.OrderBy(x => H[act(x, a)]).ThenBy(x => x))
        Console.WriteLine($"    x:{x,-40} x*a*x^-1 = a^{H[act(x, a)],-2} = {act(x, a)}");

    Console.WriteLine($"a:{a} conj:{allConjs2.Count}");
    foreach (var x in allConjs2.OrderBy(x => H[act(x, a)]).ThenBy(x => x))
        Console.WriteLine($"    x:{x,-40} x*a*x^-1 = a^{H[act(x, a)],-2} = {act(x, a)}");

    var check = allConjs2.Count == allConjs1.Length && allConjs2.SetEquals(allConjs1);
    Console.WriteLine($"a:{a} conj:{allConjs1.Length} Set Equal:{check}");
    Console.WriteLine();
    if (!check)
        throw new();
}

void TestInnerAutCnLazy(Perm a)
{
    var sn = a.Sn;
    var allConjs = InnerAutCn(a);
    var act = Group.ByConjugate(sn);
    var H = Group.Cycle(sn, a);

    var nb = 0;
    Console.WriteLine($"a:{a}");
    foreach (var x in allConjs)
    {
        ++nb;
        Console.WriteLine($"    x:{x,-40} x*a*x^-1 = a^{H[act(x, a)],-2} = {act(x, a)}");
    }

    Console.WriteLine($"sn:{sn} a:{a} total:{nb}");
    Console.WriteLine();
}

void TestInnerAut(Perm a, Perm b)
{
    var sn = a.Sn;
    var allConjs1 = InnerAut(a, b).ToArray();
    var act = Group.ByConjugate(sn);
    var allConjs2 = sn.Where(g => act(g, a).Equals(b)).ToHashSet();
    var digits = -(sn.N + 1).SeqLazy().Sum(e => $"{e},  ".Length);
    var fmt = $"    x:{{0,{digits}}} = {{1,{digits}}} x*a*x^-1 = b = {{2,{digits}}} {{3}}";

    Console.WriteLine($"a:{a} b:{b} conj:{allConjs1.Length}");
    foreach (var x in allConjs1.OrderBy(x => x.Order).ThenBy(x => x.Orbits.Length).ThenBy(x => x))
        Console.WriteLine(fmt, x, $"[{x.Table.Glue(", ")}]", act(x, a), act(x, a).Equals(b));

    Console.WriteLine($"a:{a} b:{b} conj:{allConjs2.Count}");
    foreach (var x in allConjs2.OrderBy(x => x.Order).ThenBy(x => x.Orbits.Length).ThenBy(x => x))
        Console.WriteLine(fmt, x, $"[{x.Table.Glue(", ")}]", act(x, a), act(x, a).Equals(b));

    var check = allConjs2.Count == allConjs1.Length && allConjs2.SetEquals(allConjs1);
    Console.WriteLine($"a:{a} b:{b} conj:{allConjs1.Length} Set Equal:{check}");
    Console.WriteLine();
    if (!check)
        throw new();
}

(HashSet<Perm> elements, List<Perm> uniqueGenerators) UniqueGenerators(ConcreteGroup<Perm> g, Perm[] generators)
{
    HashSet<Perm> tmpElements = new() { g.Neutral() };
    List<Perm> uniqueGenerators = new();
    var fp = generators.GroupBy(e => e.Orbits.Select(f => f.ToXSet()).ToXSet())
        .ToDictionary(
            e => e.Key,
            e => e.OrderByDescending(r => r.Order).ThenBy(r => r.Orbits.Length).ToArray());
    // fp.Remove(g.Neutral().Orbits.Select(f => f.ToXSet()).ToXSet());
    var fpG = new Dictionary<XSet<XSet<int>>, Perm>();
    while (g.Count() != tmpElements.Count)
    {
        var xset = fpG.Keys.Aggregate(new XSet<XSet<int>>(), (acc, be) => acc.X.Concat(be.X).ToXSet());
        var fp1 = fp.OrderBy(e => e.Key.Grid2D(xset).Select(f => f.t1.Intersect(f.t2).ToXSet()).ToXSet())
            .ThenByDescending(e => e.Value[0].Order)
            .ThenBy(e => e.Value[0].Orbits.Length)
            .ToArray();
        var elt = fp1.Length == 0 ? generators.First() : fp1.First().Value.Last();

        uniqueGenerators.Add(elt);
        var xelt = elt.Orbits.Select(f => f.ToXSet()).ToXSet();
        fpG[xelt] = elt;
        fp.Remove(xelt);
        tmpElements = Group.GenerateElements(g, tmpElements, uniqueGenerators);
        foreach (var xelt0 in tmpElements.Select(elt0 => elt0.Orbits.Select(f => f.ToXSet()).ToXSet()))
            fp.Remove(xelt0);
    }

    return (tmpElements, uniqueGenerators);
}

HashSet<T> GenerateElementsLimited<T>(IGroup<T> bg, T[] generators, int limits) where T : struct, IElt<T>
{
    var elements = new HashSet<T>([bg.Neutral()]);
    var q = new Queue<T>(elements);
    HashSet<T> generatedElements = new HashSet<T>(elements);
    while (q.Count != 0)
    {
        var e1 = q.Dequeue();
        foreach (var e2 in generators)
        {
            var e3 = bg.Op(e2, e1);
            if (generatedElements.Add(e3))
                q.Enqueue(e3);

            if (limits < generatedElements.Count)
                return [];
        }
    }

    return generatedElements;
}

var (success, total, maxTests) = (0, 0, 0);

void SearchAutPerm(ConcreteGroup<Perm> G)
{
    ++total;
    var sn = G.Neutral().Sn;
    var Gset = G.ToSet();
    DisplayGroup.HeadElements(G);
    DisplayGroup.Generators(G);

    var gens0 = G.GroupBy(e => e.Orbits.Select(f => f.ToXSet()).ToXSet()).Select(e => e.First()).ToArray();
    var gensG = UniqueGenerators(G, gens0.OrderByDescending(r => r.Order).ThenBy(r => r.Orbits.Length).ToArray())
        .uniqueGenerators.OrderByDescending(e => e.Order).ToArray();
    gensG.Println(l => $"{l.Order,-4} {l}");

    var actSet = Group.ByConjugateSet(sn);
    var autG = Group.AutomorphismGroup(G);
    DisplayGroup.HeadOrders(autG, false);

    var gensAut = Group.UniqueGenerators(autG, autG.OrderByDescending(e => autG.ElementsOrders[e]).ToArray()).uniqueGenerators;
    Console.WriteLine($"gensAut:{gensAut.Count} in {sn}");

    var autGcount = autG.Count();
    var autEltsOrd = autG.ElementsOrdersList().Order().ToArray();

    var gensImg = gensG.Index().ToDictionary(
        e => e.Index,
        e => autG.Where(f => Perm.TypeEquals(f[e.Item], e.Item))
            .Select(f => f[e.Item]).Distinct()
            .OrderBy(f => Perm.Distance(e.Item, f))
            .ToArray()
    );
    var map = gensAut.Count.SeqLazy()
        .Select(i => (a: gensG[i % gensG.Length], b: gensImg[i % gensG.Length][i / gensG.Length]))
        .ToArray();
    Console.WriteLine($"map:{map.OrderByDescending(e => e.a.Order).Glue(", ")}");

    var gensInnerAut = map.OrderByDescending(e => e.a.Order)
        .Select(e => InnerAut(e.a, e.b).Where(g => actSet(g, Gset).Equals(Gset)))
        .ToArray()
        .MultiLoop()
        .Select(l => l.ToArray());

    var nbTests = 0;
    foreach (var gens1 in gensInnerAut)
    {
        ++nbTests;
        var H1 = GenerateElementsLimited(sn, gens1, autGcount);
        if (H1.Count == autGcount && autEltsOrd.SequenceEqual(H1.Select(e => e.Order).Order()))
        {
            var H2 = Group.Generate(autG.Name, sn, gens1);
            gens1.Println(e => $"{e.Order,-4} {e}", $"check #{nbTests} {H2.ShortName}");
            if (H2.IsIsomorphicTo(autG))
            {
                DisplayGroup.HeadOrdersGenerators(H2);
                Console.WriteLine($"nbTests:{nbTests}");
                Console.WriteLine();
                ++success;
                maxTests = int.Max(maxTests, nbTests);
                return;
            }
        }
    }

    Console.Beep();
    Console.WriteLine($"candidates:{nbTests}");
    Console.WriteLine($"######## NOT FOUND Aut({G.Name}) for {G.ShortName}");
    Console.WriteLine();
    maxTests = int.Max(maxTests, nbTests);
}

(ConcreteGroup<Perm> autPg, Dictionary<Perm, int> idx) AutRegPerm(ConcreteGroup<Automorphism<Perm>> aut)
{
    var Dom = aut.Neutral().Domain;
    var sn = Dom.Neutral().Sn;
    if (sn.N != Dom.Count())
        throw new();

    var idx = Dom.Index().ToDictionary(e => e.Item, e => e.Index);
    var gens = aut.GetGenerators().Select(e => e.AutMap.ToDictionary(f => idx[f.Key], f => idx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .ToArray();

    return (Group.Generate(aut.Name, sn, gens), idx);
}

void AutPermDicyclic()
{
    GlobalStopWatch.Restart();
    for (int ord = 3; ord < 10; ord++)
    {
        var (name, sn, gens) = GroupPermutationForm.DicyclicGens(ord);
        var G = Group.Generate(name, sn, gens.ToArray());
        SearchAutPerm(G);
    }

    GlobalStopWatch.Show($"Missing:{total - success}/{total} Max:{maxTests}");
}

void AutPermMetacyclic()
{
    for (int ord = 4; ord < 32; ord++)
    {
        foreach (var (m, n, r) in GroupPermutationForm.MetaCyclicSdp(ord))
        {
            if ((m, n, r) == (9, 3, 4))
                continue;
    
            var (name, sn, gens) = GroupPermutationForm.MetaCyclicGens(m, n, r);
            var G = Group.Generate(name, sn, gens.ToArray());
            SearchAutPerm(G);
        }
    }

    GlobalStopWatch.Show($"Missing:{total - success}/{total} Max:{maxTests}");
}

void FirstExamples()
{
    TestInnerAutCn(Cycles(4));
    TestInnerAutCn(Cycles(5));
    TestInnerAutCn(Cycles(6));
    TestInnerAutCn(Cycles(10));
    TestInnerAutCn(ConcatPerm(Cycles(3), Cycles(4)));

    {
        var a = ConcatPerm(Cycles(3), Cycles(2));
        var b = ConcatPerm(Cycles(2), Cycles(3));
        TestInnerAut(a, b);
    }
    {
        var a = ConcatPerm(Cycles(3), Cycles(3));
        var b = ConcatPerm(Cycles(3), Cycles(3));
        TestInnerAut(a, b);
    }
    {
        var a = ConcatPerm(Cycles(3), Cycles(5));
        var b = ConcatPerm(Cycles(5), Cycles(3));
        TestInnerAut(a, b);
    }
}

void TestAutRegularPerm<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var Gpg = G.ToPermGroup().Item1;
    DisplayGroup.HeadOrdersGenerators(Gpg);
    var autG = AutRegPerm(Group.AutomorphismGroup(Gpg)).autPg;
    DisplayGroup.HeadOrdersGenerators(autG);
}

{
    TestAutRegularPerm(FG.MetaCyclicSdpMat(4, 4, 3));
    TestAutRegularPerm(FG.DiCyclic(7));
    
    // FirstExamples();
    // AutPermDicyclic();
    // AutPermMetacyclic();
}