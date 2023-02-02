using FastGoat.Commons;

namespace FastGoat.Structures.GenericGroup;

public readonly struct ConjugacyClasses<T> where T : struct, IElt<T>
{
    public ConcreteGroup<T> Gr { get; }

    public ConjugacyClasses(ConcreteGroup<T> gr)
    {
        Gr = gr;
        Classes = Group.AllConjugacyClassesNames(gr);
        ClassName = new();
        Repr2Idx = new();
        Idx2Repr = new();
        EltRepr = new();
        ClassOrbx = new();
        ClassStabx = new();
        AllReprs = new();
        for (int i = 0; i < Classes.Length; ++i)
        {
            var (name, repr, stabx, orbx) = Classes[i];
            AllReprs.Add(repr);
            ClassName[repr] = name;
            Repr2Idx[repr] = i;
            Idx2Repr[i] = repr;
            ClassOrbx[repr] = orbx;
            ClassStabx[repr] = stabx;
            foreach (var e in orbx)
                EltRepr[e] = repr;
        }
    }

    public void Display()
    {
        var digits = Classes.Max(e => $"{e.repr}".Length);
        var fmt = $"{{0,-{digits}}}";
        foreach (var e in Classes)
        {
            Console.WriteLine(
                $"{e.name,-3} = {string.Format(fmt, e.repr)} {$"Stab({e.name})",-10}:{e.stabx.Count,-4} {$"Orb({e.name})",-10}:{e.orbx.Count,-4}  {e.orbx.Glue(", ")}");
        }
        
        Console.WriteLine($"Nb Classes:{Classes.Length}");
        Console.WriteLine();
    }
    
    private (string name,T repr, HashSet<T> stabx, HashSet<T> orbx)[] Classes { get; }

    private Dictionary<T, string> ClassName { get; }
    private Dictionary<T, int> Repr2Idx { get; }
    private Dictionary<int, T> Idx2Repr { get; }
    private Dictionary<T, T> EltRepr { get; }
    private Dictionary<T, HashSet<T>> ClassOrbx { get; }
    private Dictionary<T, HashSet<T>> ClassStabx { get; }
    private HashSet<T> AllReprs { get; }

    public string GetClassName(T e) => ClassName.ContainsKey(e) ? ClassName[e] : ClassName[EltRepr[e]];
    public string GetClassName(int i) => ClassName.ContainsKey(Idx2Repr[i]) ? ClassName[Idx2Repr[i]] : ClassName[EltRepr[Idx2Repr[i]]];
    public int GetIndex(T e) => Repr2Idx.ContainsKey(e) ? Repr2Idx[e] : Repr2Idx[EltRepr[e]];
    public T GetRepresentative(int i) => Idx2Repr[i];
    public T GetRepresentative(T e) => EltRepr[e];
    public int GetClassOrder(T e) => ClassOrbx[e].Count;
    public int GetClassCentralizerOrder(T e) => ClassStabx[e].Count;
    public IEnumerable<T> GetClassOrbx(T e) => ClassOrbx[e];
    public IEnumerable<T> GetClassOrbx(int i) => ClassOrbx[Idx2Repr[i]];
    public IEnumerable<T> GetClassStabx(T e) => ClassStabx[e];
    public IEnumerable<T> GetClassStabx(int i) => ClassStabx[Idx2Repr[i]];
    public IEnumerable<T> GetRepresentatives() => AllReprs;

    public bool ValidClasses(IEnumerable<T> classes) => AllReprs.SetEquals(classes);
}