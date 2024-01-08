using System.Text.RegularExpressions;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.UserGroup.DatabaseSmallGroups;

public readonly struct IdGroup : IEquatable<IdGroup>, IComparable<IdGroup>
{
    public int Order { get; }
    public int No { get; }

    public string Name { get; }

    public string Txt { get; }
    public SubGroupsInfos Infos { get; }

    private const string reg = @"\((\d*),(\d*)\);(.*);\((\d*),(\d*),(\d*)\)";

    public IdGroup(string txt)
    {
        Txt = txt;

        Match m = Regex.Matches(txt, reg)[0];
        var (idG, noG, nameG) = (m.Groups[1], m.Groups[2], m.Groups[3]);
        var (allSubGrG, allConjsClG, allNormsG) = (m.Groups[4], m.Groups[5], m.Groups[6]);

        (Order, No) = (int.Parse(idG.Value), int.Parse(noG.Value));
        Name = nameG.Value;
        Infos = new SubGroupsInfos(int.Parse(allSubGrG.Value), int.Parse(allConjsClG.Value), int.Parse(allNormsG.Value));
    }

    public string FullName => $"SmallGroup({Order},{No}) Name:{Name}";

    public void Show()
    {
        Console.WriteLine($"Txt:{Txt}");
        Console.WriteLine(FullName);
        Console.WriteLine(Infos);
        Console.WriteLine();
    }

    public bool Equals(IdGroup other) => string.Equals(Txt, other.Txt);

    public int CompareTo(IdGroup other) => (Id: Order, No).CompareTo((other.Order, other.No));

    public override string ToString() => Txt;
}