using FastGoat.Structures;

namespace FastGoat.UserGroup.GModuleN;

public record SysSolution<Tn, Tg>(int total, List<SysReduction> sreds, CrMap<Tn, Tg> allMaps)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>;